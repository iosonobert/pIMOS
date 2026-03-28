import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
import tomllib


def _run(cmd, cwd=None, env=None):
    result = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(
            "Command failed:\n"
            + " ".join(str(x) for x in cmd)
            + "\n\nstdout:\n"
            + result.stdout
            + "\n\nstderr:\n"
            + result.stderr
        )
    return result


def _venv_binaries(venv_dir):
    if os.name == "nt":
        py = venv_dir / "Scripts" / "python.exe"
        pip = venv_dir / "Scripts" / "pip.exe"
    else:
        py = venv_dir / "bin" / "python"
        pip = venv_dir / "bin" / "pip"
    return py, pip


def _direct_dependency_names(pyproject_path):
    data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    deps = data["project"]["dependencies"]
    names = []
    for dep in deps:
        match = re.match(r"^[A-Za-z0-9_.-]+", dep)
        if match:
            names.append(match.group(0))
    return names


def _pinned_version(pyproject_path, package_name):
    data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    deps = data["project"]["dependencies"]
    needle = f"{package_name}=="
    for dep in deps:
        if dep.startswith(needle):
            return dep.split("==", 1)[1]
    return None


@pytest.mark.install
def test_fresh_env_pip_install_and_dependency_resolution():
    repo_root = Path(__file__).resolve().parents[1]
    pyproject_path = repo_root / "pyproject.toml"
    dependency_names = _direct_dependency_names(pyproject_path)

    numpy_pin = _pinned_version(pyproject_path, "numpy")
    if numpy_pin == "1.24.2" and sys.version_info >= (3, 12):
        pytest.skip("numpy==1.24.2 is not compatible with Python 3.12 for this smoke install path")

    with tempfile.TemporaryDirectory(prefix="pimos-install-smoke-") as td:
        venv_dir = Path(td) / "venv"

        _run([sys.executable, "-m", "venv", str(venv_dir)])
        py, pip = _venv_binaries(venv_dir)

        _run([str(py), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])

        install_env = os.environ.copy()
        # Make SCM version inference deterministic for ephemeral smoke-test builds.
        install_env["SETUPTOOLS_SCM_PRETEND_VERSION"] = "0.0.0"
        install_env["SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PIMOS"] = "0.0.0"
        _run([str(pip), "install", str(repo_root)], env=install_env)

        # pip check validates installed package requirements are satisfied.
        _run([str(pip), "check"])

        # Confirm direct dependencies declared in pyproject are installed.
        for dep_name in dependency_names:
            _run([str(pip), "show", dep_name])

        # Validate installed distribution metadata and module discoverability
        # without importing pIMOS (which pulls heavy optional runtime imports).
        _run(
            [
                str(py),
                "-c",
                (
                    "import importlib.metadata as m, importlib.util as u; "
                    "print(m.version('pIMOS')); "
                    "assert u.find_spec('pIMOS') is not None"
                ),
            ]
        )
