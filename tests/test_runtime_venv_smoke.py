import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


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


@pytest.mark.runtime
def test_runtime_loader_smoke_in_fresh_venv():
    repo_root = Path(__file__).resolve().parents[1]

    with tempfile.TemporaryDirectory(prefix="pimos-runtime-smoke-") as td:
        venv_dir = Path(td) / "venv"

        _run([sys.executable, "-m", "venv", str(venv_dir)])
        py, pip = _venv_binaries(venv_dir)

        _run([str(py), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel", "pytest"])

        install_env = os.environ.copy()
        install_env["SETUPTOOLS_SCM_PRETEND_VERSION"] = "0.0.0"
        install_env["SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PIMOS"] = "0.0.0"
        _run([str(pip), "install", str(repo_root)], env=install_env)

        _run([str(py), "-m", "pytest", "-q", "tests/test_netcdf_class_loader.py"], cwd=repo_root)
