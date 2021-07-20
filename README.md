# pIMOS

# Objective:
Database for raw and processed data and python tools for reading, writing and quality controlling all moored instruments we use at UWA. 

## Database
- Hold metadata for:
    - Experiments
    - Field trips
    - Sites
    - Moorings
    - Deployent files [including URL]
    - UWA raw netcdf files [including URL]
    - UWA processed netcdf files [including URL]
- Written in PostgreSQL [Code first database already written]
- Script to populate code first database from Matt Rayson's SQLLite database written also
- Still need to write tables for raw and processed netcdf data.



[Here](https://github.com/iosonobert/pIMOS/blob/master/postgres/db_create.ipynb) is a draft notebook for generating a postgres database and adding some data.

## Python code
Code falls into 3 categories:
1. Raw file readers
1. The actual zutils.xrwrap object [common for ALL* instruments]
    - This class controls all the QAQC codes CF conventions etc. 
    - Each instrument will naturally inherit from this class, adding unique features if necessary
1. Code to parse from various data formats into the common data object
1. Processing scripts

### Raw file readers
Where possible, these should read from binary files and avoid complex dependencies. Dependence on Dolfyn for exmple is not ideal as Kilcher is always changing his API.  
- SBE56 [Done using UHDAS code]
- SBE39 [Done using UHDAS code]
- SBE37 [Done using UHDAS code]
- RDI ADCP [Done with Dolfyn]
- Vector [Done with Dolfyn and Dall's Porpoise]
  - Dall's Porpoise is itself just a rip off of Dolfyn but handles some corrupt VEC files
- Signature [Done with Dolfyn]


Profiling instruments [CTD/VMP] are only partially done, and will be the lowest priority. These will be added once the rest of the system is complete

### Parse into xrwrap object
These currently reside within the crude_readers package, but will ultimately be moved into a UWA repository. 
- SBE56 drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56.ipynb)):
    - Rawfile [UHDAS code]
    - Netcdf 
    - xr
- SBE39 drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56.ipynb)):
    - Rawfile [UHDAS code]
    - Netcdf 
    - xr
- SBE37 drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56.ipynb)):
    - Rawfile [UHDAS code]
    - Netcdf 
    - xr
- RDI ADCP drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56.ipynb)):
    - Rawfile [Dolfyn]
    - Netcdf 
    - xr
- Vector drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56.ipynb)):
    - Rawfile [Dolfyn or Dall's Porpoise]
    - Netcdf 
    - xr
- Signature:
    - Not yet in xrwrap

The code for generating the database is in [this notebook](postgres/db_create.ipynb). Ignore all the django stuff in that folder - that's another rabbit hole we don't need right now. 

Data processing examples are contained in the following notebooks:
- [SBE56 T, SBE39 TP, SBE39 T, SBE37 CTD, SBE37 CT](./notebooks/Seabird_37_39_and_56.ipynb)
- 

### Processing code
Nothing actually done yet, but the hope is to pull together at least the following:
- Ocean mooring code for vertically concatenating a thermistor string of multiple CTDs
- Despiking
- Phase unwrapping
- Acoustic instrument auto QAQC

Other repos exist for actual calculations, and these will be appropriately linked to this repo with examples.
- mrayson/sfoda
    - Signal processing/filtering/harmonic anaysis
    - Modal decomposition
- iosonobert/turbo_tools
    - Turbulence calcs
        - IDM dissipation
        - Variance method
        - Structure function method

# Dependencies
jupyter
netcdf4
Dolfyn [for now]
zutils
crude_readers [this will need to change]
psycopg2 [for connecting to postgres database]
windrose [this is not a good dependency, amend this]
