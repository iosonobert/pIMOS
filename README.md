# pIMOS

# Objective:
Pull together:

1. Database for raw and processed data 
2. Python tools for reading, quality controlling,  exporting and all moored instruments we use at UWA. 

Extend to profiling [CTD and VMP] when time allows

## Database
- Hold metadata for:
    - Experiments
    - Field trips
    - Sites
    - Moorings
    - Deployment files [including URL]
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
1. External processing libraries

### 1. Raw file readers
Where possible, these should read from binary files and avoid complex dependencies. Dependence on Dolfyn for exmple is not ideal as Kilcher is always changing his API.  
- SBE56 [Done using UHDAS code]
- SBE39 [Done using UHDAS code]
- SBE37 [Done using UHDAS code]
- RDI ADCP [Done with Dolfyn]
- Vector [Done with Dolfyn and Dall's Porpoise]
  - Dall's Porpoise is itself just a rip off of Dolfyn but handles some corrupt VEC files
- Signature [Done with Dolfyn]
- WetlabsNTU [.raw]


Profiling instruments [CTD/VMP] are only partially done, and will be the lowest priority. These will be added once the rest of the system is complete

### 2/3. Parse into xrwrap object
These currently reside within the crude_readers package, but will ultimately be moved into a UWA repository. 
- SBE56 drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56_[using_xrwrap].ipynb)):
    - Rawfile [UHDAS code]
    - Netcdf 
    - xr
- SBE39 drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56_[using_xrwrap].ipynb)):
    - Rawfile [UHDAS code]
    - Netcdf 
    - xr
- SBE37 drivers ([click here for jupyter notebook](https://github.com/iosonobert/pIMOS/blob/master/notebooks/Seabird_37_39_and_56_[using_xrwrap].ipynb)):
    - Rawfile [UHDAS code]
    - Netcdf 
    - xr
- RDI ADCP drivers ([click here for jupyter notebook](./notebooks/RDI%20ADCP%20%5BSP250%20LR%20using%20xrwrap%5D.ipynb)):
    - Rawfile [Dolfyn]
    - Netcdf 
    - xr
- Vector drivers ([click here for jupyter notebook](notebooks/Vector%20.VEC%20Read%20and%20Clean%20%5BxrWrap%5D.ipynb)):
    - Rawfile [Dolfyn or Dall's Porpoise]
    - Netcdf 
    - xr
- Signature:
    - Not yet in xrwrap
- WetlabsNTU/WetLABSFLNTU/WetlabsNTUSB
- LISST

The code for generating the database is in [this notebook](postgres/db_create.ipynb). Ignore all the django stuff in that folder - that's another rabbit hole we don't need right now. 

Data processing examples are contained in the following notebooks:
- [SBE56 T, SBE39 TP, SBE39 T, SBE37 CTD, SBE37 CT](./notebooks/Seabird_37_39_and_56_[using_xrwrap].ipynb)
- [RDI ADCP](./notebooks/RDI_ADCP_%5BSP250_LR_using_xrwrap%5D.ipynb)
- [Nortek Vector](./notebooks/Vector_.VEC_Read_and_Clean_%5Busing_xrWrap%5D.ipynb)
- [Nortek Signature](./notebooks/Nortek_Signature.ipynb)

### 4. Processing code
Nothing actually done yet, but the hope is to pull together at least the following:
- Ocean mooring code for vertically concatenating a thermistor string of multiple CTDs
- Despiking
- Phase unwrapping
- Acoustic instrument auto QAQC

## 5. External processing code

Other repos exist for actual calculations, and these will be appropriately linked to this repo with examples.
- mrayson/sfoda
    - Signal processing/filtering/harmonic anaysis
    - Modal decomposition
- mrayson/mycurrents
    - stacking ocean mooring
    - data visualisation
- iosonobert/turbo_tools
    - Turbulence calcs
        - IDM dissipation
        - Variance method
        - Structure function method
    - Turbulence QC
        - Phase unwrapping
        - Goring Nikora
- iosonobert/microstructure
- iosonobert/bert_ms_private
- Whitwell git hub?
    - Motion correction?
    - Wavelet stuff?

# Dependencies

Current state of play:

- netcdf4

- xarray
- Dolfyn [for now]
- Dall's Porpoise [needed for reading corrupt .VEC files, otherwise can use Dolfyn] 
- zutils [for convenience for now, anything here can be moved to pIMOS.utils]
- crude_readers [this will need to change]
- psycopg2 [for connecting to postgres database]
- windrose [this is not a good dependency, amend this]
- jupyter [for using examples]

- 



# Known Issues

Making a list here rather than using the issues register in case I need to kill this repo and recreate. 

- Database

  - Incomplete

- Instrument readers

  - RDI ADCP

    - Only really handles BEAM coordinates
    - Doesn't yet save all of the instrument config

  - Signature

    - Only handles BEAM coordinates
    - Doesn't yet save all of the instrument config
    - Not yet on xrwrap

  - Vector

    - Doesn't handle corrupt IMU files. Nobody has written this code. 
    - Only handles XYZ coordinates

  - Seabird 56

    - Time not reading correctly

  - Microstructure

    - Python translation of ODAS needs thorough testing
      - Is there a PhD student that can do this?
      - If not Whitwell and I will pump it out

    

- CF conventions

  - Can always put more time into this

- QC

  - Not all IMOS procedures coded
  - Not following IMOS names
  - 

