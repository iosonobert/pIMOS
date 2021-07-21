# pIMOS

UWA Ocean Dynamics focussed python coding of IMOS tools. Click here for the [IMOS Wiki](https://github.com/aodn/imos-toolbox/wiki). 

The system consists of:

1. A Deployment Database
   1. IMOS database is in MS Access
   2. The pIMOS will be in postgreSQL. It will be inspired by the IMOS database, but no requirement for a matching API/structure.  
2. Instrument specific NetCDF templates that are 'as CF Compliant as possible'.
   1. Will only code instruments used by UWA Ocean Dynamics, not coding up all of the IMOS instruments [yet]
3. Parsers for each instrument file into compliant
   1. pIMOS contains a series of classes and superclasses which wrap xarray datasets and handle all of the conventions and metadata tracking
4. Pre-processing 
   1. IMOS call things like bin mapping and coordinate transforms ['pre-processing'](https://github.com/aodn/imos-toolbox/wiki/PPRoutines#coordinate-transformation-from-enu-to-beam-of-nortek-adcp-velocity-data---adcpnortekvelocityenu2beampp)
   2. This will all exist in pIMOS [so far as it is relevant to UWA instruments]. 
      1. The xarray Dataset wrapper will handle tracking of any pre-processing in the `Comments` global/variable attributes of the netCDF dataset.
      2. Ultimate goal to have the terminology and method nams as close to the IMOS terminology as possible.
5. Quality Control
   1. IMOS outlines many [QC procedures](https://github.com/aodn/imos-toolbox/wiki/QCProcedures)
   2. pIMOS will implement all that are relevant to UWA instruments. 
      1. The xarray Dataset wrapper will handle tracking of any QC in the `Comments` global/variable attributes of the netCDF dataset.
      2. Ultimate goal to have the terminology and method nams as close to the IMOS terminology as possible.

**Initial focus will be on moored instruments. Extend to profiling [Solander CTD and VMP] when time allows.

## Database
- Hold metadata for:
    - Experiments
    - Field trips
    - Sites
    - Moorings
    - Deployment files [including URL]
    - UWA raw netcdf files [including URL]
    - UWA processed netcdf files [including URL]
    - LISST ** Not done yet
- Written in PostgreSQL [Code first database already written]
- Script to populate code first database from Matt Rayson's SQLLite database written also
- Still need to write tables for raw and processed netcdf data.

[Here](notebooks/postgreSQL_db_create.ipynb) is a draft notebook for generating a postgres database and adding some data. Ignore all the django stuff in that folder - that's another rabbit hole we don't need right now. 

## Python code
Code falls into 5 categories:
1. Raw file readers
1. The actual zutils.xrwrap object [common for ALL* instruments]
    - This class controls all the QAQC codes CF conventions etc. 
    - Each instrument will naturally inherit from this class, adding unique features if necessary
1. Code to parse from various data formats into the common data object
1. 'Pre-processing' code
1. QC code
1. External processing libraries

### 1. Raw file readers
Where possible, these should read from binary files and avoid complex dependencies or other human interference. Dependancies should be reduced where possible. If code is borrowed, acknowledge it. Dependence on Dolfyn for exmple is not ideal as Kilcher is always changing his API, but the author must be acknowledged.

- SBE56 [Done using UHDAS code]
- SBE39 [Done using UHDAS code]
- SBE37 [Done using UHDAS code]
- RDI ADCP [Done with Dolfyn]
- Vector [Done with Dolfyn and Dall's Porpoise]
  - Dall's Porpoise is itself just a rip off of Dolfyn but handles some corrupt VEC files
- Signature [Done with Dolfyn]
- WetlabsNTU [reads .raw files with pandas]

** Profiling instruments [CTD/VMP] are only partially done, and will be the lowest priority. These will be added once the rest of the system is complete

### 2/3. Parse into xrwrap object
[Here](xarray_wrapper_intro.ipynb) is a basic intro to the xarray wrapper and how to interface with the underlying NetCDF file.

Instrument specific examples are covered in the following notebooks:

- [SBE56 T, SBE39 TP, SBE39 T, SBE37 CTD, SBE37 CT](./notebooks/xarray_wrapper_Seabird_37_39_and_56.ipynb)
- [RDI ADCP](./notebooks/xarray_wrapper_RDI_LONGRANGER_ADCP.ipynb)
- [Nortek Vector](./notebooks/xarray_wrapper_Nortek_Vector.ipynb)
- [Nortek Signature](./notebooks/Nortek_Signature.ipynb)
- [WetLABS NTU](./notebooks/xarray_wrapper_WetlabsNTU.ipynb)
- LISST

Only the Seabird example runs through every file [for KISSME] on a loop. The others are just individual examples. 

### 4/5. Preprocessing/QC code
To the extent that this exists it is all in those notebooks above. We have:
- Despiking [Goring Nikora]
- Phase unwrapping
- Bin Mapping
- ADV and ADCP coordinate transforms
- Acoustic instrument QAQC [incomplete]

Down the line we should add:

- Ocean mooring code for vertically concatenating a thermistor string of multiple CTDs

- Motion correction
- More auto QC routines with IMOS style names

** NOTE. A point where we may differ from IMOS. We will often 'replace' data removed in QC, rather than just flagging it. Need to think how this is handled. Maybe we should have a

1. 'modified' flag, or;
2. Just add the comment to the nc that the change was made, but not actually flag exactly which variables were changed. 

## 6. External processing code

Other repos exist for actual calculations, and these will be appropriately linked to this repo with examples.
- mrayson/sfoda
    - Signal processing/filtering/harmonic analysis
    - Modal decomposition
- mrayson/mycurrents
    - stacking ocean mooring
    - xarray visualisation
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

# Process Levels

Another thing IMOS defines are [Process(ing) Levels](https://github.com/aodn/imos-toolbox/wiki/ProcessLevelsAndQC). In a nutshell, this is how much processing has been performed, and the levels are:

- FV00 Raw Data
- FV01 Quality Controlled Data
- FV02 Derived Products
- FV03 Interpreted Products
- FV04 Knowledge Products

We must decide:

1. Whether we agree with these designations for our purposes [e.g. I'm not sure how you can call something QCd until there has been some attempt to interpret it]
2. How far I'm personally intended to take this archiving
3. What to do, going forward, with everyone's Derived Products.

# Dependencies

Current state of play:

- netcdf4

- xarray
- Dolfyn [for now]
- Dall's Porpoise [needed for reading corrupt .VEC files, otherwise can just use Dolfyn for everything] 
- zutils [for convenience for now, anything here can be moved to pIMOS.utils]
- psycopg2 [for connecting to postgres database]
- windrose [this is not a good dependency, amend this]
- jupyter [not really a dependency but the examples use jupyter to help render in github]


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

