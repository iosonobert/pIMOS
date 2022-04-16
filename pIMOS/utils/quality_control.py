import numpy as np
import datetime
"""

The following is the list of IMOS QC procedures for moored instrments (at the time of writing).

Moored time series
    *imosImpossibleDateQC* (compulsory)
    *imosImpossibleLocationSetQC) (compulsory)
    *imosInOutWaterQC* (compulsory)
    *imosGlobalRangeQC* (compulsory)
    imosRegionalRangeQC (optional)
    *imosImpossibleDepthQC* (compulsory)
    pimosVerticalSpikeQC
    pimosDensityInversionSetQC
    imosRateOfChangeQC (optional)
    imosStationarityQC (optional)
    *imosSalinityFromPTQC* (compulsory)
    *imosSideLobeVelocitySetQC* (compulsory)
    *imosTiltVelocitySetQC* (compulsory)
    *imosHorizontalVelocitySetQC* (compulsory)
    *imosVerticalVelocityQC* (compulsory)
    imosSurfaceDetectionByDepthSetQC (optional)
    imosEchoIntensitySetQC (optional)
    imosEchoIntensityVelocitySetQC (optional)
    *imosCorrMagVelocitySetQC* (compulsory)
    imosPercentGoodVelocitySetQC (optional)
    imosErrorVelocitySetQC (optional)
    imosTier2ProfileVelocitySetQC (optional)
    *imosHistoricalManualSetQC* (compulsory)
    imosTimeSeriesSpikeQC (optional)

Placeholders are created here for each of these, but they are not necessarily implemented. 

For ADCP Data IMOS Extensively references Crout et al. (2005) "REAL-TIME OIL PLATFORM OCEAN CURRENT DATA IN THE GULF OF MEXICO:
AN IOOS INDUSTRY PARTNERSHIP SUCCESS STORY".

"""

def pimosImpossibleDateQC():
    raise(NotImplementedError)

def pimosImpossibleLocationSetQC():
    raise(NotImplementedError)

def pimosInOutWaterQC(rr, mooring, db_data, year_1=1990, year_n=2200, delete_raw=False, experiment='kissme', recovered='kissme_recovery'):
    """
    This differs somewhat from IMOS. It takes a mooring table with mooring in/out dates and QCs accordingly. 
    """
    
    df = db_data['possible_mooring_dates']
    
    # Hard Code This For Now
    df = df.loc[strcmpi(df['Recovered'].values, recovered)]
    df = df.loc[strcmpi(df['Experiment'].values, experiment)]
    
    df = df.loc[strcmpi(df['Mooring'].values, mooring)]

    if df.shape[0] == 0:
        raise(Exception("Could not find dates for this mooring. {} {} {}". format(experiment, recovered, mooring)))

    start_date = df['StartDate'].values[0]
    end_date = df['EndDate'].values[0]

    # Start
    rr.update_qc_flag('*', 
                  'time',
                  np.datetime64(datetime.datetime(year_1, 1, 1)),
                  start_date,
                  1,
                  'pimosInOutWaterQC start date: {}'.format(start_date), delete_raw=delete_raw)

    #End
    rr.update_qc_flag('*', 
                      'time',
                      end_date,
                      np.datetime64(datetime.datetime(year_n, 1, 1, 1)),
                     1,
                     'pimosInOutWaterQC end date: {}'.format(end_date), delete_raw=delete_raw)
                     
def pimosGlobalRangeQC():
    """
    imosGlobalRangeQC checks that the parameter's values are within a global valid range. This must be strored in a config file or db. 
    """
    raise(NotImplementedError)

def pimosRegionalRangeQC():
    """
    imosRegionalRangeQC checks that the variables values are within a regional range which is site specific. 
    This test is available at one's discretion who is fully aware of the regional ranges applied to his datasets.
    This must be strored in a config file or db. 
    """
    raise(NotImplementedError)

def pimosImpossibleDepthQC():
    """
    imosImpossibleDepthQC checks that the pressure/depth data are given within an acceptable range. 
    Basically, the test checks that the instrument is close to its nominal depth, allowing enough flexibility to cope with possible 
    mismatch between nominal and actual depth when deployed, and knock down events related to strong current events.
    """
    raise(NotImplementedError)

def pimosVerticalSpikeQC():
    """
    imosVerticalSpikeQC is looking for spikes among vertically adjacent triplet of samples. 
    This test is for profile mode only and is available at one's discretion who is fully aware of the algorithm and flags used.
    """ 

def pimosDensityInversionSetQC():
    """
    imosDensityInversionSetQC flags bad any PSAL + TEMP + CNDC + DENS value that shows an inversion in density or an increase/decrease 
    in density that is greater than a threshold set in imosDensityInversionSetQC.txt.
    This test is for profile mode only and is available at one's discretion who is fully aware of the algorithm and flags used.
    """

def pimosRateOfChangeQC():
    """
    imosRateOfChangeQC is looking for important changes in values among temporally adjacent triplet (or pair) of samples. 
    This test is mostly useful for time series data and is available at one's discretion who is fully aware of the algorithm and flags used.
    """

    """
    SEEMS TO OVERLAP WITH THE SPIKE DETECTION ALGORITHM.
    """

    """
    I THINK WE'D IMPLEMENT THIS AS A GORING NIKORA TYPE APPROACH.
    """

    raise(NotImplementedError)

def pimosStationarityQC():
    """
    imosStationarityQC checks that consecutive equal values of a parameter don't exceed a maximum number. 
    This test is available at one's discretion who is fully aware of the algorithm and flags used.
    """

    """
    THIS COULD BE EXTENDED TO A PROPER RUNS TEST.
    """

    raise(NotImplementedError)

def pimosSalinityFromPTQC():
    raise(NotImplementedError)

def pimosSideLobeVelocitySetQC():
    """
    imosSideLobeVelocitySetQC finds for each timestamp what is the HEIGHT_ABOVE_SENSOR from which any data sample is contaminated by 
    side lobes effect (close to the surface or bottom). Impacted data samples are then flagged accordingly. 
    This test can be applied on any ADCP data providing we know the beam angle of the ADCP transducer. It is based on Nortek's documentation.
    """

    """
    Seems to overlap appreciably with imosSurfaceDetectionByDepthSetQC 
    """

    raise(NotImplementedError)

def pimosTiltVelocitySetQC(rr, thresh_1):
    """
    imosTiltVelocitySetQC is a two-level test which flags velocity data depending on the tilt values. 
    Most ADCPs see their performance degraded with increasing tilt.

    For RDI ADCPs, compass measurements are affected and fail to meet specifications when a first tilt threshold is exceeded. 
    When the second threshold is exceeded coordinates transform and bin-mapping operations are also affected.

    For Nortek ADCPs, velocity data accuracy fails to meet specifications when the first tilt threshold is exceeded. 
    When the second is exceeded then velocity data becomes unreliable.

    For each sample, the tilt is computed using the following formula : TILT = acos(sqrt(1 - sin(ROLL)^2 - sin(PITCH)^2))

    2 thresholds are used: Probably_good_data and Bad_data_that_are_potentially_correctable.
    DEFAULTS:
        INSTRUMENT          1ST_THESH       2ND_THESH
        sentinel	        15	            22	       
        monitor	            15	            22	       
        longranger	        15	            50	       
        long ranger	        15	            50	       
        quartermaster	    15	            50	       
        nortek	            20	            30	       
        """
    
    raise(NotImplementedError)

def pimosTiltVelocitySimpleQC(rr, thresh_1=22, method='max', flag_name='qc_velocity'):
    """
    Here we use a much simpler single threshold. The object must have a  "_calc_tilt" method that returns the instrument tilt. 
    """
    ds = rr.ds

    # This function must be defined for the instrument at hand
    tilt = rr._calc_tilt(method='max')

    logical_index = tilt > thresh_1
    
    rr.update_qc_flag_logical(flag_name, 
                           'time', 
                           logical_index, 
                           1,
                          'pimosTiltVelocitySimpleQC threshold: {}'.format(thresh_1))

def pimosHorizontalVelocitySetQC():
    
    raise(NotImplementedError)

def pimosVerticalVelocityQC():
    raise(NotImplementedError)

def pimosSurfaceDetectionByDepthSetQC():
    """
    Seems to overlap appreciably with pimosSideLobeVelocitySetQC 
    """
    raise(NotImplementedError)

def pimosEchoIntensitySetQC():
    """
    imosEchoIntensitySetQC looks at the difference in echo intensity between adjacent vertical bins. 
    When this difference is greater than a certain threshold, the bin is marked as bad.

    The test is more flexible than the original EchoIntensityVelocitySetQC (see below), in the sense that the test can be executed only from a pre-defined 
    bin-number - bound_by_index (e.g. mark only from bin number=10), be clipped to depth levels - bound_by_depth (mark only above/below =100m), 
    and if we should expand/propagate a bad bin flag to other bins further away from the instrument - propagate flag.

    This test is available at one's discretion who is fully aware of the algorithm and demands expertise to set relevant thresholds for each dataset.

    """
    raise(NotImplementedError)

def pimosEchoIntensityVelocitySetQC():
    """
    imosEchoIntensityVelocitySetQC looks at the difference in echo intensity between adjacent vertical bins. 
    When this difference is greater than a certain threshold for at least 3 beams, this means an obstacle (surface or bottom) is detected and from there, 
    any further bin is flagged bad. 
    
    This test is also available at ones discretion who is fully aware of the algorithm and demands expertise to set relevant thresholds for each dataset.
    """
    raise(NotImplementedError)

def pimosEchoIntensitySimpleQC(rr, thresh_1=20, echo_name='echo', beam_name='beam', flag_name='qc_velocity'):
    """
    Minimum echo amplitude test. This is not in the IMOS toolbox. 
    """

    ds = rr.ds
    
    echo = ds[echo_name]

    beam_index = get_dim_index(echo, beam_name)
    echo = echo.min(axis=beam_index)

    logical_index = echo < thresh_1
    
    rr.update_qc_flag_logical(flag_name, 
                           'time', 
                           logical_index, 
                           1,
                          'pimosEchoIntensitySimpleQC threshold: {}'.format(thresh_1))

def pimosFishDetectionQC(rr, thresh_1=20, echo_name='echo', beam_name='beam', flag_name='qc_velocity'):
    """
    ** Not exactly in the IMOS toolbox

    Run a fish detectrion algorithm that compares echo amplitude across beams. 
    """

    ds = rr.ds
    echo = ds[echo_name]

    beam_index = get_dim_index(echo, beam_name)

    emax = np.max(echo.values, axis=beam_index)
    emin = np.min(echo.values, axis=beam_index)

    logical_index = emax-emin>40

    rr.update_qc_flag_logical(flag_name, 
                           'time', 
                           logical_index, 
                           1, 
                           comment='pimosFishDetectionQC threshold: {}'.format(thresh_1))
                           

    # raise(NotImplementedError)

def pimosFishDetectionBeamwiseQC(rr, thresh_1=20, echo_name='echo', axis=2):
    """
    ** Not in the IMOS toolbox

    Run a fish detectrion algorithm that compares echo amplitude across beams. 
    """
    raise(NotImplementedError)

def pimosFishDetectionVerticalQC(rr, thresh_1=20, echo_name='echo', axis=2):
    """
    ** Not in the IMOS toolbox

    Run a fish detectrion algorithm that compares echo amplitude along beams. 
    """
    raise(NotImplementedError)

def pimosCorrMagVelocitySetQC(rr, thresh_1=60, corr_name='corr', beam_name='beam', flag_name='qc_velocity'):
    """
    imosCorrMagVelocitySetQC test checks that there is sufficient signal to noise ratio to obtain good quality data via the measure of a 
    pulse-to-pulse correlation in a ping. At each bin, if at least 2 beams see their correlation value greater than a threshold value then 
    the sample at that bin passes the test.

    IMOS States: "This test raises a warning if applied without CMAGn being vertically bin-mapped." 
    This warning is not implemented - I'm not sure why it needs to coma after bin mapping. Also ADVs aren't bin mapped.
    
    """
    
    ds = rr.ds
    
    corr = ds[corr_name]
    beam_index = get_dim_index(corr, beam_name)
    
    corr = ds[corr_name].min(axis=beam_index)

    logical_index = corr < thresh_1
    
    rr.update_qc_flag_logical(flag_name, 
                           'time', 
                           logical_index, 
                           1,
                          'pimosCorrMagVelocitySetQC threshold: {}'.format(thresh_1))

def pimosPercentGoodVelocitySetRDIQC(rr, percent_good_name='percent_good', thresh_1 = 0, flag_name='qc_velocity'):
    """
    imosPercentGoodVelocitySetQC test checks that the percentage of good 3 beams and 4 beams solutions is greater than a certain threshold. 
    A percentage of bad solution is based on low correlation and fish detection. 
    This test is available at ones discretion who is fully aware of the algorithm and demands expertise to set relevant thresholds for each dataset.
    
    PERG1 + PERG4 > pgood (performed on the 1st and 4th value only: represents the percentage good of 3 and 4 beam solutions in earth coordinates configuration)
    
    With pgood default value of 50% suggested Janet Sprintal for the RDI Workhorse Quatermaster ADCP 300kHz. However relevant threshold value needs to be set for each dataset! It is stored in imosPercentGoodVelocitySetQC.txt.

    This test is part of a procedure implemented by the U.S. National Data Buoy Center, MI, and detailed in Crout et al. (2005). Thanks to Rebecca Cowley (CSIRO Hobart) this procedure has been corrected from a mis-interpretation, putting more faith into the original document from RDI Darryl R Symonds (2006).

    This same methodology has been applied for the Quality Control of the mooring data collected during the INSTANT (International Nusantara STratification ANd Transport program) deployments in the outflow straits - Lombok Strait and Timor Passage. In addition, the determination of the parameters used in this methodology has been described by Janet Sprintal (2007) in this document (see Appendix 2). Its result specified the default parameters values used in the QC procedure.
    """

    """
    I have changed this name to be RDI specific.
    """

    percent_good = rr.ds[percent_good_name].values
    percent_good_34 = percent_good[:, :, 0] + percent_good[:, :, 3]

    logical_index = percent_good_34 < thresh_1

    rr.update_qc_flag_logical(flag_name, 
                            'time', 
                            logical_index, 
                            1,
                            'pimosPercentGoodVelocitySetRDIQC threshold: {}'.format(thresh_1))


def pimosErrorVelocitySetQC(rr, errvel_name='errvel', thresh_1=0.25, flag_name='qc_velocity'):
    """
    imosErrorVelocitySetQC test checks that the horizontal error velocity (difference between two independant estimates, basically two pair of beams) 
    is smaller than a certain threshold so that the assumption of hozontal flow homogeneity is reasonable. 
    
    This test is only available in 4 beams solution and shouldn't flag as bad any 3 beam solution. 
    
    This test is available at ones discretion who is fully aware of the algorithm and demands expertise to set relevant thresholds for each dataset.
    """
    
    ds = rr.ds
    errvel = ds[errvel_name]

    logical_index = errvel > thresh_1
    
    rr.update_qc_flag_logical(flag_name, 
                           'time', 
                           logical_index, 
                           1,
                          'pimosErrorVelocitySetQC threshold: {}'.format(thresh_1))

def pimosTier2ProfileVelocitySetQC():
    """
    imosTier2ProfileVelocitySetQC test checks that more than 50% of the bins in the water column (imosEchoIntensityVelocitySetQC test gives where the 
    surface/bottom starts) in a profile have passed all of the 5 other velocity tests otherwise the entire profile fails. 

    This test is available at ones discretion who is fully aware of the algorithm and demands expertise to set relevant thresholds for each dataset.
    """
    raise(NotImplementedError)

def pimosHistoricalManualSetQC():
    """
    Anytime a manual QC is performed, a .mqc file is created/updated with a new information about what has been flagged where and with with which flag and comment. The file is at the same location as the input file and has the same radical name.

    imosHistoricalManualSetQC QC procedure reads any existing .mqc file so that any manual QC is automatically re-applied.
    """
    raise(NotImplementedError)

def pimosTimeSeriesSpikeQC():
    """
    IMOS discuss many many spike detection methods, none of which being Goring Nikora 2002.
    
    https://github.com/aodn/imos-toolbox/wiki/QCProcedures
    """
    raise(NotImplementedError)


def strcmpi(lst, string):
    """
    Copy of the matlab function.
    Returns list.
    """
    rtn = [i.lower() == string.lower() for i in lst]
    return rtn
    pass

def get_dim_index(array, dim_name):
    """
    Find which axis of a DataArray corresponds to a certain dimension name.
    """
        
    coord_names =[x for x in array.coords ]
    dim_index =[i for i in np.arange(0, len(coord_names)) if coord_names[i]==dim_name]    
    
    if len(dim_index) == 1:
        dim_index=dim_index[0]
    else:
        dim_index = None
        
    return dim_index


class NotImplementedError(Exception):
    pass