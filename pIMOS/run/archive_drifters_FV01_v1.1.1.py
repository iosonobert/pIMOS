import os
import argparse
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pIMOS.xrwrap.drifter as drifter
import pIMOS.utils.UWA_archive_utils as ai 


def main(output_dir, exp):
    """
    Process drifter data from NAHYEAH and generate outputs and visualizations.
    
    Parameters
    ----------
    output_dir : str
        Directory to save output files
    exp : str
        Experiment name
    """


    csv_name = 'UWA_Carthe_Drifters_Rottnest_latest.csv'
    url = 'https://nahyeah.mywire.org/files/' + csv_name
    df_raw = pd.read_csv(url)
    df_raw.to_csv(os.path.join(output_dir, csv_name), index=False)

    # Save new CSV based on unique entries in DeviceName and DeviceDateTime
    unikeys = df_raw['DeviceName'].drop_duplicates()
    drifter_files = []
    for uk in unikeys:
        df_uni = df_raw[df_raw['DeviceName']==uk]
        df_uni = df_uni[~df_uni.index.duplicated(keep='first')]
        df_uni = df_uni.sort_index()
        df_uni.to_csv(os.path.join(output_dir, uk + '.csv'))
        drifter_files.append(os.path.join(output_dir, uk + '.csv'))
        
    # Read text file to pandas
    nms = ['Latitude','Longitude','GpsQuality','GPS Confidence']
    drogue_depth = 0.3


    for ff in drifter_files:
        
        rr, ds = drifter.from_csv(ff,
                                timevar='DeviceDateTime',
                                keepcols=nms,
                                parse_dates=True,
                                index_col=5)
                                # encoding='UTF-16 LE')

        rr.align_names({'Latitude':'latitude', 'Longitude':'longitude', 'GPS Confidence':'GPSConfidence'})

        rr.update_attributes_with_dict({'nominal_latitude':np.nan, 'nominal_longitude':np.nan,\
                                        'nominal_site_depth':drogue_depth, 'project':exp,\
                                        'instrument_serial_number':os.path.split(ff)[1][:-4],\
                                        'is_profile_data':1})

        rr.drop_pIMOS_coords() # These are reassigned with every call to update_attributes_with_dict

        # Set lat,lon = (0,0) to nan
        rr.ds['latitude'] = rr.ds['latitude'].where(rr.ds['latitude']!=0)
        rr.ds['longitude'] = rr.ds['longitude'].where(rr.ds['longitude']!=0)
        rr.calc_easting_northing('longitude', 'latitude')
        # rr.despike_drifter(window=3, spike_dist=400, reinstate=0.2,\
        #                     timecut=np.timedelta64(1, 'm'), qc_var=rr.ds['GPSConfidence']==0)
        rr.despike_drifter(window=3, spike_dist=100, reinstate=0.1,\
                            timecut=None, qc_var=rr.ds['GPSConfidence']==0)
        
        rr.assign_nickname('Drifter') # Needed to export
        ai.pIMOS_export(rr, output_dir)

        # Print summary figures
        ats = ai.nonempty_attrs(rr)
        ats.pop('institution')
        rr.calc_distance()
        rr.calc_raw_speed()
        rr.associate_qc_flag('raw_speed', 'position')
        fig = ai.plot_temp(rr, None, None, ats, variable='raw_speed')
        png_path = rr.fullpath_last_export[0:-3] +'_raw_speed.png'
        fig.savefig(png_path, dpi=300)
        plt.close(fig)

        fig = ai.plot_temp(rr, None, None, ats, variable='north_vel')
        png_path = rr.fullpath_last_export[0:-3] +'_processed_north_vel.png'
        fig.savefig(png_path, dpi=300)
        plt.close(fig)
        
        fig = ai.plot_temp(rr, None, None, ats, variable='east_vel')
        png_path = rr.fullpath_last_export[0:-3] +'_processed_east_vel.png'
        fig.savefig(png_path, dpi=300)
        plt.close(fig)

        # Print map figure
        fig, ax = rr.plot_drifter(attrs=ats)
        png_path = rr.fullpath_last_export[0:-3] +'_map.png'
        fig.savefig(png_path, dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process NAHYEAH drifter data')
    parser.add_argument('output_dir', help='Directory to save output files')
    parser.add_argument('exp', default='NAHYEAH26', help='Experiment name (e.g., NAHYEAH26)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    main(args.output_dir, args.exp)
