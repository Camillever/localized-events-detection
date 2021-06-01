""" Module to specify all the seismograms worthy to be use
You can modify the conditions in this module"""

import os

from obspy import UTCDateTime

from locevdet.utils import get_info_from_mseedname

def list_seismograms(folder_path:str):
    """  Give the list of all seismograms from the given folder path.
    Conditions could be given as excluding some stations or periods for some stations
    Args:
        folder_path : directory path of seismograms
    Returns : 
        List of seismogram's filenames
    """

    all_seismogram = [
            filename for filename in os.listdir(folder_path)
        if not filename.startswith('G_RER') 
        if not filename.startswith('PF_TTR')  
        ]
        
    # Remove specific station on some period
    filenames = []
    all_seismogram_reduced = []
    for filename in all_seismogram:
        starttime = UTCDateTime(get_info_from_mseedname(filename)['starttime'])
        if filename.startswith('PF_NSR') and starttime < UTCDateTime("2020-02-07T10:32:10"):
            filenames.append(filename)
            
        else:
            all_seismogram_reduced.append(filename)
    return all_seismogram_reduced
