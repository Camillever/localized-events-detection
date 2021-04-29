""" Module for utilities functions """

from typing import Tuple, List
from obspy import UTCDateTime

def clean_utc_str(utc_datetime:UTCDateTime) -> str:
    """
    To reduce the length of date (in UTCDateTime) and return a string of the date
    (Useful for nomenclature/name of file)

    Args:
        utc_datetime: date in UTCDateTime (YYYY-MM-DDThh:mm:ss.xxxx)

    Returns:
        A string of the date as "YYYY-MM-DDThh-mm-ss"

    """
    return str(utc_datetime).split('.')[0].replace(':','-')


def ref_duration_to_time_window(time_reference:UTCDateTime, duration:float,
        time_offset:Tuple[int]=(0, 0)) -> Tuple[float]:
    """
    TODO

    Args:
        time_reference: TODO
        duration: TODO
        time_offset: TODO

    Returns:
        TODO
    """
    start_time = time_reference - time_offset[0]
    end_time = time_reference + duration + time_offset[1]
    return start_time, end_time

def get_info_from_mseedname(filename:str, type_info:str)-> str:
    """
    Args:
        filename : Name of the mseed file
            with the nomenclature : '{network}_{station}_{starttime}_{endtime}'
        info_type : String of type of info 
            Possibility of values :
                - 'network'
                - 'station'
                - 'starttime' (of the mseed signal)
                - 'endtime' (of the mseed signal)
                - 'periodtime' ('starttime_endtime' of the mseed signal)
    Returns:
        The given type info (string)
    """
    if info_type == 'network':
        return str(filename.split('_')[0])
    elif info_type == 'station':
        return str(filename.split('_')[1])
    elif info_type == 'starttime':
        return str(filename.split('_')[2])
    elif info_type == 'endtime':
        return str(filename.split('_')[3])  
    elif info_type == 'periodtime':
        starttime = get_info_from_mseedname(filename, 'starttime')
        endtime = get_info_from_mseedname(filename, 'endtime')
        return "_".join((starttime, endtime))

def get_starttime_trainwave(filename_trainwave:str):
    """
    Give the start time string of the trainwave dictionary name

    Args:
        filename_trainwave: file name of the trainwave dictionary

    Returns:
        String of the start time of this trainwave

    """
    return filename_trainwave.split('_')[-1]


def freq_band_interest(matlab_folder_path:str, save_fig_path:str=None):
    """ TODO


    """
    import os
    import numpy as np
    from mat4py import loadmat

    files = [
        f for f in os.listdir(matlab_folder_path)
        if f.endswith('.mat')
    ]
    all_central_frq = []
    all_fmin = []
    for file in files:
        file_path = os.path.join(matlab_folder_path, file)
        mat = loadmat(file_path)

        # Compile all values of fmin and fc
        central_frq = mat['wavetrains']['frwidth']
        fmin = mat['fmin']
        
        all_fmin.append(fmin)
        for fq in central_frq:
            all_central_frq.append(fq)

    # Calculation
    mean_fmin = np.mean(all_fmin)
    fc = np.max(all_central_frq)

    # Histograms
    hist_band_freq(fmin=mean_fmin, fc=fc, save_fig_path=save_fig_path, show=False)
    
    return float(format(mean_fmin, '.2f')), float(format(fc, '.2f'))

def ti_matlab_results(matlab_folder_path, filename_matlab:str):
    """
    TODO
    """
    import os
    from mat4py import loadmat

    file_path = os.path.join(matlab_folder_path, filename_matlab)
    mat = loadmat(file_path)

    ti = mat['wavetrains']['ti'][0]
    return ti

def ti_to_utcdatetime(filename_matlab:str, ti:float):
    """
    Conversion of trainwave start in seconds to UTCDateTime

    Args:
        filename_matlab : matlab file name 
            with nomenclature 'Ret_{network}_{station}_{starttime}_{endtime}'
        ti : Start of the wave train in seconds refered from the start of the seismogram (filename)

    Returns :
        ti_utc_datetime : Start of the wave train in UTCDateTime
    """
    starttime = UTCDateTime(get_info_from_matname(filename_matlab, 'starttime'))
    ti_utc_datetime = starttime + ti

    return ti_utc_datetime

def get_info_from_matname(filename_mat:str, info_type:str='starttime')-> str:
    """
    Args:
        filename_mat : Name of the matlab file 
            with the nomenclature : 'Ret_{network}_{station}_{starttime}_{endtime}.mat'
        info_type : String of type of info 
            Possibility of values :
                - 'network'
                - 'station'
                - 'starttime' (of the mseed signal)
                - 'endtime' (of the mseed signal)
                - 'periodtime' ('starttime_endtime' of the mseed signal)
    Returns:
        The given type info (string)
    """
    if info_type == 'network':
        return str(filename_mat.split('_')[1])
    elif info_type == 'station':
        return str(filename_mat.split('_')[2])
    elif info_type == 'starttime':
        return str(filename_mat.split('_')[3])
    elif info_type == 'endtime':
        endtime = filename_mat.split('_')[4]
        return str(endtime.split('.')[0])  # To remove '.mat'
    elif info_type == 'periodtime':
        starttime = get_info_from_matname(filename_mat, 'starttime')
        endtime = get_info_from_matname(filename_mat, 'endtime')
        return "_".join((starttime, endtime))