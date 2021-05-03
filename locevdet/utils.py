""" Module for utilities functions """
import os
import numpy as np
import pandas as pd
from typing import Tuple, List
from obspy import UTCDateTime



def kurtosis_panda(trace, win_kurt_init):
    """ Apply recursive kurtosis calculation on the given trace and returns a normalized kurtosis matrix.
    (See : https://pandas.pydata.org/pandas-docs/version/0.25.3/reference/api/pandas.core.window.Rolling.kurt.html )

    Args:
        trace : Trace object to append to this RtTrace
        win_kurt_init : window length in seconds for the kurtosis (shift length ?)

    Returns:
        Normalized npdarray of the kurtosis matrix (trace divised by the max of the trace).    
    """

    data_pd = pd.Series(trace.data)
    kurt = data_pd.rolling(win_kurt_init).kurt()

    # Conversion to numpy array and remove 'Nan'
    kurt_np = kurt.to_numpy()
    kurt_np = np.nan_to_num(kurt_np)

    # Normalization
    kurt_max = np.max(kurt_np)
    kurt_norm = []
    if kurt_max != 0:
        kurt_norm = kurt_np/ kurt_max   
    return kurt_norm

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


def ref_duration_to_time_window(time_reference:UTCDateTime, duration:float=0, 
    endtime:UTCDateTime=None,
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

    if endtime is None:
        end_time = time_reference + duration + time_offset[1]
    else:
        end_time = endtime + time_offset[1]
    return start_time, end_time

def get_info_from_mseedname(filename:str, info_type:str)-> str:
    """
    Args:
        filename : Name of the mseed file
            with the nomenclature : '{network}_{station}_{starttime}_{endtime}'
        info_type : String of type of information 
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

def rooling_max(signal, win_lenght=None):
    if win_lenght is None:
        win_lenght = len(signal) // 20

    max_signal = -np.inf * np.ones_like(signal)
    for i in range(len(signal)):
        if len(signal) - i < win_lenght:
            max_signal[i] = np.max(signal[i:])
        else:
            max_signal[i] = np.max(signal[i:i+win_lenght])

    return max_signal




from locevdet.event import *
from locevdet.visualisation.plots import hist_band_freq

def freq_band_interest(eventlist:EventList, save_fig_path:str=None, show:bool=True):
    """ TODO


    """

    all_central_frq = []
    all_fmin = []

    for event in eventlist:
        for _, trainwave in event.trainwaves.items():
            if trainwave.matlab_data is not None : 
                all_fmin.append(trainwave.matlab_data['trainwave']['fmin'])
                all_central_frq.append(trainwave.matlab_data['trainwave']['centralfrequency'])

    mean_fmin = np.mean(all_fmin)
    fc = np.max(all_central_frq)

    # Histograms
    if save_fig_path is not None or show is True:
        hist_band_freq(
            all_fmin, 
            all_central_frq, 
            fmin=mean_fmin, 
            fc=fc, 
            save_fig_path=save_fig_path, 
            show=show)
        
    return float(format(mean_fmin, '.2f')), float(format(fc, '.2f'))

# def ti_to_utcdatetime(filename_matlab:str, ti:float):
#     """
#     Conversion of trainwave start in seconds to UTCDateTime

#     Args:
#         filename_matlab : matlab file name 
#             with nomenclature 'Ret_{network}_{station}_{starttime}_{endtime}'
#         ti : Start of the wave train in seconds refered from the start of the seismogram (filename)

#     Returns :
#         ti_utc_datetime : Start of the wave train in UTCDateTime
#     """
#     starttime = UTCDateTime(get_info_from_matname(filename_matlab)['starttime'])
#     ti_utc_datetime = starttime + ti

#     return ti_utc_datetime

# def similarity_ti(matlab_folder_path:str, mseeds_path:str, stations:list):
#     """TODO
#     """
#     all_data_mat = [ 
#         mat for mat in os.listdir(matlab_folder_path)
#         if mat.endswith('.mat')
#         ]
 
#     all_seismogram = os.listdir(mseeds_path)
#     all_ti_from_mat = []

#     for matname in all_data_mat:
#         print(matname)
#         ti = ti_matlab_results(matlab_folder_path, matname)
#         print("ti :", ti)
#         ti_utc_datetime = ti_to_utcdatetime(matname, ti)
#         all_ti_from_mat.append(ti_utc_datetime)
#     return all_ti_from_mat

