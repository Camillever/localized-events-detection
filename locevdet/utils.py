""" Module for utilities functions """
import os
from typing import Tuple
import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset

from sklearn.linear_model import LinearRegression
import matplotlib.dates as dates

def kurtosis_norm(trace, win_kurt_init):
    """ Apply recursive kurtosis calculation on the given trace
    and returns a normalized kurtosis matrix.
    (See : https://pandas.pydata.org/pandas-docs/version/0.25.3/reference/api/pandas.core.window.Rolling.kurt.html )

    Args:
        trace : Trace object to append to this RtTrace
        win_kurt_init : slidding time window in seconds for the kurtosis

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

def starttimes_trigger_by_kurtosis(matrix_kurtosis, thr_on, thr_off):
    """
    Calculate all possible starts for a trainwave from the given kurtosis matrix of this trainwave
    and threshold to trigger the event

    Args :
        matrix_kurtosis : Array of the kurtosis of the trace
        thr_on : threshold on to detect the start of the event
        thr_off : threshold off to detect the end of the event

    Returns :
        The list of all starts of the given trace
    """
    max_kurtosis = np.max(matrix_kurtosis)
    triggertime_trainwaves_kurtosis = trigger_onset(matrix_kurtosis, thr_on*max_kurtosis, thr_off*max_kurtosis)
    all_starttimes = triggertime_trainwaves_kurtosis[:,0].copy()
    return all_starttimes


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
    Calculate the start and the end of a seismogram (permit to cut before downloading)

    Args:
        time_reference: Start global of the event in UTCDateTime
        duration: duration of the event in seconds
        time_offset: Tuple of time (seconds) to add before and after the time_reference

    Returns:
        The Tuple : start and end in UTCDateTime
    """
    start_time = time_reference - time_offset[0]

    if endtime is None:
        end_time = time_reference + duration + time_offset[1]
    else:
        end_time = endtime + time_offset[1]
    return start_time, end_time

def get_info_from_mseedname(filename:str)-> dict:
    """
    Args:
        filename : Name of the mseed file
            with the nomenclature : '{network}_{station}_{starttime}_{endtime}'
    Returns:
        The dictionary of all information given by the seismogram's filename

        NB : periodtime is a string as "{str(starttime)}_{str(endtime)}"
    """
    title_mseed = filename.split('_')
  
    seismogram_info = {
        'network': title_mseed[0],
        'station': title_mseed[1],
        'starttime': UTCDateTime(title_mseed[2]),
        'endtime': UTCDateTime(title_mseed[3]),
        'periodtime': '_'.join((str(title_mseed[2]),str(title_mseed[3])))
    }
    return seismogram_info

def localisation(filename:str):
    """
    Give WGS94 (longitude, latitude, z) and RGR92 (x, y, z) coordinates
    of the station from a given MSEED filename

    Args:
        filename : name of the MSEED file with the nomenclature :
            {network}_{station}_{starttime}_{endtime}

    Returns:
        The dictionary containing location informations of the given filename

    """
    positions_station = pd.read_csv(
        os.path.join('Caracteristics_stations', 'positions_stations.csv'),
        sep=';',
        index_col='Station')
    station = get_info_from_mseedname(filename)['station']
    location = {
        "latitude" : positions_station.loc[str(station), "Latitude"],
        "longitude" : positions_station.loc[str(station), "Longitude"],
        "x" : positions_station.loc[str(station), "x"],
        "y" : positions_station.loc[str(station), "y"],
        "z": positions_station.loc[str(station), "z"]
    }
    return location

def get_starttime_trainwave(filename_trainwave:str):
    """
    Give the start time string of the trainwave dictionary name

    Args:
        filename_trainwave: file name of the trainwave dictionary

    Returns:
        String of the start time of this trainwave

    """
    return filename_trainwave.split('_')[-1]

def rolling_max(signal, win_lenght=None):
    """ TODO """
    if win_lenght is None:
        win_lenght = len(signal) // 20

    max_signal = -np.inf * np.ones_like(signal)
    for i in range(len(signal)):
        if len(signal) - i < win_lenght:
            max_signal[i] = np.max(signal[i:])
        else:
            max_signal[i] = np.max(signal[i:i+win_lenght])

    return max_signal

def linear_regression_on_dates(X, Y):
    """ TODO """
    X_train = np.array(dates.date2num(X)).reshape(-1, 1)
    Y_train = np.array(dates.date2num(Y)).reshape(-1, 1)
    modeleReg=LinearRegression()
    modeleReg.fit(X_train, Y_train)
    Y_pred = modeleReg.predict(X_train)
    return X_train, Y_train, Y_pred, modeleReg

def skipper(fname:str, header:bool=False):
    """ Permit to skip header on a txt file

    Args:
        fname : directory path of the txt file
        header : inform if the txt file has header or not
    Returns:
        read the txt file from row without comment or description (header)
    """
    with open(fname) as fin:
        no_comments = (line for line in fin if not line.lstrip().startswith('#'))
        if header:
            next(no_comments, None) # skip header
        for row in no_comments:
            yield row
