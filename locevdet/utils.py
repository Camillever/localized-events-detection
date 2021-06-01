""" Module for utilities functions """
import os
import numpy as np
import pandas as pd
from typing import Tuple
from obspy import UTCDateTime

from sklearn.linear_model import LinearRegression
import matplotlib.dates as dates


def kurtosis_norm(trace, win_kurt_init):
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
        dictionary containing location informations of the given filename

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
