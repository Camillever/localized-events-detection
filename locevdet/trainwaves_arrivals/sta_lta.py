""" Module for sta-lta method on seismograms to detect start and end times of trainwaves """

import os

import numpy as np

from tqdm import tqdm
from typing import List

from obspy.core import Stream, read
from obspy.signal.trigger import coincidence_trigger
from obspy import UTCDateTime

from locevdet.utils import get_info_from_mseedname, clean_utc_str, get_starttime_trainwave
from locevdet.event import Event
from locevdet.stations import STATIONS_NETWORKS

def stalta_per_event_coincidence_trigger(folder_in:str,
        freqmin:int, freqmax:int,
        nsta_time:int, nlta_time:int, thr_on:float, thr_off:float,
        trigger_type:str="classicstalta", thr_coincidence_sum:int=1,
        save_path:str=None, format_save:str='JSON', override:bool=False, sampling_rate=100):
    """
    Save dictionaries for each trainwaves detected
    containing global start time and order of station's arrivals.

    Dectection of arrivals is calculated by the given STA-LTA coincidence_trigger method.
    Nb: Each seismograms can be filtered before STA-LTA

    Args :
        folder_in: directory path of downloaded and processed MSEED seismograms
        save_path: directory path to save the dictionnaries
        freqmin: Minimum frequency for the bandpass filter
        freqmax: Maximum frequency for the bandpass filter
        trigger_type:
        nsta_time:
        nlta_time :
        thr_on : threshold for switching trigger on
        thr_off: threshold for switching trigger off

    Returns :
        Save in the choosen format the dictionary 'trainwave'.
        'trainwave' dictionary contains:
                            - The global start time of the trainwave
                            - The order of stations which have detected first arrivals
    """
    # Compile the names of all seismograms
    all_seismogram = [
        filename for filename in os.listdir(folder_in)
    ]

    # List all perioformat = os.ds of events
    period_filename = list(set(get_info_from_mseedname(filename, 'periodtime') 
        for filename in all_seismogram))

    if save_path is None:
        save_path = os.getcwd()
    save_path_format = os.path.join(save_path, format_save)
    os.makedirs(save_path_format, exist_ok=True)

    event_list = []

    for period_name in tqdm(period_filename):
        stream = Stream()
        for mseed_name in all_seismogram:
            build_period_stream(period_name, folder_in, mseed_name,
                stream, filter_band=(freqmin, freqmax), sampling_rate=sampling_rate)

        # Event detection with STA-LTA
        nsta = int(nsta_time * sampling_rate)
        nlta = int(nlta_time * sampling_rate)
        triggered_events = coincidence_trigger(trigger_type, thr_on, thr_off, stream=stream,
            thr_coincidence_sum=thr_coincidence_sum, nsta=nsta, nlta=nlta)

        stream_stations = [trace.stats.station for trace in stream]
        # Add every single event in this period
        for triggered_event in triggered_events:
            start_global = UTCDateTime(triggered_event['time'])

            traces = [
                stream[stream_stations.index(station)]
                for station in triggered_event['stations']
            ]
            stations = [
                STATIONS_NETWORKS[trace.stats.network][trace.stats.station]
                for trace in traces
            ]
            event = Event(
                start_global=start_global,
                stations=stations
            )
            event.add_trainwaves(stream)
            trainwave_filename = '_'.join(('trainwave', clean_utc_str(start_global)))
            trainwave_filepath = os.path.join(save_path_format, trainwave_filename)
            event.save(trainwave_filepath, format_save="JSON", override=override)
            event_list.append(event)
    
    event_list.save(trainwave_filepath)
    number_trainwaves = len(os.listdir(save_path_format))
    return number_trainwaves


def build_period_stream(period_name, folder_in, mseed_name, stream, filter_band, sampling_rate):
    if period_name == get_info_from_mseedname(mseed_name, 'periodtime'):
        filepath = os.path.join(folder_in, mseed_name)
        seismogram = read(filepath)
        for _,seismo in enumerate(seismogram):
            if seismo.stats.component == 'Z':
                trace = seismo
                if trace.stats.sampling_rate != sampling_rate:
                    return
        trace.filter('bandpass', freqmin=filter_band[0], freqmax=filter_band[1])
        stream += trace

def trainwaves_too_close_remove(folder_dict:str, intervall_time:float=0):
    """
    Delete trainwaves dictionaries too close of each others in the given intervall time
    (Each supplement trainwave is considered as a duplicate)

    Args:
        folder_dict: Directory path of the trainwaves dictionnaries
        intervall_time: TODO

    Returns:
        List of dictionaries of deleted files

    """
    dict_path = os.listdir(folder_dict)
    dict_removed=[]

    for dictionary in dict_path :
        utc_starttime_dict = UTCDateTime(get_starttime_trainwave(dictionary))

        for other_dict in dict_path:

            if other_dict == dictionary:
                continue

            utc_starttime_otherdict = UTCDateTime(get_starttime_trainwave(other_dict))

            if np.abs(utc_starttime_dict - utc_starttime_otherdict) <=intervall_time:
                os.remove(os.path.join(folder_dict,dictionary))
                dict_path.remove(dictionary)

                dict_removed.append(dictionary)
    return dict_removed


def false_trainwaves(seismograms_path:str,trainwaves_path:str, 
    nlta_time:int, time_offset:int)-> List[str]:
    
    """ Delete false start time detected by STA-LTA during the LTA time-window for each seismogram 
    
    Args:
        seismograms_path: Directory path of MSEEDS seismograms
        trainwaves_path : Directory path of trainwaves dictionaries
        nlta_time : Time length of the LTA rolling window
        time_offset: Offset in seconds added to the nlta_time 
            to be certain to count false trainwaves closed to this period nlta.

    Returns:
        List of removed false trainwaves dictionaries
    """

    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
    ]

    periods = list(set(get_info_from_mseedname(filename, 'periodtime') for filename in all_seismogram))

    all_trainwaves = [
        trainwave_name for trainwave_name in os.listdir(trainwaves_path)
    ]
    dict_removed=[]
    for period in periods: 
        starttime_period = UTCDateTime(period.split('_')[0])
        prohibited_period = [starttime_period, starttime_period + nlta_time +time_offset]
        
        for trainwave_name in all_trainwaves:
            starttime_trainwave = UTCDateTime(get_starttime_trainwave(trainwave_name))
            if starttime_trainwave >= prohibited_period[0] and starttime_trainwave <= prohibited_period[1]:
                os.remove(os.path.join(trainwaves_path, trainwave_name))
                all_trainwaves.remove(trainwave_name)
                dict_removed.append(trainwave_name)
    return dict_removed
