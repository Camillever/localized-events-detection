""" Module for sta-lta method on seismograms to detect start and end times of trainwaves """

import os

from typing import List, Tuple

import numpy as np

from tqdm import tqdm

from obspy.core import Stream, read
from obspy.signal.trigger import coincidence_trigger
from obspy import UTCDateTime

from locevdet.utils import get_info_from_mseedname
from locevdet.event import Event, EventList
from locevdet.stations import STATIONS_NETWORKS

def stalta_detect_events(folder_in:str, all_seismogram:List[str],
        freqmin:int, freqmax:int,
        nsta_time:int, nlta_time:int, thr_on:float, thr_off:float,
        trigger_type:str="classicstalta", thr_coincidence_sum:int=1,
        sampling_rate=100) -> EventList:
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
        trigger_type: TODO
        nsta_time: TODO
        nlta_time : TODO
        thr_on : threshold for switching trigger on
        thr_off: threshold for switching trigger off

    Returns :
        Class 'EventList' which compile all events detected by stations
    """
    # Compile the names of all seismograms

    period_filename = list(set(get_info_from_mseedname(filename, 'periodtime')
        for filename in all_seismogram))
    period_filename.sort()

    event_list = EventList()
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

            event = Event(start_global=start_global, stations=stations)
            event.add_trainwaves(stream)
            event_list.append(event)

    return event_list


def build_period_stream(period_name:str, folder_in:str, mseed_name:str, 
        stream:Stream, filter_band:Tuple[float], sampling_rate:float):
    """
    Add traces (vertical component) in the stream for a given period time

    Args:
        period_name : String of the start and end time of traces in the stream with
        folder_in : Directory path of the folder containing every Mseed seismograms
        mseed_name : Name of the Mseed file
        stream : Stream containing each trace 
        filter_band : [frequency minimum, frequency maximum], in Hertz, of the bandpass filter
        sampling_rate : Sampling rate of all trace in the stream (in Hertz)

    Returns:
        Modify the stream (by adding the traces for the given period time)
    
    """
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

def remove_too_close_trainwaves(event_list:EventList, minimum_time:float=0)-> List[str]:
    """
    Delete events too closed of each others ( with the minimum time)

    Args:
        event_list: Class which is the list of the class Event (compil each event detected)
                    See locevdet.event python file
        minimum_time: Minimum time where two trainwaves are considered close enough (in seconds)

    Returns:
        List of all event (class) which have been removed

    """
    events_removed = []

    for event in event_list:
        utc_start = UTCDateTime(event.start_global)
        for other_event in event_list:
            utc_start_other = UTCDateTime(other_event.start_global)            
            if event != other_event and np.abs(utc_start - utc_start_other) <= minimum_time:
                event_list.remove(other_event)
                events_removed.append(other_event.start_global)
    return events_removed

def remove_border_stalta_false_trainwaves(event_list:EventList, nlta_time:int, 
        uncertainty_ratio:float=0.05)-> List[str]:
    
    """ TODO
    Class 'EventList' which compile all events detected by stations, is modified : 
        without false trainwaves

    Args:
        event_list : Class which is the list of the class Event (compil each event detected)
                    See locevdet.event python file
        nlta_time : Time length of the LTA rolling window
        uncertainty_ratio : TODO

    Returns:
        List of all event (class) which have been removed
    """
    events_removed = []

    for event in event_list:
        utc_start = UTCDateTime(event.start_global)
        for _, trainwave in event.trainwaves.items():
            start_of_one_rawsignal = UTCDateTime(trainwave.trace.stats.starttime)
            break

        uncertainty_time = uncertainty_ratio * nlta_time
        limit_time = start_of_one_rawsignal + nlta_time + uncertainty_time

        if utc_start >= start_of_one_rawsignal and utc_start <= limit_time :
            event_list.remove(event)
            events_removed.append(utc_start)
    return events_removed
