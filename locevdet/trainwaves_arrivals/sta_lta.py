""" Module for sta-lta method on seismograms to detect start and end times of trainwaves """

import os

from typing import List, Tuple

import numpy as np

from tqdm import tqdm

from obspy.core import Stream, read
from obspy.signal.trigger import coincidence_trigger
from obspy import UTCDateTime

from locevdet.utils import get_info_from_mseedname
from locevdet.event import Event
from locevdet.eventlist import EventList
from locevdet.stations import STATIONS_NETWORKS
from locevdet.examples.casse_riv_est.download import apodisation
from locevdet.descriptors import envelope_fct, snr_calculation_fct

def stalta_detect_events(folder_in:str, all_seismogram:List[str],
        freqmin:int, freqmax:int,
        nsta_time:int, nlta_time:int, thr_on:float, thr_off:float,
        minimum_time:float, uncertainty_ratio:float=5/100,
        trigger_type:str="classicstalta", thr_coincidence_sum:int=1,
        sampling_rate=100, low_snr_remove:bool=False, **kwargs) -> EventList:
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
        low_snr_remove : if True, calculate the signal to noise ratio (snr) of each trace
                and remove events with low snr
                !Need: rolling_max_window, general_thr, pre_trigger and post_trigger to trim traces!

        kwargs: other arguments as:
            - rolling_max_window : length of the moving window (in seconds) to calculate envelope
            - pre_trigger : duration before the start_global to cut the trace, in seconds
                    (see trim_trace function from waveform_processing.py )
            - post_trigger : duration after the start_global to cut the trace, in seconds
                    (see trim_trace function from waveform_processing.py )

    Returns :
        Class 'EventList' which compile all events detected by stations
    """
    # Compile the names of all seismograms
    lowsnr_removed = []
    period_filename = list(set(get_info_from_mseedname(filename)['periodtime']
        for filename in all_seismogram))
    period_filename.sort()

    event_list = EventList()
    for period_name in tqdm(period_filename):
        stream_traces = Stream()
        stream_traces_filt = Stream()
        for mseed_name in all_seismogram:
            build_period_stream(period_name, folder_in, mseed_name,
                stream_traces, stream_traces_filt, 
                filter_band=(freqmin, freqmax), sampling_rate=sampling_rate)
        
        # Event detection with STA-LTA
        nsta = int(nsta_time * sampling_rate)
        nlta = int(nlta_time * sampling_rate)
        triggered_events = coincidence_trigger(trigger_type, thr_on, thr_off, stream=stream_traces_filt,
            thr_coincidence_sum=thr_coincidence_sum, nsta=nsta, nlta=nlta)

        stream_stations = [trace.stats.station for trace in stream_traces]
        # Add every single event in this period
        

        for triggered_event in triggered_events:
            start_global = UTCDateTime(triggered_event['time'])
            traces = [
                stream_traces[stream_stations.index(station)]
                for station in triggered_event['stations']
            ]

            stations = [
                STATIONS_NETWORKS[trace.stats.network][trace.stats.station]
                for trace in traces
            ]

            event_accepted = True
            # Calculate snr before creation of events
            if low_snr_remove is True:
                general_thr = kwargs.get('general_thr', 3)
                all_snr = []
                stations_order = []
                for num, trace in enumerate(traces):
                    stations_order.append(stations[num].name)
                    rolling_max_window = kwargs.get('rolling_max_window', 0)

                    envelope_trim_filt = envelope_fct(
                        trace_type='trimmed_filtered', 
                        rolling_max_window=rolling_max_window,
                        trace=trace, 
                        freqmin = freqmin,
                        freqmax = freqmax,
                        pre_trigger = kwargs.get('pre_trigger', 10),
                        post_trigger = kwargs.get('post_trigger', 50),
                        start_global = start_global)

                    envelope_filt = envelope_fct(
                        trace_type='trace_filtered', 
                        rolling_max_window=rolling_max_window, 
                        trace=trace, freqmin = freqmin, freqmax = freqmax)

                    _, snr = snr_calculation_fct(
                        envelope_trim_filt=envelope_trim_filt, envelope_filt=envelope_filt)
                    all_snr.append(snr)

                if all(i <= general_thr for i in all_snr) :
                    event_accepted = False
                    lowsnr_removed.append(start_global)

                elif 'PER' in stations_order:
                    index_per = stations_order.index('PER')
                    snr_per = all_snr[index_per]
                    all_snr_copy = all_snr.copy()
                    all_snr_copy.pop(index_per)

                    if snr_per > general_thr and all(i <= general_thr for i in all_snr_copy):
                        event_accepted = False
                        lowsnr_removed.append(start_global)

            starttime = stream_traces[0].stats.starttime
            endtime = stream_traces[0].stats.endtime

            # To not save event which are too close of borders (affected by apodization)
            apodization_time_restricted = apodisation*(endtime-starttime)

            if start_global <= (endtime - apodization_time_restricted) and \
                start_global >= (starttime + apodization_time_restricted) and \
                    event_accepted is True:
                event = Event(start_global=start_global, stations=stations)
                event.add_trainwaves_from_stalta(stream_traces)
                event_list.append(event)

                for _, trainwave in event.trainwaves.items():
                    station_name = trainwave.station.name
                    trainwave.freqmin_interest = freqmin
                    trainwave.freqmax_interest = freqmax
                    if low_snr_remove is True:
                        index_station = stations_order.index(station_name)
                        snr_station = all_snr[index_station]
                        trainwave.snr = snr_station


    # Remove duplicate and false events from event_list
    duplicate_ev_removed = remove_too_close_trainwaves(event_list, minimum_time)
    false_ev_removed = remove_border_stalta_false_trainwaves(event_list, nlta_time, uncertainty_ratio)

    return event_list, duplicate_ev_removed, false_ev_removed, lowsnr_removed


def build_period_stream(period_name:str, folder_in:str, mseed_name:str, 
        stream_traces:Stream, stream_traces_filt:Stream, 
        filter_band:Tuple[float], sampling_rate:float):
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
    if period_name == get_info_from_mseedname(mseed_name)['periodtime']:
        filepath = os.path.join(folder_in, mseed_name)
        seismogram = read(filepath)
        

        for _,seismo in enumerate(seismogram):
            if seismo.stats.component == 'Z':
                trace = seismo
                stream_traces += trace
                if trace.stats.sampling_rate != sampling_rate:
                    return
        trace_copy = trace.copy()
        trace_filt = trace_copy.filter('bandpass', freqmin=filter_band[0], freqmax=filter_band[1])
        stream_traces_filt += trace_filt

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
