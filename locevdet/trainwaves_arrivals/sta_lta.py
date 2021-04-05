""" Module for sta-lta method on seismograms to detect start and end times of trainwaves """

import os
import json

import numpy as np

from tqdm import tqdm
from obspy.core import Stream, read
from obspy.signal.trigger import coincidence_trigger
from obspy import UTCDateTime


from locevdet.utils import get_period, clean_utc_str, get_starttime_trainwave

def stalta_per_event_coincidence_trigger(folder_in:str,
        freqmin:int, freqmax:int,
        nsta_time:int, nlta_time:int, thr_on:float, thr_off:float,
        trigger_type:str="classicstalta", thr_coincidence_sum:int=1,
        save_path:str=None, format_save:str='JSON', override:bool=False):
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

    # List all periods of events
    period_filename = get_period(all_seismogram)

    if save_path is None:
        save_path = os.getcwd()
    save_path_format = os.path.join(save_path, format_save)
    os.makedirs(save_path_format, exist_ok=True)

    # The loop
    for period_name in tqdm(period_filename):
        stream = Stream()

        for filename in all_seismogram :
            if period_name == '_'.join([filename.split('_')[2],filename.split('_')[3]]) :

                filepath = os.path.join(folder_in, filename)
                seismogram = read(filepath)

                for _,seismo in enumerate(seismogram):
                    if seismo.stats.component == 'Z':
                        trace = seismo

                #Filter Band-pass
                trace.filter('bandpass', freqmin=freqmin, freqmax=freqmax)

                stream += trace

        nsta = int(nsta_time*trace.stats.sampling_rate)
        nlta = int(nlta_time*trace.stats.sampling_rate)

        trig = coincidence_trigger(trigger_type, thr_on, thr_off, stream=stream,
            thr_coincidence_sum=thr_coincidence_sum, nsta=nsta, nlta=nlta)

        for i,_ in enumerate(trig): # Event i in the number of global detected trainwaves
            starttime_global = UTCDateTime(trig[i]['time'])
            trainwave = {
                "starttime_global" : str(starttime_global),
                "order_arrivals_detected" : trig[i]['stations']
            }

            trainwave_filename = '_'.join((
                'trainwave',
                clean_utc_str(starttime_global)
            ))

            trainwave_filepath = os.path.join(save_path_format, trainwave_filename)

            # Save in the given format
            if override or not os.path.isfile(save_path_format):
                if format_save == 'JSON':
                    with open(trainwave_filepath, 'w') as content:
                        json.dump(trainwave, content,  indent=1)

    number_trainwaves = len(os.listdir(save_path_format))
    return(print(f"{number_trainwaves} dictionnaires ont été importés"))

def trainwaves_too_close_remove(folder_dict:str, intervall_time:float=0):
    """
    Delete trainwaves dictionaries too close of each others in the given intervall time
    (Each supplement trainwave is considered as a duplicate)

    Args:
        folder_dict: Directory path of the trainwaves dictionnaries
        intervall_time: TODO

    Returns:
        Delete duplicate trainwaves
        and returns the number and names of deleted files
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
    return(print(f"{len(dict_removed)} dictionnaires ont été supprimés, dont : {dict_removed}"))

def false_trainwaves(seismograms_path:str, trainwaves_path:str, nlta_time:int, time_offset:int):
    """ TODO  """

    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
    ]

    periods = get_period(all_seismogram)

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
    return(print(f"{len(dict_removed)} dictionnaires ont été supprimés, dont : {dict_removed}"))
