""" Module Kurtosis to detect specific and more precise start times of trainwaves, per station """

import os
import json

import matplotlib.pyplot as plt
import numpy as np

from obspy import UTCDateTime
from obspy.core import read
from obspy.realtime.signal import kurtosis
from obspy.signal.trigger import trigger_onset

from locevdet.utils import get_starttime_trainwave, get_period
from locevdet.waveform_processing import trim_trace

# Utils 

def starttimes_trigger_by_kurtosis(matrix_kurtosis, thr_on, thr_off):
    """ TODO
    
    """
    max_kurtosis = np.max(matrix_kurtosis)
    triggertime_trainwaves_kurtosis = trigger_onset(matrix_kurtosis, thr_on*max_kurtosis, thr_off*max_kurtosis) 
    all_starttimes = triggertime_trainwaves_kurtosis[:,0].copy()
    return all_starttimes

def kurtosis_norm(trace, win_kurt_init):
    """ TODO
    
    """
    matrix_kurtosis = kurtosis(trace, win_kurt_init)
    norm = np.linalg.norm(matrix_kurtosis)
    kurt_norm = matrix_kurtosis/norm
    return kurt_norm


# View all seismograms 
from locevdet.visualisation.plots import kurt_param_sliders_per_trace

def kurtosis_for_all_seismograms(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float,
    win:float=3, format_save:str="JSON"):
    """ TODO
    
    """
    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
        if not filename.startswith('G_RER') 
        if not filename.startswith('PF_TTR')
    ]
    
    all_trainwaves = [
        trainwave_filename for trainwave_filename in os.listdir(trainwaves_path)
    ]
    # Initial value of thresholds (25th and 75th purcentiles defined by RangeSlider)
    thr_on_init = 0.75
    thr_off_init = 0.25

    for trainwave in all_trainwaves:
        trainwave_path = os.path.join(trainwaves_path, trainwave)
        
        with open(trainwave_path) as json_file:
            trainwave_data = json.load(json_file)
            start_global = UTCDateTime(trainwave_data['starttime_global'])

            for filename in all_seismogram:
                period = get_period([filename])
                seismo_start = UTCDateTime(period[0].split('_')[0])
                seismo_end = UTCDateTime(period[0].split('_')[1])

                for station in trainwave_data['order_arrivals_detected']:
                    if start_global > seismo_start and start_global < seismo_end and filename.split('_')[1] == station:
                        filepath = os.path.join(seismograms_path, filename)
                        seismogram = read(filepath)

                        for _,seismo in enumerate(seismogram):
                            if seismo.stats.component == 'Z':
                                trace = seismo

                        trace_trim = trim_trace(trace, start_global, pre_trigger, post_trigger)

                        # Kurtosis
                        kurt_norm = kurtosis_norm(trace_trim, win)                    
                        all_starttimes_init = starttimes_trigger_by_kurtosis(kurt_norm, thr_on_init, thr_off_init)

                        kurt_param_sliders_per_trace(
                            trace_trim, 
                            start_global, 
                            kurt_norm, 
                            all_starttimes_init, 
                            thr_off_init, 
                            thr_on_init, 
                            win
                        )

                        
                        # # Add this matrix into trainwave dictionary
                        # # try:
                        # #     trainwave_data["starttime_kurtosis"][station] = matrix_kurtosis
                        # # except KeyError:
                        # #     trainwave_data["starttime_kurtosis"] = {station : matrix_kurtosis}

                        # # print(trainwave_data)

                        # # if format_save == 'JSON':
                        # #     with open(trainwave_path, 'w') as content:
                        # #         json.dump(trainwave_data, content,  indent=1)       

