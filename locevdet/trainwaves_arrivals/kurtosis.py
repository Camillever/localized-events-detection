""" Module Kurtosis to detect specific and more precise start times of trainwaves, per station """

import os
import json

import matplotlib.pyplot as plt

import numpy as np

from obspy import UTCDateTime
from obspy.core import read
from obspy.realtime.signal import kurtosis

from locevdet.utils import get_starttime_trainwave, get_period


def kurtosis_for_all_seismograms(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float,
    freqmin:float, freqmax:float,
    win:float=3, format_save:str="JSON"):
    """ TODO
    
    """

    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
    ]

    all_trainwaves = [
        trainwave_filename for trainwave_filename in os.listdir(trainwaves_path)
    ]

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

                        # Filter
                        trace.filter('bandpass', freqmin=freqmin, freqmax=freqmax)

                        trace_copy = trace.trim(
                            start_global - pre_trigger,
                            start_global + post_trigger,
                            nearest_sample=True
                        )

                        # Kurtosis
                        matrix_kurtosis = kurtosis(trace_copy, win)
                        
                        # print(np.max(matrix_kurtosis))
                        plt.title(f"Fenêtre de visualisation : {start_global - pre_trigger} - {start_global + post_trigger}")
                        ax1 = plt.subplot(211)
                        ax1.set_title(f'Signal en fonction du temps - {trace.stats.station} ')
                        ax1.plot(trace_copy.times(), trace_copy.data)
                        ax1.set_xlabel('temps (secondes)')
                        ax2 = plt.subplot(212)
                        ax2.set_title(f'Kurtosis - Fenêtre glissante : {win} secondes')
                        ax2.plot(matrix_kurtosis)
                        plt.tight_layout()
                        plt.show()
                        matrix_tolist_kurtosis = matrix_kurtosis.tolist()

                        # Add this matrix into trainwave dictionary
                        try:
                            trainwave_data["matrix_kurtosis"][station] = matrix_kurtosis
                        except KeyError:
                            trainwave_data["matrix_kurtosis"] = {station : matrix_kurtosis}

                        # print(trainwave_data)

                        # if format_save == 'JSON':
                        #     with open(trainwave_path, 'w') as content:
                        #         json.dump(trainwave_data, content,  indent=1)       


        
