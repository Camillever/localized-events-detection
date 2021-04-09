""" Module Kurtosis to detect specific and more precise start times of trainwaves, per station """

import os
import json

import matplotlib.pyplot as plt

import numpy as np

from matplotlib.widgets import Slider, Button
from obspy import UTCDateTime
from obspy.core import read
from obspy.realtime.signal import kurtosis
from obspy.signal.trigger import trigger_onset

from locevdet.utils import get_starttime_trainwave, get_period
from locevdet.visualisation.plots import update


def kurtosis_for_all_seismograms(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float,
    thr_on:float, thr_off:float,
    freqmin:float, freqmax:float,
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

    for trainwave in all_trainwaves:
        trainwave_path = os.path.join(trainwaves_path, trainwave)
        fig = plt.figure(1)
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
                        

                        trace_copy = trace.trim(
                            start_global - pre_trigger,
                            start_global + post_trigger,
                            nearest_sample=True
                        )
                        
                        trace_filtered = trace_copy.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
                        
                        # Kurtosis
                        matrix_kurtosis = kurtosis(trace_copy, win)
                        max_kurtosis = np.max(matrix_kurtosis)
                        triggertime_trainwaves_kurtosis = trigger_onset(matrix_kurtosis, thr_on*max_kurtosis, thr_off*max_kurtosis)
                        
                        all_starttimes = triggertime_trainwaves_kurtosis[:,0].copy()
                        # start_specific = lambda starttime : starttime
                        
                        #Plot
                        plt.title(f"Fenêtre de visualisation : {start_global - pre_trigger} - {start_global + post_trigger}")

                        ax1 = plt.subplot(211)

                        ax1.set_title(f'Signal en fonction du temps - {trace.stats.station} ')
                        ax1.plot(trace_copy.times(), trace_copy.data, color='grey')
                        trace_line = ax1.plot(trace_copy.times(), trace_filtered.data)
                        
                        ymin1, ymax1 = ax1.get_ylim()
                        ax1.axvline(start_global-trace.stats.starttime, ymin1, ymax1, color='brown') # Start global
                        ax1.axvline()
                        if len(triggertime_trainwaves_kurtosis)!=0 :
                            ax1.vlines(all_starttimes*trace.stats.delta, ymin1, ymax1, color='red') # Each detected start events
                        ax1.set_xlabel('temps (secondes)')

                        ax2 = plt.subplot(212, sharex=ax1)

                        ax2.set_title(f'Kurtosis - Décalage de la fenêtre : {win} secondes')
                        ax2.plot(trace_copy.times(), matrix_kurtosis)
                        xmin, xmax = ax2.get_xlim()
                        ax2.axhline(thr_on*max_kurtosis, xmin, xmax, color='red', linestyle='--')
                        ax2.axhline(thr_off*max_kurtosis, xmin, xmax, color='blue', linestyle='--')
                        ymin2, ymax2 = ax2.get_ylim()
                        ax2.axvline(start_global-trace.stats.starttime,ymin2, ymax2, color='brown')

                        # Sliders
                        ## ( https://matplotlib.org/stable/gallery/widgets/slider_demo.html)

                        axtrig_on = plt.axes([0.1, 0.25, 0.0025, 0.63], facecolor=axcolor)
                        trig_on_slider = Slider(
                            ax=axtrig_on, 
                            label='Seuil de déclenchement (%)',
                            valmin=0, 
                            valmax=1, 
                            valinit=thr_on
                        )

                        
                        # trig_on_slider.on_changed(update)

                        plt.tight_layout()
                        plt.show()

                        # Add this matrix into trainwave dictionary
                        # try:
                        #     trainwave_data["starttime_kurtosis"][station] = matrix_kurtosis
                        # except KeyError:
                        #     trainwave_data["starttime_kurtosis"] = {station : matrix_kurtosis}

                        # print(trainwave_data)

                        # if format_save == 'JSON':
                        #     with open(trainwave_path, 'w') as content:
                        #         json.dump(trainwave_data, content,  indent=1)       


        
