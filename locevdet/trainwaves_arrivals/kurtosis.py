""" Module Kurtosis to detect specific and more precise start times of trainwaves, per station """

import os
import json

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import Slider, RangeSlider, Button
import mpl_interactions.ipyplot as iplt
from mpl_interactions import interactive_axhline

from obspy import UTCDateTime
from obspy.core import read
from obspy.realtime.signal import kurtosis
from obspy.signal.trigger import trigger_onset

from locevdet.utils import get_starttime_trainwave, get_period
# from locevdet.visualisation.plots import update_trig

def starttimes_trigger_by_kurtosis(matrix_kurtosis, max_kurtosis, thr_on, thr_off):
    triggertime_trainwaves_kurtosis = trigger_onset(matrix_kurtosis, thr_on*max_kurtosis, thr_off*max_kurtosis) 
    all_starttimes = triggertime_trainwaves_kurtosis[:,0].copy()
    return all_starttimes




def kurtosis_for_all_seismograms(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float,
    freqmin:float, freqmax:float,
    thr_on:float=0.15, 
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

    thr_off = thr_on-0.05

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
                        
                        
                        
                        # Kurtosis
                        matrix_kurtosis = kurtosis(trace_copy, win)
                        print(matrix_kurtosis)
                        max_kurtosis = np.max(matrix_kurtosis)
                        print("max :", max_kurtosis)
                    
                        all_starttimes = starttimes_trigger_by_kurtosis(matrix_kurtosis, max_kurtosis, thr_on, thr_off)
                        
                        
                        #Plot
                        plt.title(f"Fenêtre de visualisation : {start_global - pre_trigger} - {start_global + post_trigger}")

                        ax1 = plt.subplot(411)
                        ax1.set_title(f'Signal en fonction du temps - {trace.stats.station} ')
                        trace_line = ax1.plot(trace_copy.times(), trace_copy.data, color='grey')
                        trace_line = ax1.plot(trace_copy.times(), trace_copy.data)

                        
                        
                        ymin1, ymax1 = ax1.get_ylim()
                        ax1.set_xlabel('temps (secondes)')

                        ax2 = plt.subplot(412, sharex=ax1)
                        ax2.set_title(f'Kurtosis - Décalage de la fenêtre : {win} secondes')
                        ax2.plot(trace_copy.times(), matrix_kurtosis)
                        xmin, xmax = ax2.get_xlim()
                        ymin2, ymax2 = ax2.get_ylim()

                        

                        # Start global of trainwave
                        ax1.axvline(start_global-trace.stats.starttime,ymin1, ymax1, color='brown')
                        ax2.axvline(start_global-trace.stats.starttime,ymin2, ymax2, color='brown')
                            
                        
                        if len(all_starttimes)!=0 :
                            for _, starttime in enumerate(all_starttimes):
                                index_times = trace.times()[int(starttime)]
                                start_vlineax1 = ax1.axvline(index_times, ymin1, ymax1, color='red') # Each detected start events
                                start_vlineax2 = ax2.axvline(index_times, ymin2, ymax2, color='red')

                        kurt_on2 = ax2.axhline(thr_on*max_kurtosis, 
                            xmin, 
                            xmax, 
                            color='red', 
                            linestyle='--')

                        # SLIDERS
                        ## Filter 
                        
                        axfreq = plt.axes([0.25, 0.2, 0.65, 0.03])
                        freq_slider = RangeSlider(
                            ax=axfreq, 
                            label='Frequency (Hz)',
                            valmin=0.05, 
                            valmax=50,
                            valinit=[freqmin, freqmax]
                        )

                        def update_filter(val):
                            trace_filtered = trace_copy.filter('bandpass', freqmin=val[0], freqmax=val[1])
                            trace_line.set_ydata(trace_filtered.data)                                                        
                            fig.canvas.draw_idle()

                        # Window for kurtosis
                        axwin = plt.axes([0.2, 0.1, 0.65, 0.03])
                        win_slider = Slider(
                            ax=axwin, 
                            label='Kurtosis : Window shift (s)',
                            valmin=0.02, 
                            valmax=0.05,
                            valinit=win
                        )

                        def update_win(val):
                            win = val
                            matrix_kurtosis = kurtosis(trace_filtered, win)
                            

                        ## Trigger on
                        
                        axtrig_on = plt.axes([0.15, 0.1, 0.65, 0.03])
                        trig_slider = RangeSlider(
                            ax=axtrig_on, 
                            label='Threshold (%)',
                            valmin=0, 
                            valmax=1,
                            valinit=[thr_off, thr_on] 
                        )

                        
                        def update_trig(val):
                            thr_off_val = val[0]*max_kurtosis
                            thr_on_val = val[1]*max_kurtosis
                            kurt_on2.set_ydata(thr_on_val)

                            # Triggered start times
                            matrix_kurtosis = kurtosis(trace_filtered, win)
                            max_kurtosis = np.max(matrix_kurtosis)
                            all_starttimes = starttimes_trigger_by_kurtosis(matrix_kurtosis, max_kurtosis, val[0], val[1])
                            for _, starttime in enumerate(all_starttimes):
                                index_times = trace.times()[int(starttime)]
                                start_vlineax1.set_xdata(index_times)
                                start_vlineax2.set_xdata(index_times)
                                                        
                            fig.canvas.draw_idle()

                        win_slider.on_changed(update_win)
                        trig_slider.on_changed(update_trig)
                        freq_slider.on_changed(update_filter)

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


        
