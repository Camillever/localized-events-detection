""" Module Kurtosis_ test """

import os
import json

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import RangeSlider, Slider

class RangeSlider(RangeSlider):

    def __init__(self, ax, label, valmin, valmax, **kwargs):
        self.val = (valmin, valmax)
        super().__init__(ax, label, valmin, valmax, **kwargs)

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

def kurtosis_norm(trace, win_kurt_init:float):
    """ Apply recursive kurtosis calculation on the given trace and returns a normalized kurtosis matrix.
    (See : https://docs.obspy.org/packages/autogen/obspy.realtime.signal.kurtosis.html )

    Args:
        trace : Trace object to append to this RtTrace
        win_kurt_init : window length in seconds for the kurtosis (shift length ?)
    
    Returns:
        Normalized (by numpy.linalg.norm) npdarray of the kurtosis matrix.    
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
        
        numero = 0

        with open(trainwave_path) as json_file:
            trainwave_data = json.load(json_file)
            start_global = UTCDateTime(trainwave_data['starttime_global'])
            lists_stations = trainwave_data['order_arrivals_detected'].copy()
            if 'TTR' in lists_stations :
                lists_stations.remove('TTR')
            elif 'RER' in lists_stations :
                lists_stations.remove('RER')
            print(lists_stations)
            fig, ax = plt.subplots(nrows=2, ncols=len(lists_stations))
            numero = -1

            for filename in all_seismogram:
                period = get_period([filename])
                seismo_start = UTCDateTime(period[0].split('_')[0])
                seismo_end = UTCDateTime(period[0].split('_')[1])
                
                for station in lists_stations:
                    
                    
                    if start_global > seismo_start and start_global < seismo_end and filename.split('_')[1] == station:
                        numero += 1
                        filepath = os.path.join(seismograms_path, filename)
                        seismogram = read(filepath)

                        for _,seismo in enumerate(seismogram):
                            if seismo.stats.component == 'Z':
                                trace = seismo

                        trace_trim = trim_trace(trace, start_global, pre_trigger, post_trigger)

                        # Kurtosis
                        kurt_norm = kurtosis_norm(trace_trim, win)                    
                        all_starttimes_init = starttimes_trigger_by_kurtosis(kurt_norm, thr_on_init, thr_off_init)
                        
                        
                        # print(numero)
                        # kurt_param_sliders_per_trace(
                        #     trace_trim, 
                        #     start_global, 
                        #     kurt_norm, 
                        #     all_starttimes_init, 
                        #     thr_off_init, 
                        #     thr_on_init, 
                        #     win,
                        #     fig,
                        #     numero, 
                        #     len(lists_stations)
                        # )
                        # start_name_trace = clean_utc_str(trace.stats.starttime)
                        # end_name_trace = clean_utc_str(trace.stats.endtime)

                        max_kurtosis = np.max(kurt_norm)

                        # fig = plt.figure(f"{trace.stats.station} : {start_name_trace} - {end_name_trace}")
                        # fig.suptitle(f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}")

                        ax[numero] = fig.add_subplot(2,nb_stations, numero)
                        # ax1 = plt.subplot(211)
                        ax[numero].set_title(f'Signal en fonction du temps - {trace_trim.stats.station} ')
                        ax[numero].plot(trace_trim.times(), trace_trim.data, color='grey')
                        trace_line = ax[numero].plot(trace_trim.times(), trace_trim.data)
                        ax[numero].set_xlabel('temps (secondes)')
                        
                        # Start global of trainwave
                        ymin1, ymax1 = ax[numero].get_ylim()
                        offset = start_global-trace_trim.stats.starttime
                        ax[numero].axvline(offset, ymin1, ymax1, color='brown')  # Is not working ??
                        
                        # ax2 = fig.add_subplot(2,nb_stations, nb_stations + numero_col_subplot)
                        ax[nb_stations + numero] = plt.subplot(212, sharex=ax[numero])
                        ax[nb_stations + numero] .set_title(f'Kurtosis - Décalage de la fenêtre : {win_init} secondes')
                        kurto_line = ax[nb_stations + numero] .plot(trace_trim.times(), kurt_norm)
                        ax[nb_stations + numero] .set_ylim(top=1, bottom=0)
                        xmin, xmax = ax[nb_stations + numero] .get_xlim()

                        # Start global of trainwave
                        ymin2, ymax2 = ax[nb_stations + numero] .get_ylim()
                        ax[nb_stations + numero] .axvline(offset, ymin2, ymax2, color='brown')

                        # Kurtosis threshold and start times triggered
                        vline = ax[nb_stations + numero] .vlines(all_starttimes_init*trace.stats.delta, ymin2, ymax2, color='red')
                        kurt_on2 = ax[nb_stations + numero] .axhline(thr_on_init*max_kurtosis, xmin, xmax, color='red',linestyle='--')

                        # SLIDERS

                        # ## Filter 
                        # axfreq = plt.axes([0.25, 0.2, 0.65, 0.03])
                        # freq_slider = RangeSlider(
                        #     ax=axfreq, 
                        #     label='Frequency (Hz)',
                        #     valmin=0.05, 
                        #     valmax=49.5,
                        #     valinit=(0.05, 49.5)
                        # )
                        
                        # ## Window for kurtosis
                        # axwin = plt.axes([0.25, 0.15, 0.65, 0.03])
                        # win_slider = Slider(
                        #     ax=axwin, 
                        #     label='Kurtosis : Window shift (s)',
                        #     valmin=0.02, 
                        #     valmax=0.05,
                        #     valinit=win_init, 
                        #     valstep=0.005
                        # )

                        # ## Trigger on
                        # axtrig_on = plt.axes([0.25, 0.1, 0.65, 0.03])
                        # trig_slider = RangeSlider(
                        #     ax=axtrig_on, 
                        #     label='Threshold (%)',
                        #     valmin=0, 
                        #     valmax=1,
                        #     valinit=(thr_off_init, thr_on_init)
                        # )

                        # ## Update sliders
                        # def update(val):
                        #     # Filter
                        #     freqmin_filt = freq_slider.val[0]
                        #     freqmax_filt = freq_slider.val[1]
                            
                        #     trace_trim_copy = trace.copy()
                        #     trace_filtered = trace_trim_copy.filter('bandpass', freqmin=freqmin_filt, freqmax=freqmax_filt)
                        #     trace_line[0].set_ydata(trace_filtered.data) 

                        #     # Kurtosis
                        #     new_win = win_slider.val
                        #     new_kurt_norm = kurtosis_norm(trace_filtered, new_win)
                            
                        #     new_max_kurtosis = np.max(new_kurt_norm)
                            
                        #     kurto_line[0].set_ydata(new_kurt_norm)
                        #     new_win_3decim = "{:.3f}".format(new_win)
                        #     ax2.set_title(f'Kurtosis - Décalage de la fenêtre : {new_win_3decim} secondes')
                            
                        #     # Threshold
                        #     new_thr_off = trig_slider.val[0]
                        #     new_thr_on = trig_slider.val[1]
                        #     thr_off_val = new_thr_off * new_max_kurtosis
                        #     thr_on_val = new_thr_on * new_max_kurtosis
                        #     kurt_on2.set_ydata(thr_on_val)

                        #     new_all_starttimes = starttimes_trigger_by_kurtosis(new_kurt_norm, new_thr_on, new_thr_off)

                        #     xvline = new_all_starttimes * trace_filtered.stats.delta
                        #     delta = trace_filtered.stats.delta
                        #     seg_new = [np.array([[x*delta, ymin2], [x*delta, ymax2]]) for x in new_all_starttimes] 
                        #     vline.set_segments(seg_new)
                            
                        #     plt.draw()
                        
                        # freq_slider.on_changed(update)
                        # win_slider.on_changed(update)
                        # trig_slider.on_changed(update)

                        plt.subplots_adjust(bottom=0.28, wspace=0.4, hspace=0.4)
                        plt.show()


                        
                        # # Add this matrix into trainwave dictionary
                        # # try:
                        # #     trainwave_data["starttime_kurtosis"][station] = matrix_kurtosis
                        # # except KeyError:
                        # #     trainwave_data["starttime_kurtosis"] = {station : matrix_kurtosis}

                        # # print(trainwave_data)

                        # # if format_save == 'JSON':
                        # #     with open(trainwave_path, 'w') as content:
                        # #         json.dump(trainwave_data, content,  indent=1)       

