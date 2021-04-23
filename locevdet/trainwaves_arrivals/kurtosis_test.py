""" Module Kurtosis_ test """

import os
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.widgets import RangeSlider, Slider

class RangeSlider(RangeSlider):

    def __init__(self, ax, label, valmin, valmax, **kwargs):
        self.val = (valmin, valmax)
        super().__init__(ax, label, valmin, valmax, **kwargs)

from obspy import UTCDateTime
from obspy.core import read, Stream
from obspy.realtime.signal import kurtosis
from obspy.signal.trigger import trigger_onset

from locevdet.utils import get_starttime_trainwave, get_period, clean_utc_str
from locevdet.waveform_processing import trim_trace
from locevdet.trainwaves_arrivals.hilbert import envelope_by_hilbert

# Utils 

def starttimes_trigger_by_kurtosis(matrix_kurtosis, thr_on, thr_off):
    """ TODO
    
    """
    max_kurtosis = np.max(matrix_kurtosis)
    # print(max_kurtosis)
    triggertime_trainwaves_kurtosis = trigger_onset(matrix_kurtosis, thr_on*max_kurtosis, thr_off*max_kurtosis) 
    # print(triggertime_trainwaves_kurtosis)
    all_starttimes = triggertime_trainwaves_kurtosis[:,0].copy()
    return all_starttimes


def kurtosis_panda(trace, win_kurt_init):
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
    # print("kurt_np :", kurt_np)

    # Normalization
    kurt_max = np.max(kurt_np)
    if kurt_max == 0:
        kurt_norm = []
        return kurt_norm
    kurt_norm = kurt_np/ kurt_max
    # print("kurt_max : ", np.max(kurt_np))

    return kurt_norm

# View all seismograms 
# from locevdet.visualisation.plots import kurt_param_sliders_per_trace

def kurtosistest_for_all_seismograms(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float,
    win:float=3):
    """ TODO
    
    """
    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
        if not filename.startswith('G_RER') 
        if not filename.startswith('PF_TTR')
        
    ]
    

    # Remove specific station on some period
    filenames = []
    all_seismogram_reduced = []
    for filename in all_seismogram:
        starttime = UTCDateTime(filename.split('_')[2])
        station = filename.split('_')[1]
        print(f" Station {station} starts : {starttime}")
        if filename.startswith('PF_NSR') and starttime < UTCDateTime("2020-02-07T10:32:10"):
            filenames.append(filename)
        else:
            all_seismogram_reduced.append(filename)
    
    print(all_seismogram)
    all_trainwaves = [
        trainwave_filename for trainwave_filename in os.listdir(trainwaves_path)
    ]
    # Initial value of thresholds (25th and 75th purcentiles defined by RangeSlider)
    thr_on_init = 0.75
    thr_off_init = 0.25

    for trainwave in all_trainwaves:
        trainwave_path = os.path.join(trainwaves_path, trainwave)
        
        numero = -1

        with open(trainwave_path) as json_file:
            trainwave_data = json.load(json_file)
            start_global = UTCDateTime(trainwave_data['starttime_global'])
            print("start_global :", start_global)
            lists_stations = trainwave_data['order_arrivals_detected'].copy()

            if 'TTR' in lists_stations :
                lists_stations.remove('TTR')
            elif 'RER' in lists_stations :
                lists_stations.remove('RER')
            elif 'NSR' in lists_stations :
                lists_stations.remove('NSR')
            nb_stations = len(lists_stations)

            print(lists_stations)

            ax = ['ax' + str(nb) for nb in range(2*nb_stations)]
            line = ['line' + str(nb) for nb in range(2*nb_stations)]
            vlines = ['vline' + str(nb) for nb in range(2*nb_stations)]
            thr_line = ['thr_line' + str(nb) for nb in range(2*nb_stations)]

            st = Stream()

            for filename in all_seismogram_reduced:
                period = get_period([filename])
                seismo_start = UTCDateTime(period[0].split('_')[0])
                seismo_end = UTCDateTime(period[0].split('_')[1])

                for station in lists_stations:
                    if start_global > seismo_start +100 and start_global < seismo_end-25 and filename.split('_')[1] == station:
                        numero += 1
                        
                        print('numero :', numero, 'and filename :', filename)
                        filepath = os.path.join(seismograms_path, filename)
                        seismogram = read(filepath)
                        fig = plt.figure(trainwave)
                        for _,seismo in enumerate(seismogram):
                            if seismo.stats.component == 'Z':
                                trace = seismo
                        st += trace
                        trace_trim = trim_trace(trace, start_global, pre_trigger, post_trigger)
                        print(trace_trim.stats.npts)
                        
                        # Kurtosis
                        kurt_norm = kurtosis_panda(trace_trim, win)
                            
                  
                        all_starttimes_init = starttimes_trigger_by_kurtosis(kurt_norm, thr_on_init, thr_off_init)

                        start_name_trace = clean_utc_str(trace.stats.starttime)
                        end_name_trace = clean_utc_str(trace.stats.endtime)

                        max_kurtosis = np.max(kurt_norm)


                        ## SUBPLOT SEISMOGRAM ##
                        ax[numero] = fig.add_subplot(2,nb_stations, numero+1)
                        ax[numero].set_title(f'{trace_trim.stats.station}')
                        ax[numero].plot(trace_trim.times(), trace_trim.data, color='grey')
                        line[numero] = ax[numero].plot(trace_trim.times(), trace_trim.data)
                        
                        # Start global of trainwave
                        ymin1, ymax1 = ax[numero].get_ylim()
                        start_global_line1 = ax[numero].axvline(pre_trigger, -ymax1, ymax1, color='darkgreen')  # Is not working ??
                        

                        ## SUBPLOT KURTOSIS ##
                        ax[nb_stations + numero] = fig.add_subplot(2,nb_stations, nb_stations + numero+1, sharex=ax[numero])
                        ax[nb_stations + numero].set_xlabel('temps (secondes)')
                        
                        line[nb_stations + numero] = ax[nb_stations + numero].plot(trace_trim.times(), kurt_norm)
                        ax[nb_stations + numero].set_ylim(top=1, bottom=0)
                        
                        ax[nb_stations + numero].set_ylim([0,max_kurtosis+0.1])
                        xmin, xmax = ax[nb_stations + numero].get_xlim()

                        # Start global of trainwave
                        start_global_line2 = ax[nb_stations + numero].axvline(pre_trigger, 0, 1, color='darkgreen')

                        # Kurtosis threshold and start times triggered
                        all_starttimes_delta = all_starttimes_init*trace_trim.stats.delta

                        vlines[numero] = ax[numero].vlines(all_starttimes_delta, ymin1, ymax1, color='red')
                        vlines[nb_stations + numero] = ax[nb_stations + numero].vlines(all_starttimes_delta, 0, 1, color='red')

                        thr_line[nb_stations + numero] = ax[nb_stations + numero].axhline(thr_on_init*max_kurtosis, xmin, xmax, color='red',linestyle='--')
                

                ## Update sliders
                def update(val):
                    # Values of sliders
                    freqmin_filt = freq_slider.val[0]
                    freqmax_filt = freq_slider.val[1]

                    new_win = win_slider.val

                    new_thr_off = trig_slider.val[0]
                    new_thr_on = trig_slider.val[1]

                    for nb, trace in enumerate(st) :
                        trace_trim = trim_trace(trace, start_global, pre_trigger, post_trigger)
                        trace_trim_copy = trace_trim.copy()
                        trace_filtered = trace_trim_copy.filter('bandpass', freqmin=freqmin_filt, freqmax=freqmax_filt)
                        
                        new_kurt_norm = kurtosis_panda(trace_filtered, new_win)
                         
                        new_max_kurtosis = np.max(new_kurt_norm)

                        thr_off_val = new_thr_off * new_max_kurtosis
                        thr_on_val = new_thr_on * new_max_kurtosis

                        new_all_starttimes = starttimes_trigger_by_kurtosis(new_kurt_norm, new_thr_on, new_thr_off)
                        tr_delta = trace_filtered.stats.delta
                        

                        # #Update of plots
                        line[nb][0].set_ydata(trace_filtered.data)
                        line[nb + nb_stations][0].set_ydata(new_kurt_norm)
                        # new_win_3decim = "{:.3f}".format(new_win)
                        # ax2.set_title(f'Kurtosis - Décalage de la fenêtre : {new_win_3decim} s')

                        ax[nb + nb_stations].set_ylim([0,new_max_kurtosis+0.1])

                        xvline = new_all_starttimes * tr_delta
                        seg_new = [np.array([[x*tr_delta, 0], [x*tr_delta, 1]]) for x in new_all_starttimes] 


                        thr_line[nb + nb_stations].set_ydata(thr_on_val)
                        vlines[nb].set_segments(seg_new)
                        vlines[nb + nb_stations].set_segments(seg_new)
                        

                    fig.canvas.draw_idle()
        # SLIDERS
            if len(st)!=0:
                ## Filter 
                axfreq = fig.add_axes([0.25, 0.2, 0.65, 0.03])
                freq_slider = RangeSlider(
                    ax=axfreq, 
                    label='Frequency (Hz)',
                    valmin=0.05, 
                    valmax=49.5,
                    valinit=(0.05, 49.5), 
                    valstep=0.05
                )
                
                ## Window for kurtosis
                axwin = fig.add_axes([0.25, 0.15, 0.65, 0.03])
                win_slider = Slider(
                    ax=axwin, 
                    label='Kurtosis : Window shift (s)',
                    valmin=4, 
                    valmax=500,
                    valinit=win, 
                    valstep=2
                )

                ## Trigger on
                axtrig_on = fig.add_axes([0.25, 0.1, 0.65, 0.03])
                trig_slider = RangeSlider(
                    ax=axtrig_on, 
                    label='Threshold (%)',
                    valmin=0, 
                    valmax=1,
                    valinit=(thr_off_init, thr_on_init)
                )
                
                if len(st)!=0:
                    freq_slider.on_changed(update)
                    win_slider.on_changed(update)
                    trig_slider.on_changed(update)
                    fig.suptitle(f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}")
                    print(f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}")
                    plt.subplots_adjust(bottom=0.28, wspace=0.4, hspace=0.4)
                    plt.show()

                    print('Done')

###############################################################################################

def kurtosis_starts_extraction(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float,
    thr_on:float, win:float,
    freqmin:float=0.05, freqmax:float=50,
    save_path:str=None, format_save:str='JSON', override:bool=False):

    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
        if not filename.startswith('G_RER') 
        if not filename.startswith('PF_TTR')  
    ]
    
    # Remove specific station on some period
    filenames = []
    all_seismogram_reduced = []
    for filename in all_seismogram:
        starttime = UTCDateTime(filename.split('_')[2])
        station = filename.split('_')[1]
        if filename.startswith('PF_NSR') and starttime < UTCDateTime("2020-02-07T10:32:10"):
            filenames.append(filename)
        else:
            all_seismogram_reduced.append(filename)

    if save_path is None:
        save_path = os.getcwd()
    save_path_format = os.path.join(save_path, format_save)
    save_path_png = os.path.join(save_path, 'PNG')
    os.makedirs(save_path_format, exist_ok=True)
    os.makedirs(save_path_png, exist_ok=True)
    
    all_trainwaves = [
        trainwave_filename for trainwave_filename in os.listdir(trainwaves_path)
    ]
    for trainwave in all_trainwaves:
        trainwave_path = os.path.join(trainwaves_path, trainwave)
        
        numero = -1

        with open(trainwave_path) as json_file:
            trainwave_data = json.load(json_file)
            start_global = UTCDateTime(trainwave_data['starttime_global'])
            print("start_global :", start_global)
            lists_stations = trainwave_data['order_arrivals_detected'].copy()
            print('before :',lists_stations)

            # Copy
            import copy
            new_trainwave = copy.deepcopy(trainwave_data)
            print(new_trainwave)

            lists_stations_reduced =[]
            station_removed = []
            for station in lists_stations:
                if station == 'TTR' :
                    station_removed.append(station)
                elif station == 'RER':
                    station_removed.append(station)
                elif station == 'NSR' and start_global < UTCDateTime("2020-02-07T10:32:10"):
                    station_removed.append(station)
                else : 
                    lists_stations_reduced.append(station)
            nb_stations = len(lists_stations_reduced)

            print('after :', lists_stations_reduced)
            print('stations removed :', station_removed)
    
            ax = ['ax' + str(nb) for nb in range(2*nb_stations)]
            line = ['line' + str(nb) for nb in range(2*nb_stations)]
            vlines = ['vline' + str(nb) for nb in range(2*nb_stations)]
            thr_line = ['thr_line' + str(nb) for nb in range(2*nb_stations)]

            st = Stream()

            for filename in all_seismogram_reduced:
                period = get_period([filename])
                seismo_start = UTCDateTime(period[0].split('_')[0])
                seismo_end = UTCDateTime(period[0].split('_')[1])
                file_station = filename.split('_')[1]

                for station in lists_stations_reduced:
                    if start_global > seismo_start and start_global < seismo_end and file_station == station:
                        lists_stations_reduced.remove(file_station)
                        numero += 1
                        
                        print('numero :', numero, 'and filename :', filename)
                        filepath = os.path.join(seismograms_path, filename)
                        seismogram = read(filepath)
                        fig = plt.figure(trainwave)
                        fig.set_size_inches((20, 10), forward=False)
                    
                        for _,seismo in enumerate(seismogram):
                            if seismo.stats.component == 'Z':
                                trace = seismo
                                
                        
                        trace_trim = trim_trace(trace, start_global, pre_trigger, post_trigger)
                        trace_filtered = trace_trim.filter('bandpass', freqmin=freqmin, freqmax=freqmax)

                        # Hilbert

                        envelope, _ = envelope_by_hilbert(trace_filtered.data)
                        
                        # Kurtosis
                        kurt_norm = kurtosis_panda(trace_filtered, win)
                        if len(kurt_norm) != 0:
                            st += trace
                            thr_off = 0.25           
                            all_starttimes = starttimes_trigger_by_kurtosis(kurt_norm, thr_on, thr_off)
                            print(all_starttimes)

                            start_name_trace = clean_utc_str(trace.stats.starttime)
                            end_name_trace = clean_utc_str(trace.stats.endtime)

                            max_kurtosis = np.max(kurt_norm)

                            # Add features to the new dictionary

                            
                            # all_starttimes_array = np.array(all_starttimes)
                            # pos = (np.abs(all_starttimes_array-myNumber)).argmin()
                            # start_kurt = all_starttimes[pos]
                                
                            
                            # # Add this matrix into trainwave dictionary
                        # # try:
                        # #     new_trainwave["kurtosis"][file_station] = start_kurt
                        # # except KeyError:
                        # #     new_trainwave["kurtosis"] = {file_station : start_kurt}



                            ## SUBPLOT SEISMOGRAM ##
                            ax[numero] = fig.add_subplot(2,nb_stations, numero+1)
                            ax[numero].set_title(f'{trace_filtered.stats.station}', fontsize=15, fontweight='bold')
                            ax[numero].plot(trace_filtered.times(), trace_trim.data, color='grey', linewidth=1.2)
                            ax[numero].plot(trace_filtered.times(), envelope, color='orange', linewidth=2)
                            line[numero] = ax[numero].plot(trace_filtered.times(), trace_filtered.data, linewidth=1.2)
                            
                            # Start global of trainwave
                            ymin1, ymax1 = ax[numero].get_ylim()
                            start_global_line1 = ax[numero].axvline(pre_trigger, -ymax1, ymax1, color='darkgreen')  # Is not working ??
                            

                            ## SUBPLOT KURTOSIS ##
                            ax[nb_stations + numero] = fig.add_subplot(2,nb_stations, nb_stations + numero+1, sharex=ax[numero])
                            ax[nb_stations + numero].set_xlabel('temps (secondes)')
                            
                            line[nb_stations + numero] = ax[nb_stations + numero].plot(trace_filtered.times(), kurt_norm, linewidth=1.2)
                            ax[nb_stations + numero].set_ylim(top=1, bottom=0)
                            
                            xmin, xmax = ax[nb_stations + numero].get_xlim()

                            # Start global of trainwave
                            start_global_line2 = ax[nb_stations + numero].axvline(pre_trigger, 0, 1, color='darkgreen', linewidth=3)

                            # Kurtosis threshold and start times triggered
                            all_starttimes_delta = all_starttimes*trace_filtered.stats.delta

                            vlines[numero] = ax[numero].vlines(all_starttimes_delta, ymin1, ymax1, color='red', linewidth=2.5)
                            vlines[nb_stations + numero] = ax[nb_stations + numero].vlines(all_starttimes_delta, 0, 1, color='red', linewidth=2.5)

                            thr_line[nb_stations + numero] = ax[nb_stations + numero].axhline(thr_on*max_kurtosis, xmin, xmax, color='red',linestyle='--')
            if len(st)!=0:
                title = f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}\n {freqmin}-{freqmax}Hz / Fenêtre glissante de {win}s / Seuil à {thr_on} %"
                fig.suptitle(title, fontsize=20)
                print(f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}")
                # fig.subplots_adjust(bottom=0.28, wspace=0.4, hspace=0.4)
                labels = ['Début global du train d\'onde détecté par STA-LTA', 'Débuts potentiels du train d\'onde détectés par Kurtosis']
                fig.legend([start_global_line2, vlines[0]], labels, loc='upper right', fontsize=11)
                plt.tight_layout()
                # plt.show()

                # Save figures
                figname = trainwave + '.png'
                fig_save_path = os.path.join(save_path_png, figname)
                fig.savefig(fig_save_path, bbox_inches='tight')
            
                # Save specific starts into dictionary
                
                # trainwave_filepath = os.path.join(save_path_format, trainwave)
                # if override or not os.path.isfile(save_path_format):
                # if format_save == 'JSON':
                #     with open(trainwave_filepath, 'w') as content:
                #         json.dump(new_trainwave, content,  indent=1)

                         






