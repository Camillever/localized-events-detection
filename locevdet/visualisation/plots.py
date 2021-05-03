""" Module for plotting waveforms """

import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0}) # To avoid warning from matplotlib mermory

from typing import List

from matplotlib.widgets import RangeSlider, Slider

class RangeSlider(RangeSlider):

    def __init__(self, ax, label, valmin, valmax, **kwargs):
        self.val = (valmin, valmax)
        super().__init__(ax, label, valmin, valmax, **kwargs)

from obspy import read
from tqdm import tqdm

from obspy import UTCDateTime

from locevdet.utils import clean_utc_str
# from locevdet.trainwaves_arrivals.kurtosis import kurtosis_norm, starttimes_trigger_by_kurtosis

# def kurt_param_sliders_per_trace(trace, start_global:UTCDateTime, kurt_norm, 
#     all_starttimes_init:list, 
#     thr_off_init, thr_on_init, win_init, 
#     numero_col_subplot, nb_stations):
    
#     start_name_trace = clean_utc_str(trace.stats.starttime)
#     end_name_trace = clean_utc_str(trace.stats.endtime)

#     max_kurtosis = np.max(kurt_norm)

#     # fig = plt.figure(f"{trace.stats.station} : {start_name_trace} - {end_name_trace}")
#     # fig.suptitle(f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}")

#     # ax1 = fig.add_subplot(2,nb_stations, numero_col_subplot)
#     ax1 = plt.subplot(211)
#     ax1.set_title(f'Signal en fonction du temps - {trace.stats.station} ')
#     ax1.plot(trace.times(), trace.data, color='grey')
#     trace_line = ax1.plot(trace.times(), trace.data)
#     ax1.set_xlabel('temps (secondes)')
    
#     # Start global of trainwave
#     ymin1, ymax1 = ax1.get_ylim()
#     offset = start_global-trace.stats.starttime
#     ax1.axvline(offset, ymin1, ymax1, color='brown')  # Is not working ??
    
#     # ax2 = fig.add_subplot(2,nb_stations, nb_stations + numero_col_subplot)
#     ax2 = plt.subplot(212, sharex=ax1)
#     ax2.set_title(f'Kurtosis - Décalage de la fenêtre : {win_init} secondes')
#     kurto_line = ax2.plot(trace.times(), kurt_norm)
#     ax2.set_ylim(top=1, bottom=0)
#     xmin, xmax = ax2.get_xlim()

#     # Start global of trainwave
#     ymin2, ymax2 = ax2.get_ylim()
#     ax2.axvline(offset, ymin2, ymax2, color='brown')

#     # Kurtosis threshold and start times triggered
#     vline = ax2.vlines(all_starttimes_init*trace.stats.delta, ymin2, ymax2, color='red')
#     kurt_on2 = ax2.axhline(thr_on_init*max_kurtosis, xmin, xmax, color='red',linestyle='--')

#     # SLIDERS

#     ## Filter 
#     axfreq = plt.axes([0.25, 0.2, 0.65, 0.03])
#     freq_slider = RangeSlider(
#         ax=axfreq, 
#         label='Frequency (Hz)',
#         valmin=0.05, 
#         valmax=49.5,
#         valinit=(0.05, 49.5)
#     )
    
#     ## Window for kurtosis
#     axwin = plt.axes([0.25, 0.15, 0.65, 0.03])
#     win_slider = Slider(
#         ax=axwin, 
#         label='Kurtosis : Window shift (s)',
#         valmin=0.02, 
#         valmax=0.05,
#         valinit=win_init, 
#         valstep=0.005
#     )

#     ## Trigger on
#     axtrig_on = plt.axes([0.25, 0.1, 0.65, 0.03])
#     trig_slider = RangeSlider(
#         ax=axtrig_on, 
#         label='Threshold (%)',
#         valmin=0, 
#         valmax=1,
#         valinit=(thr_off_init, thr_on_init)
#     )

#     ## Update sliders
#     def update(val):
#         # Filter
#         freqmin_filt = freq_slider.val[0]
#         freqmax_filt = freq_slider.val[1]
        
#         trace_trim_copy = trace.copy()
#         trace_filtered = trace_trim_copy.filter('bandpass', freqmin=freqmin_filt, freqmax=freqmax_filt)
#         trace_line[0].set_ydata(trace_filtered.data) 

#         # Kurtosis
#         new_win = win_slider.val
#         new_kurt_norm = kurtosis_norm(trace_filtered, new_win)
        
#         new_max_kurtosis = np.max(new_kurt_norm)
        
#         kurto_line[0].set_ydata(new_kurt_norm)
#         new_win_3decim = "{:.3f}".format(new_win)
#         ax2.set_title(f'Kurtosis - Décalage de la fenêtre : {new_win_3decim} secondes')
        
#         # Threshold
#         new_thr_off = trig_slider.val[0]
#         new_thr_on = trig_slider.val[1]
#         thr_off_val = new_thr_off * new_max_kurtosis
#         thr_on_val = new_thr_on * new_max_kurtosis
#         kurt_on2.set_ydata(thr_on_val)

#         new_all_starttimes = starttimes_trigger_by_kurtosis(new_kurt_norm, new_thr_on, new_thr_off)

#         xvline = new_all_starttimes * trace_filtered.stats.delta
#         delta = trace_filtered.stats.delta
#         seg_new = [np.array([[x*delta, ymin2], [x*delta, ymax2]]) for x in new_all_starttimes] 
#         vline.set_segments(seg_new)
        
#         plt.draw()
    
#     freq_slider.on_changed(update)
#     win_slider.on_changed(update)
#     trig_slider.on_changed(update)

#     plt.subplots_adjust(bottom=0.28, wspace=0.4, hspace=0.4)
#     plt.show()



# def capture_one_event(folder_in, folder_out, start_time):
#     """
#     Returns a PNG image of one event started from the given start_time from all seismograms (MSEED) of this event into the given folder
    
#     Args :
#         folder_in : directory path of downloaded and processed MSEED files (str)
#         folder_out : directory path of folder where the png image will be saved
#         start_time : UTC start-time of the record of the event (cf nomenclatures of MSEED seismograms)
#     Returns :
#         A PNG image with all seismograms of this event with the nomenclature : start-time.png
        
#     """
#     all_seismogram = [
#         filename for filename in os.listdir(folder_in)
#         if filename.startswith("PF_") or filename.startswith("G_RER")
#     ]
#     start_time_name = str(start_time)[:-8]
    
#     fig = plt.figure(start_time_name, figsize = [25, 15])
#     filenames_event = []
#     for filename in all_seismogram :
        
#         fig.suptitle(filename.split('_')[2]+' --- '+filename.split('_')[3])
        
#         if start_time_name.replace(':', '-') == filename.split('_')[2]:
#             filepath = os.path.join(folder_in, filename)
#             filenames_event.append(filename)
#             seismogram = read(filepath)
            
#             if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_PER"):
#                 trace = seismogram[0]
                
#             elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
#                 trace = seismogram[0]
#                 trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
                
#             else:
#                 trace = seismogram[2]
            
            
#             ax = fig.add_subplot(len(stations), 1, len(filenames_event))
#             ax.plot(trace.times(), trace.data)
#             ax.set_title(filename.split('_')[0]+'_'+filename.split('_')[1], loc='right', fontweight ="bold")
            
#             plt.tight_layout()
    
#     fig.savefig(os.path.join(folder_out, start_time_name.replace(':', '-')+'.png'))
    
#     plt.show()
    
#     return (print('Figure saved into :', folder_out, 'with the name :', start_time_name.replace(':', '-')+'.png'))




# # Par station
# def captures_per_station(folder_in:str, folder_out:str, all_seismogram, networks_stations, n_cols, n_lines):
#     """
    
    
#     """

#     for network_station in networks_stations :
#         i=0
#         for seismogram_number, filename in enumerate(all_seismogram):
#             if filename.startswith(network_station):
#                 i+=1

#                 filepath = os.path.join(folder, filename)
#                 seismogram = read(filepath)


#                 if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_PER"):
#                     trace = seismogram[0] # Une seule composante et elle est verticale
#                 elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
#                     trace = seismogram[0]
#                     trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
#                 else:
#                     trace = seismogram[2] # Composant verticale

#                 dt = trace.stats.starttime.timestamp
#                 t = np.linspace(trace.stats.starttime.timestamp - dt,
#                                 trace.stats.endtime.timestamp - dt,
#                                 trace.stats.npts)

#                 fig = plt.figure(network_station, figsize = [30, 20])
#                 plt.title(network_station)
#                 plt.subplot(n_lines, n_cols, i)
#                 plt.plot(trace.times(), trace.data)
#                 plt.tight_layout()

#         fig.savefig(os.path.join(folder_out, network_station+'_processed.png'))


# n_cols = 8  # for 72 events
# n_lines = 9

# networks_stations = ['PF_FRE', 'PF_HIM', 'PF_NSR', 'PF_PER', 'PF_TTR', 'G_RER']
# folder_path = os.path.join(path, 'mseed_RESIF')
# folder_out_per_station = os.path.join(folder_path, 'captures','per_station')

# all_seismogram = [
#         filename for filename in os.listdir(folder_path)
#         if filename.startswith("PF_") or filename.startswith("G_RER")
#     ]

# # captures_per_station(folder_path, folder_out_per_station, all_seismogram, networks_stations, n_cols, n_lines)


# #Par événement 

# def get_starttime(filename: str)-> str:
#     post_name = filename.split('_')[2]
#     return post_name



# def captures_per_event(folder, folder_out, all_seismogram, networks_stations):
#     starttime_filename = []

#     HIM_filenames = [
#         filename for filename in os.listdir(folder)
#         if filename.startswith('PF_HIM') 
#     ]

#     for starttime in HIM_filenames:
#         starttime_filename.append(get_starttime(starttime))


#     for start_name in tqdm(starttime_filename):
#         filenames_event = []
#         fig = plt.figure(start_name, figsize = [25, 15])

#         for filename in all_seismogram :

#             if start_name == filename.split('_')[2]:
#                 fig.suptitle(filename.split('_')[2]+' --- '+filename.split('_')[3]) 
#                 filepath = os.path.join(folder, filename)
#                 filenames_event.append(filename)

#                 seismogram = read(filepath)


#                 if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_TTR") or filename.startswith("PF_PER"):
#                     trace = seismogram[0]
#                 elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
#                     trace = seismogram[0]
#                     trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
#                 else:
#                     trace = seismogram[2]

#                 dt = trace.stats.starttime.timestamp
#                 t = np.linspace(trace.stats.starttime.timestamp - dt,
#                                         trace.stats.endtime.timestamp - dt,
#                                         trace.stats.npts)

#                 ax = fig.add_subplot(len(networks_stations), 1, len(filenames_event))

#                 ax.plot(trace.times(), trace.data)
#                 ax.set_title(filename.split('_')[0]+'_'+filename.split('_')[1], loc='right', fontweight ="bold")
#                 plt.tight_layout()

#         fig.savefig(os.path.join(folder_out, start_name+'_processed.png'))
# #         plt.show()

# networks_stations = ['PF_FRE', 'PF_HIM', 'PF_NSR', 'PF_PER', 'PF_TTR', 'G_RER']
# folder_path = os.path.join(path, 'mseed_RESIF')
# folder_out_per_event = os.path.join(folder_path, 'captures', 'per_event')

# all_seismogram = [
#         filename for filename in os.listdir(folder_path)
#         if filename.startswith("PF_") or filename.startswith("G_RER")
#     ]

# # captures_per_event(folder_path, folder_out_per_event, all_seismogram, networks_stations)


# # Par événement (et avec spectogrammes)
# import math
# from scipy import signal

# def nextpow2(i):
#     """
#     Find 2^n that is equal to or greater than.
#     """
    
#     n = 1
#     while n < i:
#         n *= 2
        
#     return n

# def captures_per_event_spectogramms(folder, folder_out, all_seismogram, networks_stations):
    
#     starttime_filename = []
#     HIM_filenames = [
#         filename for filename in os.listdir(folder)
#         if filename.startswith('PF_HIM') 
#     ]

#     for starttime in HIM_filenames:
#         starttime_filename.append(get_starttime(starttime))
    
    
#     for start_name in tqdm(starttime_filename):
#         filenames_event = []
#         fig = plt.figure(start_name, figsize = [25, 15])

#         for filename in all_seismogram :

#             if start_name == filename.split('_')[2]:
#                 fig.suptitle(filename.split('_')[2]+' --- '+filename.split('_')[3]) 

#                 filepath = os.path.join(folder, filename)
#                 filenames_event.append(filename)

#                 seismogram = read(filepath)

#                 if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_TTR") or filename.startswith("PF_PER"):
#                     trace = seismogram[0]
#                 elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
#                     trace = seismogram[0]
#                     trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
#                 else:
#                     trace = seismogram[2]

#                 dt = trace.stats.starttime.timestamp
#                 t = np.linspace(trace.stats.starttime.timestamp - dt,
#                                         trace.stats.endtime.timestamp - dt,
#                                         trace.stats.npts)

#                 # Seismograms
#                 ax1 = fig.add_subplot(len(networks_stations), 2, len(filenames_event))


#                 ax1.plot(trace.times(), trace.data, 'black')

#                 ax1.set_title(filename.split('_')[0]+'_'+filename.split('_')[1], loc='right', fontweight ="bold")
#                 ax1.set_xlabel('Temps (s)')
#                 ax1.set_ylabel('Vitesse (m/s)')

#                 # Spectograms
#                 window_nb_ech = int(5*trace.stats.sampling_rate)
#                 noverlap = int(math.floor(window_nb_ech*90/100))
#                 nfft = int(nextpow2(window_nb_ech))
#                 signal_windowed = np.ones(int(nfft))
#                 window = signal_windowed *signal.tukey(len(signal_windowed), 0.05)
#                 number_pts_padded = int(nfft)

#                 filenames_event.append('spectro')
#                 print(len(filenames_event))
#                 ax2 = fig.add_subplot(len(networks_stations), 2, len(filenames_event))

#                 f, t, sxx = signal.spectrogram(trace.data, trace.stats.sampling_rate, window=window, noverlap=noverlap, nfft=nfft, mode='psd')
#                 sxx_log10_normalized = (np.log10(sxx) - np.min(np.log10(sxx)))/(np.max(np.log10(sxx))-np.min(np.log10(sxx)))
#                 spectro = plt.pcolor(t, f, sxx_log10_normalized>0.9 , cmap='magma')


#                 ax2.set_xlabel('Temps (s)')
#                 ax2.set_ylabel('Fréquence (Hz)')
#                 if filename.startswith("G_RER") :
#                     ax2.set_ylim([0.5, 10])
#                 else : 
#                     ax2.set_ylim([0.5, 25])

#                 fig.colorbar(spectro, ax=ax2)
#                 plt.tight_layout()


#         fig.savefig(os.path.join(folder_out, start_name+'_processed.png'))
#         plt.show() 


# folder_path = os.path.join(path, 'mseed_RESIF')
# folder_out_per_event_spectograms = os.path.join(folder_path, 'captures', 'per_event_with_spectrogram', 'max_energy')
# networks_stations = ['PF_FRE', 'PF_HIM', 'PF_NSR', 'PF_PER', 'PF_TTR', 'G_RER']

# all_seismogram = [
#         filename for filename in os.listdir(folder_path)
#         if filename.startswith("PF_") or filename.startswith("G_RER")
#     ]

# # captures_per_event_spectogramms(folder_path, folder_out_per_event_spectograms, all_seismogram, networks_stations)



def demo_con_style(ax, connectionstyle):
    x1, y1 = 0.3, 0.2
    x2, y2 = 0.8, 0.6

    # ax.plot([x1, x2], [y1, y2], ".")
    ax.annotate("",
                xy=(x1, y1), xycoords='data',
                xytext=(x2, y2), textcoords='data',
                )

    ax.text(.05, .95, connectionstyle.replace(",", ",\n"),
            transform=ax.transAxes, ha="left", va="top")

def hist_band_freq(all_fmin:List[float], all_fc:List[float], 
        fmin:float, fc:float, 
        save_fig_path:str=None, show:bool=True):
    """
    Produces two histograms showing the repartition of minmum frequencies and central frequencies
    determined by spectograms Stockwell method on seismograms from 01/02/2020 to 11/02/2020
    and plots values of fmin et fc choosen on these histograms.

    Args :
        fmin : Minimum frequency choosen for the 'freqmin' of bandpass filtering
        fc : Central frequency choosen for the 'freqmax' of bandpass filtering
        save_fig_path : Save path directory (Default None means no save)
    
    Returns :
        Histograms (can be save into the given save path)
    
    """
    fig_hist = plt.figure('Histogrammes')
    if save_fig_path is not None:
        fig_hist.set_size_inches((20, 10), forward=False)
    title = (
        "Répartition des fréquences minimums et centrales \n"
        "sur les trains d\'ondes d\'éboulements déterminées par la méthode de Stockwell \n"
        "(01/02/2020 - 11/02/2020)"
    )
    fig_hist.suptitle(title)

    ax1 = plt.subplot(121)
    ax1.hist(all_fmin, alpha=0.5)
    ax1.set_title('Répartition des fréquences minimums')
    ax1.set_xlabel('Fréquence (Hz)')
    ax1.set_ylabel('Nombre')
    ax1.axvline(fmin, color='darkred')
    _, ymax1 = ax1.get_ylim()
    fmin_format = float(format(fmin, '.2f'))
    ax1.annotate(
        f"{fmin_format}", 
        xy=(fmin, 0), 
        xytext=(fmin-0.15, ymax1-4), 
        color='darkred'
        )
    ax1.grid(True)

    ax2 = plt.subplot(122)
    ax2.hist(all_fc, alpha=0.5)
    ax2.set_title('Répartition des fréquences centrales')
    ax2.set_xlabel('Fréquence (Hz)')
    ax2.axvline(fc, color='darkred')
    _, ymax2 = ax2.get_ylim()
    fc_format = float(format(fc, '.2f'))
    ax2.annotate(
        f"{fc_format}", 
        xy=(fc, 0),
        xytext=(fc-0.5, ymax2-5), 
        color='darkred'
        )
    ax2.grid(True)

    if show is True :
        fig_hist.show()

    if save_fig_path is not None:
        histname = 'hist_fmin_fc_01-02-2020_11-02-2020.png'
        hist_path = os.path.join(save_fig_path, histname)
        fig_hist.savefig(hist_path)

from locevdet.event import Event, EventList
from obspy.core import Stream

def plot_stalta_per_event(eventlist:EventList, save_fig_path:str=None, show:bool=True):
    """ TODO """
    for event in eventlist:
        utc_start_global = event.start_global

        fig_stalta = plt.figure(f"Sta-Lta _{str(utc_start_global)}")
        if save_fig_path is not None:
            fig_stalta.set_size_inches((20, 10), forward=False)
        list_stations = event.stations
        nb_stations = len(list_stations)

        ax = ['ax' + str(nb) for nb in range(nb_stations)]
        num = -1

        for _, trainwave in event.trainwaves.items():
            num +=1

            ax[num] = fig_stalta.add_subplot(nb_stations, 1, num+1)
            ax[num].set_title(f'{trainwave.trace.stats.station}', fontsize=15, fontweight='bold')
            ax[num].plot(trainwave.trace.times(), trainwave.trace.data, color='black', linewidth=1.1)
            
            ymin1, ymax1 = ax[num].get_ylim()
            start_global = utc_start_global - trainwave.trace.stats.starttime
            start_global_line = ax[num].axvline(start_global, color='darkgreen', linewidth=2.5)

        title = (
            f" Début global d\'un événement à : {utc_start_global}"
        )
        fig_stalta.suptitle(title, fontsize=18, fontweight='bold')
        handles = [start_global_line]
        labels = [" Début global de l'événement détecté par STA-LTA "]
        fig_stalta.legend(
                    handles, 
                    labels, 
                    loc='upper right', 
                    fontsize=16, 
                    fancybox=True, 
                    shadow=True, 
                    bbox_to_anchor=(1.1, 1)
                    )
        plt.tight_layout()
        if show is True:
            fig_stalta.show()
        
        if save_fig_path is not None:
            figname = f"event_{clean_utc_str(utc_start_global)}.png"
            fig_save_path = os.path.join(save_fig_path, figname)
            fig_stalta.savefig(fig_save_path, bbox_inches='tight')

def plot_kurtosis_per_event(eventlist:EventList,
        threshold_on:float, window:float, 
        freqmin:float, freqmax:float,
        save_fig_path:str=None, show:bool=True, **kwargs):
    """ TODO """
    for event in eventlist:
        utc_start_global = event.start_global

        fig_kurt = plt.figure(f"Kurtosis _{str(utc_start_global)}")
        if save_fig_path is not None:
            fig_kurt.set_size_inches((20, 10), forward=False)
        list_stations = event.stations
        nb_stations = len(list_stations)

        ax = ['ax' + str(nb) for nb in range(2*nb_stations)]
        line = ['line' + str(nb) for nb in range(2*nb_stations)]
        vlines = ['vline' + str(nb) for nb in range(2*nb_stations)]
        thr_line = ['thr_line' + str(nb) for nb in range(2*nb_stations)]

        for num, (_, trainwave) in enumerate(event.trainwaves.items()):

            trace_trim = trainwave.trace_trimmed
            trace_filtered = trainwave.trace_filtered
            
            start_name_trace = clean_utc_str(trace_filtered.stats.starttime)
            end_name_trace = clean_utc_str(trace_filtered.stats.endtime)

            ## SUBPLOT SEISMOGRAM ##
            ax[num] = fig_kurt.add_subplot(2, nb_stations, num+1)
            ax[num].set_title(f'{trainwave.station.name}', fontsize=15, fontweight='bold')
            sismo_non_filt = ax[num].plot(trace_filtered.times(), trace_trim.data, color='grey', linewidth=1.1)
            line[num] = ax[num].plot(trace_filtered.times(), trace_filtered.data, color='black', linewidth=1.1)
            
            # Start global of trainwave
            ymin1, ymax1 = ax[num].get_ylim()
            start_global_line1 = ax[num].axvline(trainwave.pre_trigger, color='darkgreen', linewidth=3)

            ## SUBPLOT KURTOSIS ##
            ax[nb_stations + num] = fig_kurt.add_subplot(2, nb_stations, nb_stations+num+1, sharex=ax[num])
            ax[nb_stations + num].set_xlabel('temps (secondes)')
            kurt_norm = trainwave.kurtosis_data['kurtosis_matrix']
            line[nb_stations + num] = ax[nb_stations + num].plot(trace_filtered.times(), kurt_norm, linewidth=1.2)
            ax[nb_stations + num].set_ylim(top=1, bottom=0)
            xmin, xmax = ax[nb_stations + num].get_xlim()

            # Start global of trainwave
            start_global_line2 = ax[nb_stations + num].axvline(trainwave.pre_trigger, 0, 1, color='darkgreen', linewidth=3)

            # Kurtosis threshold and start times triggered
            all_starttimes_delta = trainwave.kurtosis_data['all_starttimes_delta']
            vlines[num] = ax[num].vlines(all_starttimes_delta, ymin1, ymax1, color='darkorange', linewidth=2.3)
            vlines[nb_stations + num] = ax[nb_stations + num].vlines(
                all_starttimes_delta, 
                0, 
                1, 
                color='darkorange', 
                linewidth=2.3
            )

            # Choosen start_specific
            start_specific = trainwave.kurtosis_data['start_specific'] - trace_filtered.stats.starttime
            start_specific_line = ax[num].vlines(
                start_specific, 
                ymin1, 
                ymax1, 
                color='darkred', 
                linewidth=2.5
            )
            ax[nb_stations + num].vlines(
                start_specific, 
                0, 
                1, 
                color='darkred', 
                linewidth=2.5
            )

            max_kurtosis = np.max(trainwave.kurtosis_data['kurtosis_matrix'])
            
            thr_line[nb_stations + num] = ax[nb_stations + num].axhline(
                threshold_on * max_kurtosis, 
                xmin, xmax, 
                color='red',
                linestyle='--'
            )

        title = (
            f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}\n {freqmin}-{freqmax}Hz"
            f" / Fenêtre glissante de {window}s / Seuil à {threshold_on} %"
        )
        fig_kurt.suptitle(title, fontsize=18)
        fig_kurt.subplots_adjust(top=5)

        handles = [
            line[0][0], 
            sismo_non_filt[0], 
            line[nb_stations][0], 
            start_global_line2, 
            vlines[0], 
            start_specific_line
            ]

        labels = [
            'Sismogramme filtré',
            'Sismogramme non filtré',
            'Kurtosis',  
            'Début global de l\événement détecté par STA-LTA', 
            'Débuts potentiels du train d\'onde détectés par Kurtosis', 
            'Début du train d\'onde détecté par Kurtosis '
            ]
        fig_kurt.legend(
            handles, 
            labels, 
            loc='upper right', 
            fontsize=10, 
            fancybox=True, 
            shadow=True, 
            bbox_to_anchor=(1.1, 1.1)
            )
        plt.tight_layout()

        if show is True:
            fig_kurt.show()

        if save_fig_path is not None:
            figname = f"event_{clean_utc_str(utc_start_global)}.png"
            fig_save_path = os.path.join(save_fig_path, figname)
            fig_kurt.savefig(fig_save_path, bbox_inches='tight')


def plots_comparaison_starts_HIM_FRE(eventlist:EventList, save_fig_path:str=None, show:bool=True):
    """ TODO """
    list_HIM_starts_specific =[]
    list_HIM_ti = []
    list_HIM_nsr = []

    list_FRE_starts_specific =[]
    list_FRE_ti = []
    list_FRE_nsr = []

    for event in eventlist:
        for _, trainwave in event.trainwaves.items():
            if trainwave.ti is not None:

                if trainwave.station == 'HIM':
                    list_HIM_ti.append(trainwave.ti)
                    list_HIM_starts_specific.append(trainwave.start_specific)
                    list_HIM_nsr.append(trainwave.snr)

                elif trainwave.station == 'FRE':
                    list_FRE_ti.append(trainwave.ti)
                    list_FRE_starts_specific.append(trainwave.start_specific)
                    list_FRE_nsr.append(trainwave.snr)

    fig_comp = plt.figure('01-02-2020_11-02-2020')
    myFmt = DateFormatter('%Y-%m-%d %H:%M:%S')

    if save_fig_path is not None:
            fig_comp.set_size_inches((20, 10), forward=False)
    
    ax1 = fig_comp.add_subplot(121)
    ax1.set_title('HIM', fontsize=15, fontweight='bold')
    him_sc = ax1.scatter(list_HIM_ti, list_HIM_starts_specific, c=list_HIM_nsr)
    ax1.xaxis.set_major_formatter(myFmt)

    plt.colorbar(him_sc)

    ax2 = fig_comp.add_subplot(121)
    ax2.set_title('FRE', fontsize=15, fontweight='bold')
    fre_sc = ax2.scatter(list_FRE_ti, list_FRE_starts_specific, c=list_FRE_nsr)
    ax2.xaxis.set_major_formatter(myFmt)

    fig_comp.colorbar(fre_sc)

    fig_comp.tight_layout()

    if show is True:
        fig_comp.show()
    
    if save_fig_path is not None:
        figname = "compare_ti_01-02-2020_11-02-2020.png"
        fig_save_path = os.path.join(save_fig_path, figname)
        fig_comp.savefig(fig_save_path, bbox_inches='tight')
