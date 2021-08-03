""" Module for plotting waveforms """

import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0}) # To avoid warning from matplotlib mermory

from typing import List

from matplotlib.widgets import RangeSlider, Slider

# class RangeSlider(RangeSlider):
#     """ Redfine RangeSlider Class to fix initial parameters on slider"""
#     def __init__(self, ax, label, valmin, valmax, **kwargs):
#         self.val = (valmin, valmax)
#         super().__init__(ax, label, valmin, valmax, **kwargs)

from obspy import read
from tqdm import tqdm

from obspy import UTCDateTime
from locevdet.utils import clean_utc_str

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


# Visualisation of seismograms per event
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms
from locevdet.utils import get_info_from_mseedname

def captures_per_event(mseeds_path, save_fig_path:str=None):

    all_seismograms = list_seismograms(mseeds_path)
    stations = list((get_info_from_mseedname(filename))['station']
        for filename in all_seismograms)
    stations = list(dict.fromkeys(stations))
    starts_seismograms = list(str((get_info_from_mseedname(filename))['starttime'])
        for filename in all_seismograms)
    starts_seismograms = list(dict.fromkeys(starts_seismograms))

    for start_period in tqdm(starts_seismograms):
        num = -1
        fig = plt.figure(str(start_period))
        if save_fig_path is not None:
            fig.set_size_inches((20, 10), forward=False)

        for filename in all_seismograms :
            filename_starttime = get_info_from_mseedname(filename)['starttime']
            if filename_starttime == UTCDateTime(start_period) :
                filename_network = get_info_from_mseedname(filename)['network']
                filename_station = get_info_from_mseedname(filename)['station']
                filename_endtime = get_info_from_mseedname(filename)['endtime']
                periodtime = '_'.join((str(filename_starttime), str(filename_endtime)))

                fig.suptitle(periodtime)
                filepath = os.path.join(mseeds_path, filename)
                seismogram = read(filepath)
                for _,seismo in enumerate(seismogram):
                    if seismo.stats.component == 'Z':
                        trace = seismo
                num += 1
                ax = fig.add_subplot(len(stations), 1, num+1)
                ax.plot(trace.times(), trace.data)
                filename_net_station = '_'.join((filename_network, filename_station))
                ax.set_title(filename_net_station, loc='right', fontweight ="bold")
    
        plt.tight_layout()

        if save_fig_path is not None:
            starttime_str = clean_utc_str(periodtime.split('_')[0])
            endtime_str = clean_utc_str(periodtime.split('_')[1])
            periodtime_str = '_'.join((starttime_str, endtime_str))
            figname = f"per_event_{periodtime_str}.png"
            fig_save_path = os.path.join(save_fig_path, figname)
            fig.savefig(fig_save_path, bbox_inches='tight')


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


def demo_con_style(ax, connectionstyle, colorborder):
    x1, y1 = 0.3, 0.2
    x2, y2 = 0.8, 0.6

    # ax.plot([x1, x2], [y1, y2], ".")
    ax.annotate("",
                xy=(x1, y1), xycoords='data',
                xytext=(x2, y2), textcoords='data'
                )

    ax.text(.05, .95, connectionstyle.replace(",", ",\n"),
            transform=ax.transAxes, ha="left", va="top",
            bbox=dict(boxstyle="round", fc="w", color=colorborder))


def circular_hist(ax, x, color_fill, bins=16, density=True, offset=0, gaps=True):
    """
    Produce a circular histogram of angles on ax.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').

    x : array
        Angles to plot, expected in units of radians.
    
    color_fill : RGB color
        Color to fill the bars

    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.

    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.

    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.

    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.

    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.

    bins : array
        The edges of the bins.

    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    # Wrap angles to [-pi, pi)
    x = (x+np.pi) % (2*np.pi) - np.pi

    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)

    # Compute width of each bin
    widths = np.diff(bins)

    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n

    # Plot data on ax
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths,
                     edgecolor=None, fill=True, color=color_fill, alpha=0.5, linewidth=1)

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels for area plots (they are mostly obstructive)
    if density:
        ax.set_yticks([])

    return n, bins, patches

