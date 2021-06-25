""" Module to calculate band frequency
from data collected by Stockwell spectograms method (matlab_results) """
import os

from typing import List

import numpy as np
import matplotlib.pyplot as plt


from locevdet.eventlist import EventList

def freq_band_interest(eventlist:EventList, save_fig_path:str=None, show:bool=False):
    """ Create a figure of histograms showing the repartition of seismogram's 
    minimum frequencies and central frequencies

    Args:
        eventlist: Class EventList (see locevdet.eventlist)
        save_fig_path : string of the save directory
        show : boolean to show or not the figure in a matplotlib window.
                Default : False
    Returns :
        The mean of minimum frequencies and the maximum of central frequencies
        if save_fig_path given : Save figure
        if show is True : show the figure in a matplotlib window
    """
    all_central_frq = []
    all_fmin = []

    for event in eventlist:
        for _, trainwave in event.trainwaves.items():
            if trainwave.matlab_data is not None :
                if trainwave.matlab_data['trainwave'] is not None: 
                    all_fmin.append(trainwave.matlab_data['trainwave']['fmin'])
                    all_central_frq.append(trainwave.matlab_data['trainwave']['centralfrequency'])

    mean_fmin = np.mean(all_fmin)
    fc = np.max(all_central_frq)

    # Histograms
    if save_fig_path is not None or show is True:
        hist_band_freq(
            all_fmin, 
            all_central_frq, 
            fmin=mean_fmin, 
            fc=fc, 
            save_fig_path=save_fig_path, 
            show=show)
        
    return float(format(mean_fmin, '.2f')), float(format(fc, '.2f'))

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
