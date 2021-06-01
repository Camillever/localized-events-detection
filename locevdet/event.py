""" Class for Event"""

import os
import json

from typing import List

import matplotlib.pyplot as plt

from obspy.core import Stream
from obspy.signal.trigger import classic_sta_lta

from locevdet.stations import Station
from locevdet.trainwave import Trainwave

from locevdet.utils import clean_utc_str
from locevdet.visualisation.plots import demo_con_style

class Event():

    def __init__(self, start_global, stations:List[Station]):
        self.start_global = start_global
        self.stations = stations
        self.trainwaves = {}

    def add_trainwaves(self, stream:Stream):
        stream_stations = [trace.stats.station for trace in stream]
        for station in self.stations:
            station_index = stream_stations.index(station.name)
            trace = stream[station_index]
            self.trainwaves[station] = Trainwave(trace, station, self.start_global)


    def plot(self, type_graph:str='stalta', 
            save_fig_path:str=None, show:bool=False, **kwargs):
        """
        Plots the given type of graphic

        Args:
            type_graph : Specify the type of graphic
                - 'stalta' : To plot stalta method on traces per event.
                            This method determines a global start of an event 
                            (triggerd by 2 stations in simultaneous)
                - 'kurtosis' : To plot kurtosis method on traces trimmed and filtered per event.
                            This method determines a start a the event per station (trainwave)
                            more precisely.
                - 'envelope' : To plot the envelopes of traces trimmed and filtered per event.
                            The signal to noise ratio is given and the specific end of the trainwave
                            is also calculated.
            
            save_fig_path : Directory path of the folder where the plots will be saved.
                Per default : None --> No save
            show : To show the plots in a matplotlib window if True. Per default : False.
        Returns : 
            The given type plots and save it eventually (if directory path is given)

        """
        plt.close("all")

        fig = plt.figure(f"{type_graph}")

        if show is None:
            show = save_fig_path is None

        if save_fig_path is not None:
            fig.set_size_inches((20, 10), forward=False)

        nb_stations = len(self.stations)
        for num, (_, trainwave) in enumerate(self.trainwaves.items()):
            if type_graph == 'stalta':
                trace = trainwave.trace
                times = trace.times()
                freqmin, freqmax = 0.5, 50
                ax1 = fig.add_subplot(nb_stations, 2, (num*nb_stations)+1)
            else :
                trace_trim = trainwave.trace_trimmed
                trace_filtered = trainwave.trace_filtered
                trace = trace_filtered.copy()
                times = trace_filtered.times()

                freqmin = trainwave.freqmin_interest
                freqmax = trainwave.freqmax_interest

                if type_graph == 'envelope':
                    ax1 = fig.add_subplot(nb_stations, 1, num+1)
                else:
                    ax1 = fig.add_subplot(2, nb_stations, num+1)

            start_name_trace = clean_utc_str(trace.stats.starttime)
            end_name_trace = clean_utc_str(trace.stats.endtime)
            
            ## SUBPLOT SEISMOGRAM ##

            ax1_title = f'{trainwave.station.name} {freqmin}-{freqmax}Hz'
            ax1.set_title(ax1_title, fontsize=15, fontweight='bold')
            xmin1, xmax1 = ax1.get_xlim()
            ymin1, ymax1 = ax1.get_ylim()

            label_rawtrace = "Sismogramme (0.5-50Hz)"
            label_tracetrim = f"Sismogramme filtré ({freqmin}-{freqmax}Hz)"
            label_startglobal = "Début global de l\'événement détecté par STA-LTA"

            startglobal = trainwave.start_global - trace.stats.starttime
            ax1.axvline(startglobal, ymin1, ymax1, color='darkgreen', linewidth=2.5, label=label_startglobal)

            if type_graph == 'stalta':
                ax1.plot(times, trace.data, color='black', alpha=0.3, linewidth=1.1, label=label_rawtrace)

                nsta = int(trainwave.nsta_time*trace.stats.sampling_rate)
                nlta = int(trainwave.nlta_time*trace.stats.sampling_rate)
                cft = classic_sta_lta(trace.data, nsta, nlta)

                ax2 = fig.add_subplot(nb_stations, 2, (num*nb_stations)+2)
                ax2.plot(times, cft, color='blue', alpha=0.3, linewidth=1.1, label='STA-LTA classique')
                ymin2, ymax2 = ax2.get_ylim()
                ax2.axvline(startglobal, ymin2, ymax2, color='darkgreen', linewidth=2.5, label=label_startglobal)
            
            elif type_graph == 'kurtosis':
                trace_trim = trainwave.trace_trimmed
                trim_times = trace_trim.times()
                ax1.plot(trim_times, trace_trim.data, color='grey', alpha=0.3, linewidth=1.1, label=label_rawtrace)
                ax1.plot(times, trace.data, color='black', linewidth=0.5, label=label_tracetrim)
                ymin1, ymax1 = ax1.get_ylim()

                ax2 = fig.add_subplot(2, nb_stations, nb_stations+num+1)
                ymin2, ymax2 = ax2.get_ylim()
                xmin2, xmax2 = ax2.get_xlim()
                ax2.axvline(startglobal, ymin2, ymax2, color='darkgreen', linewidth=2.5)

                kurt_norm = trainwave.kurtosis_data['kurtosis_matrix']
                ax2.plot(times, kurt_norm, color='blue', linewidth=1.2, label="Kurtosis")

                all_starttimes_delta = trainwave.kurtosis_data['all_starttimes_delta']  
                labelstarts = "Tous les débuts de trains d\'ondes détectés par Kurtosis"
                ax1.vlines(all_starttimes_delta, ymin1, ymax1, colors='darkorange', label=labelstarts)
                ax2.vlines(all_starttimes_delta, ymin2, ymax2, colors='darkorange')

                start_specific = trainwave.kurtosis_data['start_specific'] - trace.stats.starttime
                label_startspecific = "Début précis du train d\'onde détecté par Kurtosis"
                ax1.axvline(start_specific, ymin1, ymax1, color='darkred', alpha=0.5, linewidth=3, label=label_startspecific)
                ax2.axvline(start_specific, ymin2, ymax2, color='darkred', alpha=0.5, linewidth=3)

                if trainwave.matlab_data is not None :
                    ti = trainwave.matlab_data['trainwave']['initial_time']
                    ti_delta = ti - trace.stats.starttime
                    label_ti = "Début précis du train d\'onde détecté par Méthode de Stockwell"
                    ax1.axvline(ti_delta, ymin1, ymax1, color='pink', linewidth=3, alpha=0.5, label=label_ti)
                    ax2.axvline(ti_delta, ymin2, ymax2, color='pink', alpha=0.5, linewidth=3)


                window = kwargs.get('window', 150)
                window_in_time = window * trace.stats.delta
                demo_con_style(ax2, f"Fenêtre glissante : {window_in_time} s", colorborder='blue')
                threshold_on = kwargs.get('threshold_on', 0.5)
                ax2.axhline(threshold_on, 
                xmin2, xmax2,
                color='red',
                linestyle='--', 
                linewidth=1.2
                )

            else:
                ax1.plot(times, trace.data, color='grey', linewidth=0.5, label=label_tracetrim)

                rolling_max_window = kwargs.get('rolling_max_window', 1)
                envelope = trainwave.envelope(rolling_max_window)
                ax1.plot(times, envelope, color='darkblue', linewidth=2.5, label='Enveloppe du signal')
                ymin1, ymax1 = ax1.get_ylim()

                snr_title = format(trainwave.snr, '.2e')
                demo_con_style(ax1, f"snr : {snr_title}", colorborder='black')

                endspecific = trainwave.end_specific - trace.stats.starttime
                ax1.axvline(endspecific, ymin1, ymax1, color='darkred', linewidth=3)

                thr_snr_purcent = kwargs.get('rolling_max_window', 1.1)
                noise_level = trainwave.noise_level
                threshold_snr_end = thr_snr_purcent * noise_level
                ax1.axhline(threshold_snr_end, 
                xmin1, xmax1,
                color='red',
                linestyle='--', 
                linewidth=1.2
                )

        title = (
            f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}\n"
        )
        fig.suptitle(title, fontsize=18)
        fig.subplots_adjust(top=5)

        
        handles_, labels_ = ax1.get_legend_handles_labels()
        if ax2:
            handles2_, labels2_ = ax2.get_legend_handles_labels()
            handles_ += handles2_
            labels_ += labels2_
        by_label = dict(zip(labels_, handles_)) # To have no duplicates in the legend
        fig.legend(
            by_label.values(), 
            by_label.keys(),
            loc='upper right',
            fontsize=10,
            fancybox=True,
            shadow=True,
            bbox_to_anchor=(1.1, 1.1)
            )

        plt.tight_layout()

        if show:
            plt.show()

        if save_fig_path is not None:
            figname = f"event_{clean_utc_str(self.start_global)}.png"
            fig_save_path = os.path.join(save_fig_path, figname)
            fig.savefig(fig_save_path, bbox_inches='tight')


    def to_json(self):
        return {
            'start_global': str(self.start_global),
            'stations': str(self.stations)
        }

    def save(self, filepath, format_save, override=False):
        if override or not os.path.isfile(filepath):
            if format_save == 'JSON':
                with open(filepath, 'w') as content:
                    json.dump(self.to_json(), content,  indent=1)

    def __eq__(self, other):
        if isinstance(other, Event):
            return self.start_global == other.start_global
        else:
            raise TypeError("Only two Events can be compared")

    def __str__(self):
        return clean_utc_str(self.start_global)

    def __repr__(self):
        return f"Event({str(self)}, {self.trainwaves})"
