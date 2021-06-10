""" Class for Event"""

import os
import json

from typing import List

import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.backend_bases import MouseButton

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


    def plot(self, type_graph:str='stalta', seismograms_type:str='discontinuous',
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
            if type_graph == 'stalta' and seismograms_type == 'discontinuous':
                trace = trainwave.trace
                times = trace.times()
                freqmin, freqmax = 0.5, 50
                ax1 = fig.add_subplot(2, nb_stations, num+1)
            else :
                trace_trim = trainwave.trace_trimmed
                trace_filtered = trainwave.trace_filtered
                trace = trace_filtered.copy()
                times = trace_filtered.times()

                freqmin = trainwave.freqmin_interest
                freqmax = trainwave.freqmax_interest

                if type_graph == 'envelope':
                    ax1 = fig.add_subplot(nb_stations, 1, num+1)
                
                elif type_graph == 'stalta' and seismograms_type == 'continuous':
                    freqmin, freqmax = 0.5, 50
                    ax1 = fig.add_subplot(2, nb_stations, num+1)
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

                ax2 = fig.add_subplot(2, nb_stations, (num+nb_stations)+1)
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
                envelope = trainwave.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)
                ax1.plot(times, envelope, color='orange', linewidth=2, label='Enveloppe du signal')
                ymin1, ymax1 = ax1.get_ylim()

                snr_title = ("{:.2f}".format(trainwave.snr))
                demo_con_style(ax1, f"snr : {snr_title}", colorborder='black')

                start_specific = trainwave.kurtosis_data['start_specific'] - trace.stats.starttime
                label_startspecific = "Début précis du train d\'onde détecté par Kurtosis"
                ax1.axvline(start_specific, color='darkred', alpha=0.5, linewidth=3, label=label_startspecific)

                if trainwave.end_specific is not None :
                    endspecific = trainwave.end_specific - trace.stats.starttime
                    ax1.axvline(endspecific, color='darkblue', alpha=0.5, linewidth=3, label="Détection de la fin du train d\'onde")

                if trainwave.matlab_data is not None :
                    ti = trainwave.matlab_data['trainwave']['initial_time']
                    ti_delta = ti - trace.stats.starttime
                    label_ti = "Début précis du train d\'onde détecté par Méthode de Stockwell"
                    ax1.axvline(ti_delta, color='pink', linewidth=3, alpha=0.8, label=label_ti)

                thr_snr = kwargs.get('rolling_max_window', 1.1)
                noise_level = trainwave.noise_level
                ax1.axhline(noise_level, xmin1, xmax1,
                color='purple',
                linestyle='--', 
                linewidth=1,
                label="Niveau de bruit"
                )

                signal_level = noise_level*trainwave.snr
                ax1.axhline(signal_level, 
                color='purple',
                linestyle='--', 
                linewidth=1.2, 
                label="Niveau de signal"
                )

                import numpy as np
                delta = trainwave.trace_filtered.stats.delta
                index_start_global = int((trainwave.start_global - trainwave.trace_filtered.stats.starttime) / delta)
                threshold_on = np.mean(envelope[index_start_global - 1 : index_start_global + 1])
                ax1.axhline(threshold_on,
                color='red',
                linestyle='--', 
                linewidth=1.5
                )

                threshold_snr_end = thr_snr * noise_level
                ax1.axhline(threshold_snr_end,
                color='darkblue',
                linestyle='--', 
                linewidth=1.5
                )

        title = (
            f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}\n"
        )
        fig.suptitle(title, fontsize=18)
        fig.subplots_adjust(top=5)

        
        handles_, labels_ = ax1.get_legend_handles_labels()
        if type_graph != 'envelope':
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

    def pick_arrivals_manually(self):
        plt.close("all")

        fig = plt.figure("pick_arrivals_per_event")
        nb_stations = len(self.stations)
        ax = ["ax" + str(i) for i in range(nb_stations)]
        axvline = ["axvline" + str(i) for i in range(nb_stations)]
        for num, (_, trainwave) in enumerate(self.trainwaves.items()):
            freqmin = trainwave.freqmin_interest
            freqmax = trainwave.freqmax_interest

            ax[num] = fig.add_subplot(nb_stations, 1, num+1)
            Artist.set_picker(ax[num], True)
            ax_title = f'{trainwave.station.name} {freqmin}-{freqmax}Hz'
            ax[num].set_title(ax_title, fontsize=15, fontweight='bold')

            trace_filtered = trainwave.trace_filtered
            ax[num].plot(trace_filtered.times(), trace_filtered)

            start_name_event = clean_utc_str(trace_filtered.stats.starttime)
            end_name_event = clean_utc_str(trace_filtered.stats.endtime)

            label_startglobal = "Début global de l\'événement détecté par STA-LTA"
            startglobal = trainwave.start_global - trace_filtered.stats.starttime
            ax[num].axvline(startglobal, color='darkgreen', linewidth=2.5, label=label_startglobal)

            start_specific = trainwave.kurtosis_data['start_specific'] - trace_filtered.stats.starttime
            label_startspecific = "Début précis du train d\'onde détecté par Kurtosis"
            ax[num].axvline(start_specific, color='darkred', alpha=0.5, linewidth=3, label=label_startspecific)
            
            axvline[num] = ax[num].axvline(x=0., color="red", linestyle = 'dashed')

            if trainwave.matlab_data is not None :
                all_ti = trainwave.matlab_data['trainwave']['all_initial_times']
                # for ti in all_ti:
                #     tis_delta = ti - trace_filtered.stats.starttime
                #     ax[num].axvline(tis_delta, color='purple', linewidth=3, alpha=0.8)
                ax[num].vlines(all_ti, color='purple')
                
                all_td = trainwave.matlab_data['trainwave']['all_central_times']
                # for td in all_td:
                #     tds_delta = td - trace_filtered.stats.starttime
                #     ax[num].axvline(tds_delta, color='blue', linewidth=3, alpha=0.8)
                ax[num].vlines(all_td, color='blue')


                if trainwave.matlab_data['trainwave']['initial_time'] != 'None':
                    ti = trainwave.matlab_data['trainwave']['initial_time']
                    ti_delta = ti - trace_filtered.stats.starttime
                    label_ti = "Début précis du train d\'onde détecté par Méthode de Stockwell"
                    ax[num].axvline(ti_delta, color='pink', linewidth=3, alpha=0.8, label=label_ti)


        
        def on_move(event):
            for n, axe in enumerate(ax):
                if axe == event.inaxes:
                    axvline[n].set_data([event.xdata, event.xdata], [0, 1])
                    plt.draw()                    
                
        def on_pick(event):
            if event.button is MouseButton.LEFT:
                for n, axe in enumerate(ax):
                # For infomation, print which axes the click was in
                    if axe == event.inaxes:
                        for num, (_, trainwave) in enumerate(self.trainwaves.items()):
                            if num == n:
                                utc_starttime = trainwave.trace_filtered.stats.starttime
                                print("utc pick :", utc_starttime + event.xdata)
                                print(f"{trainwave.station.name} - {trainwave.start_specific_manual}")
                                trainwave.start_specific_manual = utc_starttime + event.xdata
                                print(f"{trainwave.start_specific_manual}")
        plt.connect('motion_notify_event', on_move)
        fig.canvas.mpl_connect('button_press_event', on_pick)
        
        title = (f"Fenêtre de visualisation : {start_name_event} - {end_name_event}")
        fig.suptitle(title, fontsize=18)
        plt.tight_layout()
        handles_, labels_ = ax[0].get_legend_handles_labels()
        by_label = dict(zip(labels_, handles_))
        fig.legend(
            by_label.values(), 
            by_label.keys(),
            loc='upper right',
            fontsize=10,
            fancybox=True,
            shadow=True,
            bbox_to_anchor=(1.1, 1.1)
            )

        plt.show()
        
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
