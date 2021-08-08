""" Class for Event"""

import os
import json

import numpy as np

from typing import List

import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.backend_bases import MouseButton
from matplotlib.widgets import RangeSlider, Slider

from obspy.core import Stream, read
from obspy.signal.trigger import classic_sta_lta

from locevdet.stations import Station
from locevdet.trainwave import Trainwave

from locevdet.utils import clean_utc_str, get_info_from_mseedname, kurtosis_norm
from locevdet.waveform_processing import trim_trace
from locevdet.visualisation.plots import demo_con_style

class Event():

    def __init__(self, start_global, stations:List[Station]):
        self.start_global = start_global
        self.stations = stations
        self.trainwaves = {}

    def add_trainwaves_from_stalta(self, stream:Stream):
        stream_stations = [trace.stats.station for trace in stream]
        for station in self.stations:
            station_index = stream_stations.index(station.name)
            trace = stream[station_index]
            self.trainwaves[station] = Trainwave(trace, station, self.start_global)
        
    def add_station_trainwaves(self, mseeds_folder:str, network:str, station:str):
        from locevdet.stations import STATIONS
        for station_class in STATIONS:
            if station_class.name == station:
                self.stations.append(STATIONS[STATIONS.index(station_class)])
                station_class_remind = station_class

        seismograms_station = [
            mseed for mseed in os.listdir(mseeds_folder)
            if mseed.startswith(f"{network}_{station}")
        ]
        start_global = self.start_global

        for mseed in seismograms_station :
            mseed_info = get_info_from_mseedname(mseed)
            if start_global > mseed_info['starttime'] \
                and start_global < mseed_info['endtime'] :

                filepath = os.path.join(mseeds_folder, mseed)
                seismogram = read(filepath)
                for _,seismo in enumerate(seismogram):
                    if seismo.stats.component == 'Z':
                        trace = seismo
                        self.trainwaves[station] = Trainwave(
                            trace, station_class_remind, self.start_global
                        )
    
    def fc_td_clustering(self, save_fig_path:str):
        import chart_studio
        chart_studio.tools.set_credentials_file(username='Camillever', api_key='m3tWSC2sd64HFTeAkcX9')
        from sklearn.cluster import KMeans
        from plotly.plotly import plot, iplot

        plt.close("all")
        # fig_cluster, ax = plt.subplots(1, 2)
        fig_cluster = plt.figure('clustering_fc_td')
        
        if save_fig_path is not None:
            fig_cluster.set_size_inches((20, 10), forward=False)

        stations_order = {
            "HIM": (1, 'o'),
            "FRE": (2, '.'),
            "PER": (3, 's'),
            "NSR": (4, '*'),
            "RER": (5, 'x')
        }

        start_global_str = clean_utc_str(self.start_global)
        
        all_td_fc = []
        for _, (_, trainwave) in enumerate(self.trainwaves.items()):
            # td_fc_list = []
            if trainwave.matlab_data is not None:
                station_criterion = stations_order[trainwave.station.name]
                for n_fc, fc in enumerate(trainwave.matlab_data['trainwave']['all_centralfrequency']):
                    td = trainwave.matlab_data['trainwave']['all_central_times'][n_fc] - \
                        trainwave.trace.stats.starttime
                    # td_fc_list.append([td, fc, station_criterion[0]])
                    all_td_fc.append([td, fc, station_criterion[0]])
                print('all_td_fc :', all_td_fc)
                # ax[0].scatter(np.array(td_fc_list)[0:,0], np.array(td_fc_list)[0:,1], \
                #     marker=station_criterion[1])
        
        if len(all_td_fc )!=0:
            print("np.array(all_td_fc)[0:,0] :", np.array(all_td_fc)[0:,0])
            scatter = dict(
                mode = "markers",
                name = "y",
                type = "scatter3d",
                x = np.array(all_td_fc)[0:,0], y = np.array(all_td_fc)[0:,1], \
                    z = np.array(all_td_fc)[0:,2],
                marker = dict( size=2, color="rgb(23, 190, 207)" )
            )
            clusters = dict(
                alphahull = 7,
                name = "y",
                opacity = 0.1,
                type = "mesh3d",
                x = np.array(all_td_fc)[0:,0], y = np.array(all_td_fc)[0:,1], \
                    z = np.array(all_td_fc)[0:,2]
            )
            layout = dict(
                title = '3d point clustering',
                scene = dict(
                    xaxis = dict( zeroline=False ),
                    yaxis = dict( zeroline=False ),
                    zaxis = dict( zeroline=False ),
                )
            )
            fig = {"data":[scatter, clusters], "type":layout }
            plot(fig, filename='3d point clustering')
            plt.show()

            # if len(all_td_fc)!=0:
            #     all_td_fc_array = np.array(all_td_fc)
            #     print('all_td_fc_array :',all_td_fc_array)
            #     data_set = np.dstack((all_td_fc_array[0:,0],all_td_fc_array[0:,1]))
            #     data_set = data_set[0]
            #     model = KMeans(3).fit(data_set)

            #     for point in data_set:
            #         if model.predict(point.reshape(1,-1)) == [0]:
            #             ax[1].scatter(point[0], point[1], c='b')
            #         elif model.predict(point.reshape(1,-1)) == [1]:
            #             ax[1].scatter(point[0], point[1], c='g')
            #         elif model.predict(point.reshape(1,-1)) == [2]:
            #             ax[1].scatter(point[0], point[1], c='r')
        
            #     for center in model.cluster_centers_:
            #         ax[1].scatter(center[0],center[1], marker='+')
            #     plt.show()

            # if save_fig_path is not None:
            #     figname = f"event_{start_global_str}.png"
            #     fig_save_path = os.path.join(save_fig_path, figname)
            #     fig_cluster.savefig(fig_save_path, bbox_inches='tight')


    def kurt_param_sliders(self, thr_on_init, thr_off_init, win_init):
        """ Tool to adjust kurtosis' parameters manually with sliders.
        The objectif is to have pertinent event's start.
        
        Args:
            thr_on_init : Initial threshold on to detect the start of the event
            thr_off_init : Initial threshold off to detect the end of the event
            win_init : Initial slidding time window in seconds for the kurtosis
        
        Returns :
            Interactive Matplotlib window with sliders of kurtosis' parameters
        """
        plt.close("all")

        nb_stations = len(self.stations)
        fig, ax = plt.subplots(2, nb_stations)

        vlines_seismo = ["vline_seismo" + str(i) for i in range(nb_stations)]
        vlines = ["vline" + str(i) for i in range(nb_stations)]
        kurto_lines = ["kurto_line"+ str(i) for i in range(nb_stations)]
        kurt_on2s = ["kurt_on2"+ str(i) for i in range(nb_stations)]

        for num, (_, trainwave) in enumerate(self.trainwaves.items()):
            kurtosis_data = trainwave.kurtosis(win_init, thr_on_init, thr_off_init)

            ax[num].set_title(f"{trainwave.station.name}")
            label_trace = "Sismogramme (2 -10 Hz)"
            times = trainwave.trace_filtered.times()
            data_trace = trainwave.trace_filtered.data
            ax[num].plot(times, data_trace, color='black', \
                label=label_trace)
            ax[num].set_xlabel('temps (secondes)')

            startglobal = trainwave.start_global - trainwave.trace.stats.starttime
            label_startglobal = "Début global de l\'événement détecté par STA-LTA"
            ymin1, ymax1 = ax[num].get_ylim()
            ymin2, ymax2 = ax[num + nb_stations].get_ylim()
            ax[num].axvline(startglobal, ymin1, ymax1, color='darkgreen', linewidth=2.5, \
                label=label_startglobal)
            ax[num + nb_stations ].axvline(startglobal, ymin2, ymax2, color='darkgreen', \
                linewidth=2.5)
            
            ax[num + nb_stations].set_title(f'Kurtosis - Décalage de la fenêtre : {win_init} s')
            kurt_norm = kurtosis_data["kurtosis_matrix"]
            kurto_lines[num] = ax[num + nb_stations].plot(times, kurt_norm, color='blue', \
                label='Kurtosis')
            ax[num + nb_stations].set_ylim(top=1, bottom=0)

            # Kurtosis threshold and start times triggered
            all_starttimes_init = kurtosis_data["all_starttimes_delta"]
            vlines_seismo[num] = ax[num].vlines(all_starttimes_init, ymin1, ymax1, color='orange')
            vlines[num] = ax[num + nb_stations].vlines(all_starttimes_init, ymin2, ymax2, color='orange')
            max_kurtosis = np.max(kurt_norm)
            xmin2, xmax2 = ax[num + nb_stations].get_xlim()
            kurt_on2s[num] = ax[num + nb_stations].axhline(thr_on_init*max_kurtosis, xmin2, xmax2, \
                color='red',linestyle='--')
            
            start_name_trace = clean_utc_str(trainwave.trace_filtered.stats.starttime)
            end_name_trace = clean_utc_str(trainwave.trace_filtered.stats.endtime)

            title = (f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}\n")
            fig.suptitle(title, fontsize=18)
        
        # SLIDERS
        ## Window for kurtosis
        axwin = plt.axes([0.25, 0.15, 0.65, 0.03])
        win_slider = Slider(
            ax=axwin, 
            label='Kurtosis : Window shift (s)',
            valmin=0.02, 
            valmax=0.05,
            valinit=win_init, 
            valstep=0.005
        )

        ## Trigger on
        axtrig_on = plt.axes([0.25, 0.1, 0.65, 0.03])
        trig_slider = RangeSlider(
            ax=axtrig_on, 
            label='Threshold (%)',
            valmin=0, 
            valmax=1,
            valinit=(thr_off_init, thr_on_init)
        )

        ## Update sliders
        def update(val):
            new_win = win_slider.val

            new_thr_off = trig_slider.val[0]
            new_thr_on = trig_slider.val[1]
            
            for num, (_, trainwave) in enumerate(self.trainwaves.items()):
                # Kurtosis 
                new_kurt_norm = kurtosis_norm(trainwave.trace_filtered, new_win)
                new_max_kurtosis = np.max(new_kurt_norm)
                kurto_lines[num][0].set_ydata(new_kurt_norm)
                new_win_3decim = "{:.3f}".format(new_win)
                ax[num + nb_stations].set_title(f'Kurtosis - \
                    Décalage de la fenêtre : {new_win_3decim} s')
                
                # Threshold
                thr_off_val = new_thr_off * new_max_kurtosis
                thr_on_val = new_thr_on * new_max_kurtosis
                kurt_on2s[num].set_ydata(thr_on_val)

                kurtosis_data = trainwave.kurtosis(new_win, new_thr_on, new_thr_off)
                new_all_starttimes = kurtosis_data["all_starttimes_delta"]
                ymin2, ymax2 = ax[num + nb_stations].get_ylim()
                seg_new1 = [np.array([[x, ymin1], [x, ymax1]]) for x in new_all_starttimes] 
                seg_new2 = [np.array([[x, ymin2], [x, ymax2]]) for x in new_all_starttimes] 
                vlines_seismo[num].set_segments(seg_new1)
                vlines[num].set_segments(seg_new2)

                plt.draw()
        win_slider.on_changed(update)
        trig_slider.on_changed(update)

        plt.subplots_adjust(bottom=0.28, wspace=0.4, hspace=0.4)
        plt.show()

    def plot(self, type_graph:str='stalta', seismograms_type:str='discontinuous',
            save_fig_path:str=None, show:bool=False, override:bool=False, **kwargs):
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
            override: False per default to avoid override existing folder
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
        freqmin_raw, freq_max_raw = 0.5, 50

        for num, (_, trainwave) in enumerate(self.trainwaves.items()):
            if type_graph == 'stalta' :
                freqmin = kwargs.get('freqmin', 0.5)
                freqmax = kwargs.get('freqmax', 50)
                trace = trainwave.trace.copy().filter('bandpass', freqmin=freqmin, freqmax=freqmax)
                times = trace.times()
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
                
                else:
                    ax1 = fig.add_subplot(2, nb_stations, num+1)

            ## SUBPLOT SEISMOGRAM ##

            ax1_title = f'{trainwave.station.name} {freqmin}-{freqmax}Hz'
            ax1.set_title(ax1_title, fontsize=15, fontweight='bold')
            xmin1, xmax1 = ax1.get_xlim()
            ymin1, ymax1 = ax1.get_ylim()

            label_rawtrace = f"Sismogramme ({freqmin} - {freqmax}Hz)"
            label_tracetrim = f"Sismogramme filtré ({freqmin}-{freqmax}Hz)"
            label_startglobal = "Début global de l\'événement détecté par STA-LTA"

            
            if seismograms_type != 'continuous': 
                startglobal = trainwave.start_global - trace.stats.starttime
                ax1.axvline(startglobal, ymin1, ymax1, color='darkgreen', linewidth=2.5, label=label_startglobal)

            if type_graph == 'stalta':
                nsta = int(trainwave.nsta_time*trace.stats.sampling_rate)
                nlta = int(trainwave.nlta_time*trace.stats.sampling_rate)

                if trainwave.snr is None:
                    rolling_max_window = kwargs.get('rolling_max_window', 100)
                    trainwave.snr_calculation(rolling_max_window=rolling_max_window)
                snr_title = ("{:.2f}".format(trainwave.snr))
                demo_con_style(ax1, f"snr : {snr_title}", colorborder='black')

                if seismograms_type == 'continuous':
                    pre_trigger = kwargs.get('pre_trigger', 10)
                    post_trigger = kwargs.get('post_trigger', 50)

                    trace_trim_stalta = trim_trace(trace, trainwave.start_global, pre_trigger+nlta, post_trigger)

                    trace_trim = trim_trace(trace, trainwave.start_global, pre_trigger, post_trigger)
                    times = trace_trim.times()
                    trace = trace_trim
                    trainwave.trace_filtered = trace

                    startglobal = trainwave.start_global - trace.stats.starttime
                    ax1.axvline(startglobal, ymin1, ymax1, color='darkgreen', linewidth=2.5, label=label_startglobal)

                    cft_trim = classic_sta_lta(trace_trim_stalta.data, nsta, nlta)
                    cft = cft_trim[-len(trace_trim):]

                else: 
                    cft = classic_sta_lta(trace.data, nsta, nlta)

                ax1.plot(times, trace.data, color='black', alpha=0.3, linewidth=1.1, label=label_rawtrace)
                ax2 = fig.add_subplot(2, nb_stations, (num+nb_stations)+1)
                ax2.plot(times, cft, color='blue', alpha=0.3, linewidth=1.1, label='STA-LTA classique')
                ymin2, ymax2 = ax2.get_ylim()
                ax2.axvline(startglobal, ymin2, ymax2, color='darkgreen', linewidth=2.5, label=label_startglobal)
                thr_on = kwargs.get('thr_on', 2.9)
                ax2.axhline(thr_on,
                color='red',
                linestyle='--', 
                linewidth=1,
                label="Seuil de détection par STA-LTA"
                )

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
                linewidth=1.2,
                label='Seuil de détection du début de l\'événement'
                )

            else:
                ax1.plot(times, trace.data, color='grey', linewidth=0.5, label=label_tracetrim)

                rolling_max_window = kwargs.get('rolling_max_window', 4)
                envelope = trainwave.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)
                ax1.plot(times, envelope, color='orange', linewidth=2, label='Enveloppe du signal')
                ymin1, ymax1 = ax1.get_ylim()

                snr_title = ("{:.2f}".format(trainwave.snr))
                demo_con_style(ax1, f"snr : {snr_title}", colorborder='black')

                if trainwave.kurtosis_data is not None :
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

                thr_snr = kwargs.get('thr_snr_purcent', 1.1)
                noise_level = trainwave.noise_level
                # ax1.axhline(noise_level, xmin1, xmax1,
                # color='purple',
                # linestyle='--', 
                # linewidth=1,
                # label="Niveau de bruit"
                # )

                # signal_level = noise_level*trainwave.snr
                # ax1.axhline(signal_level, 
                # color='purple',
                # linestyle='--', 
                # linewidth=1.2, 
                # label="Niveau de signal"
                # )

                import numpy as np
                delta = trainwave.trace_filtered.stats.delta
                index_start_global = int((trainwave.start_global - trainwave.trace_filtered.stats.starttime) / delta)
                threshold_on = np.mean(envelope[index_start_global - 1 : index_start_global + 1])
                ax1.axhline(threshold_on,
                color='red',
                linestyle='--', 
                linewidth=1.5,
                label="Seuil de détection du début d\'un événement"
                )

                threshold_snr_end = thr_snr * noise_level
                ax1.axhline(threshold_snr_end,
                color='darkblue',
                linestyle='--', 
                linewidth=1.5,
                label="Seuil de détection de fin d\'un événement"
                )
            
            start_name_trace = clean_utc_str(trace.stats.starttime)
            end_name_trace = clean_utc_str(trace.stats.endtime)

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
            if override or not os.path.isfile(fig_save_path):
                fig.savefig(fig_save_path, bbox_inches='tight')
        fig.clear()
        plt.close(fig)

    def pick_arrivals_manually(self):
        plt.close("all")

        nb_stations = len(self.stations)
        fig, ax = plt.subplots(nb_stations-1, 1)
        # ax = ["ax" + str(i) for i in range(nb_stations)]
        axvline = ["axvline" + str(i) for i in range(nb_stations)]
        for num, (_, trainwave) in enumerate(self.trainwaves.items()):
            if trainwave.station.name != 'RER':
                freqmin = trainwave.freqmin_interest
                freqmax = trainwave.freqmax_interest

                # ax[num] = fig.add_subplot(nb_stations, 1, num+1)
                Artist.set_picker(ax[num], True)
                ax_title = f'{trainwave.station.name} {freqmin}-{freqmax}Hz'
                ax[num].set_title(ax_title, fontsize=15, fontweight='bold')

                trace_filtered = trainwave.trace_filtered
                    
                ax[num].plot(trace_filtered.times(), trace_filtered)

                start_name_event = clean_utc_str(trace_filtered.stats.starttime)
                end_name_event = clean_utc_str(trace_filtered.stats.endtime)

                label_startglobal = "Début global de l\'événement détecté par STA-LTA"
                startglobal = trainwave.start_global - trace_filtered.stats.starttime
                ax[num].axvline(
                    startglobal, color='darkgreen', linewidth=2.5, label=label_startglobal
                )

                start_specific = trainwave.kurtosis_data['start_specific'] - trace_filtered.stats.starttime
                label_startspecific = "Début précis du train d\'onde détecté par Kurtosis"
                ax[num].axvline(
                    start_specific, 
                    color='darkred', alpha=0.5, 
                    linewidth=3, label=label_startspecific
                )
                
                axvline[num] = ax[num].axvline(x=0., color="red", linestyle = 'dashed')

                if trainwave.matlab_data is not None :
                    all_ti_utc = trainwave.matlab_data['trainwave']['all_initial_times']
                    all_ti = [ti - trace_filtered.stats.starttime for ti in all_ti_utc]
                    ymin, ymax = ax[num].get_ylim()
                    ax[num].vlines(all_ti, ymin, ymax, color='purple')
                    
                    all_td_utc = trainwave.matlab_data['trainwave']['all_central_times']
                    all_td = [td - trace_filtered.stats.starttime for td in all_td_utc]
                    ax[num].vlines(all_td, ymin, ymax, color='blue')

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
