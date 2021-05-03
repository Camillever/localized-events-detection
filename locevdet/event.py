""" Classes Trainwave / Event / Eventlist"""

import os

import json
import pickle

from typing import List
from copy import copy

import numpy as np
import matplotlib.pyplot as plt

from obspy import UTCDateTime
from obspy.core import Stream
from obspy.signal.trigger import trigger_onset

from mat4py import loadmat
from scipy.signal import hilbert

from locevdet.stations import Station
from locevdet.waveform_processing import trim_trace
from locevdet.utils import clean_utc_str, rooling_max, kurtosis_panda

class Trainwave():

    def __init__(self, trace, station:Station, start_global, **kwargs):
        self.trace = trace
        self.station = station
        self.start_global = UTCDateTime(start_global)

        self.pre_trigger = kwargs.get('pre_trigger', 5)
        self.post_trigger = kwargs.get('pre_trigger', 30)
        self.trace_trimmed = trim_trace(self.trace.copy(), self.start_global,
            self.pre_trigger, self.post_trigger)

        self.freqmin_interest = kwargs.get('freqmin_interest', 0.05)
        self.freqmax_interest = kwargs.get('freqmax_interest', 50)
        self.trace_filtered = self.trace_trimmed.copy().filter('bandpass',
            freqmin=self.freqmin_interest, freqmax=self.freqmax_interest)

        # Kurtosis
        self.kurtosis_params = kwargs.get('kurtosis_params', None)
        self.kurtosis_data = kwargs.get('kurtosis_data', None)

        # Envelope, SNR and trainwave's end detection
        self.snr = kwargs.get('snr', None)
        self.end_specific = kwargs.get('end_specific', None)

        # Matlab other variables
        self.matlab_data = kwargs.get('matlab_data', None)

    def kurtosis(self, window, threshold_on, threshold_off=0.25, freqmin=0.05, freqmax=50):
        kurt_norm = kurtosis_panda(self.trace_filtered, window)
        kurtosis_data = None
        if len(kurt_norm) != 0:
            max_kurtosis = np.max(kurt_norm)
            triggertime_trainwaves_kurtosis = trigger_onset(kurt_norm,
                threshold_on*max_kurtosis, threshold_off*max_kurtosis)
            all_starttimes = triggertime_trainwaves_kurtosis[:,0]
            all_starttimes_delta = all_starttimes * self.trace_filtered.stats.delta
            all_starttimes_from_start = [
                self.trace_filtered.stats.starttime + starts
                for starts in all_starttimes_delta
            ]

            specific_start_utc = min(
                all_starttimes_from_start,
                key=lambda x:abs(x-self.start_global)
            )

            kurtosis_data = {
                'start_specific': specific_start_utc,
                'all_starttimes_delta': all_starttimes_delta,
                'kurtosis_matrix': kurt_norm
            }

        self.kurtosis_data = kurtosis_data
        return kurtosis_data

    def envelope(self, rolling_max_window_purcent:float=0):
        analytic_signal = hilbert(self.trace_filtered)
        envelope = np.abs(analytic_signal)
        if rolling_max_window_purcent > 0:
            rolling_max_window = int(rolling_max_window_purcent * len(self.trace_filtered))
            envelope = rooling_max(envelope, rolling_max_window)
        return envelope

    def __repr__(self):
        return f"Trainwave{self.trace}"

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

    def plot_envelope(self, save_fig_path:str=None, show:bool=None):
        plt.close("all")
        fig = plt.figure(f"enveloppe _{str(self.start_global)}")

        if show is None:
            show = save_fig_path is None

        if save_fig_path is not None:
            fig.set_size_inches((20, 10), forward=False)

        nb_stations = len(self.stations)
        for num, (_, trainwave) in enumerate(self.trainwaves.items()):

            trace_trim = trainwave.trace_trimmed
            trace_filtered = trainwave.trace_filtered
            start_name_trace = clean_utc_str(trace_filtered.stats.starttime)
            end_name_trace = clean_utc_str(trace_filtered.stats.endtime)
            freqmin = trainwave.freqmin_interest
            freqmax = trainwave.freqmax_interest

            times = trace_filtered.times()

            ## SUBPLOT SEISMOGRAM ##
            ax = fig.add_subplot(1, nb_stations, num+1)
            ax_title = f'{trainwave.station.name} {freqmin}-{freqmax}Hz'
            ax.set_title(ax_title, fontsize=15, fontweight='bold')
            line_raw = ax.plot(times, trace_trim.data, color='grey', alpha=0.3, linewidth=1.1)
            line_filt = ax.plot(times, trace_filtered.data, color='black', linewidth=0.5)

            # Envelope
            envelope = trainwave.envelope(rolling_max_window_purcent=0.01)
            env_line = ax.plot(times, envelope, color='purple', linewidth=2.5)

            # Start global of trainwave
            start_global_line = ax.axvline(trainwave.pre_trigger, color='darkgreen', linewidth=3)

            # Start specific
            # start_specific = trainwave.kurtosis_data['start_specific'] \
            #     - trace_filtered.stats.starttime
            # start_specific_line = ax.axvline(start_specific, color='darkred', linewidth=3)

            # Threshold SNR
            # thr_line[nb_stations + num] = ax[nb_stations + num].axhline(
            #     threshold_on * max_kurtosis, 
            #     xmin, xmax, 
            #     color='red',
            #     linestyle='--'
            # )

        title = (
            f"Fenêtre de visualisation : {start_name_trace} - {end_name_trace}\n"
        )
        fig.suptitle(title, fontsize=18)
        fig.subplots_adjust(top=5)

        handles = [
            line_raw[0],
            line_filt[0],
            start_global_line,
            # start_specific_line,
            env_line[0],
            ]

        labels = [
            'Sismogramme non filtré',
            'Sismogramme filtré',
            'Début global de l\'événement détecté par STA-LTA',
            # 'Début du train d\'onde détecté par Kurtosis',
            'Enveloppe de chaque train d\'onde'
            ]

        fig.legend(
            handles,
            labels,
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

class EventList(list):

    def __init__(self, events:List[Event]=None):
        events = events if events is not None else []
        super().__init__(events)
        if events is None:
            self.events = []
        else:
            for event in events:
                if not isinstance(event, Event):
                    raise TypeError('EventList must be composed of Event objects')
            self.events = events

    def to_json(self):
        return {str(event): event.to_json() for event in self.events}

    def save(self, filepath, format_save='PICKLE', override=False):
        if override or not os.path.isfile(filepath):
            if format_save == 'JSON':
                with open(filepath, 'w') as content:
                    json.dump(self.to_json(), content,  indent=1)
            elif format_save == 'PICKLE':
                with open(filepath, 'wb') as content:
                    pickle.dump(self, content)


    def set_matlab_data(self, matlab_folder_path, max_time_difference=5):
        matlab_filepaths = [
            os.path.join(matlab_folder_path, filename)
            for filename in os.listdir(matlab_folder_path)
            if filename.endswith('.mat')
        ]

        for matlab_filepath in matlab_filepaths:
            
            title_data = matlab_filepath.split('_')
            
            matlab_data = {
                'network': title_data[3],
                'station': title_data[4],
                'starttime': UTCDateTime(title_data[5]),
                'endtime': UTCDateTime(title_data[6].split('.')[0]),
            }
            matlab_data['periodtime'] = "_".join((str(matlab_data['starttime']), str(matlab_data['endtime'])))
            matlab_starttime = UTCDateTime(matlab_data['starttime'])

            for event in self:
                for _, trainwave in event.trainwaves.items():
                    starttime = UTCDateTime(trainwave.trace.stats.starttime)
                    same_starttime = abs(starttime - matlab_starttime) 
                    if same_starttime < max_time_difference and trainwave.station.name == matlab_data['station']:
                        mat = loadmat(matlab_filepath)
                        matlab_data['trainwave'] = {
                            'radial_signal': mat["wavetrains"]['r'][0],
                            'initial_time': mat["wavetrains"]['Tiwp'][0],
                            'fmin': mat["fmin"],
                            'centralfrequency': mat["wavetrains"]['centralfrequency'][0],
                            'duration': mat["wavetrains"]['duration'][0]
                        }
                        trainwave.matlab_data = matlab_data
                        
