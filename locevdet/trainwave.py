""" Class for Trainwave """

import numpy as np

from copy import copy
from typing import Tuple

from obspy import UTCDateTime

from obspy.signal.trigger import trigger_onset
from scipy.signal import hilbert

from locevdet.stations import Station

from locevdet.waveform_processing import trim_trace
from locevdet.utils import rolling_max, kurtosis_norm

class Trainwave():

    def __init__(self, trace, station:Station, start_global, **kwargs):
        self.trace = trace
        self.station = station

        self.start_global = UTCDateTime(start_global)
        self.nsta_time = kwargs.get('nsta_time', 1)
        self.nlta_time = kwargs.get('nlta_time', 30)

        self.pre_trigger = kwargs.get('pre_trigger', 10)
        self.post_trigger = kwargs.get('pre_trigger', 50)
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
        self.noise_level = kwargs.get('noise_level', None)
        self.all_endtimes_delta = kwargs.get('all_endtimes_delta', None)
        self.end_specific = kwargs.get('end_specific', None)

        # Matlab other variables
        self.matlab_data = kwargs.get('matlab_data', None)

    def kurtosis(self, window, threshold_on, threshold_off=0.25):
        kurt_norm = kurtosis_norm(self.trace_filtered, window)
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

    def envelope(self, rolling_max_window:float=0):
        analytic_signal = hilbert(self.trace_filtered)
        envelope = np.abs(analytic_signal)
        print("hilbert envelope :", envelope.size)
        if rolling_max_window > 0:
            envelope = rolling_max(envelope, rolling_max_window)
        return envelope
    
    def snr_calculation(self, rolling_max_window:float=0, 
            time_intervall_inspect:Tuple[float]=[1,10], 
            window_inspect:float=5):

        envelope = self.envelope(rolling_max_window)

        delta = self.trace_filtered.stats.delta
        index_start_global = (self.start_global - self.trace_filtered.stats.starttime) / delta
        window_inspect_npts = window_inspect / delta
        npts_intervall_inspect_noise = time_intervall_inspect[0] / delta
        npts_intervall_inspect_signal = time_intervall_inspect[1] / delta

        end_index_noise_window = int(index_start_global - npts_intervall_inspect_noise)
        start_index_noise_window = int(end_index_noise_window - window_inspect_npts)

        noise_level = np.mean(envelope[start_index_noise_window : end_index_noise_window])
        self.noise_level = noise_level

        start_index_signal_window = int(index_start_global + npts_intervall_inspect_signal)
        end_index_signal_window = int(start_index_signal_window + window_inspect_npts)

        signal_level = np.mean(envelope[start_index_signal_window : end_index_signal_window])

        snr = signal_level / noise_level
        self.snr = snr

        return snr

    def endtime_detection(self, 
            rolling_max_window:float=0, 
            time_intervall_inspect:Tuple[float]=[1,10],
            thr_snr_purcent:float=1.1):
        
        envelope = self.envelope(rolling_max_window)
        print("envelope :", envelope)
        print("envelope size :",envelope.size)
        delta = self.trace_filtered.stats.delta
        index_start_global = int((self.start_global - self.trace_filtered.stats.starttime) / delta)
        index_inspect_signal = int(time_intervall_inspect[1] / delta)

        threshold_snr_end = thr_snr_purcent * self.noise_level
        thrsedhold_on = envelope[index_start_global + index_inspect_signal]
        print("thrsedhold_on :", thrsedhold_on)

        triggersnr_samples_detection = trigger_onset(envelope, thrsedhold_on, threshold_snr_end)

        all_endtimes = triggersnr_samples_detection[:,1]
        all_endtimes_delta = all_endtimes * self.trace_filtered.stats.delta
        all_endtimes_utc = [
            self.trace_filtered.stats.starttime + ends
            for ends in all_endtimes_delta
        ]
        self.all_endtimes_delta = all_endtimes_delta
        print(all_endtimes_utc)

        end_specific = UTCDateTime(all_endtimes_utc[0])
        print(end_specific)
        self.end_specific = end_specific

        return all_endtimes_utc, end_specific

    def __repr__(self):
        return f"Trainwave{self.trace}"
    
