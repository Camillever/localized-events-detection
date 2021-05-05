""" Class for Trainwave """

import numpy as np

from copy import copy

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

    def envelope(self, rolling_max_window_purcent:float=0):
        analytic_signal = hilbert(self.trace_filtered)
        envelope = np.abs(analytic_signal)
        if rolling_max_window_purcent > 0:
            rolling_max_window = int(rolling_max_window_purcent * len(self.trace_filtered))
            envelope = rolling_max(envelope, rolling_max_window)
        return envelope

    def __repr__(self):
        return f"Trainwave{self.trace}"
