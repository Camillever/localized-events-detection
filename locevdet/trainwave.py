""" Class for Trainwave """

import numpy as np

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

        # if start picked manually
        self.start_specific_manual = kwargs.get('start_specific_manually', None)

        # Envelope, SNR and trainwave's end detection
        self.snr = kwargs.get('snr', None)
        self.noise_level = kwargs.get('noise_level', None)
        self.all_endtimes_delta = kwargs.get('all_endtimes_delta', None)
        self.end_specific = kwargs.get('end_specific', None)
        self.form_ratio = kwargs.get('form_ratio', None)
        self.duration = kwargs.get('duration', None)

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

    def envelope(self, trace_type:str='trimmed_filtered', rolling_max_window:float=0):
        if trace_type == 'trimmed_filtered':
            analytic_signal = hilbert(self.trace_filtered)
        elif trace_type == 'trace_filtered':
            trace = self.trace.copy()
            freqmin = self.freqmin_interest
            freqmax = self.freqmax_interest
            trace_filtered = trace.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
            analytic_signal = hilbert(trace_filtered)
        envelope = np.abs(analytic_signal)
        if rolling_max_window > 0:
            envelope = rolling_max(envelope, rolling_max_window)
        return envelope
    
    def snr_calculation(self, rolling_max_window:float=0):
        """ TODO """
        envelope = self.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)
        envelope_trace = self.envelope(trace_type='trace_filtered', rolling_max_window=rolling_max_window)

        noise_level = np.quantile(envelope_trace, 0.5)
        self.noise_level = noise_level
        signal_level = np.max(envelope)

        snr = signal_level / noise_level
        self.snr = snr

        return snr

    def endtime_detection(self, 
            rolling_max_window:float=0, 
            time_restricted:float=5,
            time_inspect_startglobal:float=1,
            thr_snr_purcent:float=1.1):
        """ TODO """
        envelope = self.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)

        delta = self.trace_filtered.stats.delta
        index_start_global = int((self.start_global - self.trace_filtered.stats.starttime) / delta)
        index_inspect_signal = int(time_inspect_startglobal / delta)
        thrsedhold_on = np.quantile(envelope[index_start_global - index_inspect_signal : index_start_global + index_inspect_signal], 0.9)

        threshold_snr_end = thr_snr_purcent * self.noise_level
        triggersnr_samples_detection = trigger_onset(envelope, thrsedhold_on, threshold_snr_end)

        all_endtimes = triggersnr_samples_detection[:,1]
        all_endtimes_delta = all_endtimes * self.trace_filtered.stats.delta
        all_endtimes_utc = [
            self.trace_filtered.stats.starttime + ends
            for ends in all_endtimes_delta
            if self.trace_filtered.stats.starttime + ends > self.start_global + time_restricted and \
                 self.trace_filtered.stats.starttime + ends > self.kurtosis_data['start_specific'] + time_restricted
        ]
        self.all_endtimes_delta = all_endtimes_delta

        if len(all_endtimes_utc) != 0:
            end_specific = UTCDateTime(all_endtimes_utc[0])
            self.end_specific = end_specific
        else :
            end_specific = None
        return all_endtimes_utc, end_specific
    
    def form_ratio_and_duration(self, rolling_max_window):
        if self.end_specific is not None and self.kurtosis_data['start_specific'] is not None:
            envelope = self.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)
            max_level = np.max(envelope)
            max_level_index = np.where(envelope == np.max(envelope))
            times = self.trace_trimmed.times()
            max_level_time = []

            for index in max_level_index[0]:
                utc_max_level_time = UTCDateTime(self.trace_filtered.stats.starttime + times[index])
                if utc_max_level_time > self.kurtosis_data['start_specific'] \
                    and utc_max_level_time < self.end_specific :
                    max_level_time = [utc_max_level_time]

            if len(max_level_time) != 0:
                duration = self.end_specific - self.kurtosis_data['start_specific']
                start_to_max = self.end_specific - max_level_time[0]
                self.duration = duration
                self.form_ratio = start_to_max / duration
                return start_to_max / duration, duration

    def __repr__(self):
        return f"Trainwave{self.trace}"
    
