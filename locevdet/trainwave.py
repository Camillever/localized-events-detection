""" Class for Trainwave """

import numpy as np

from obspy import UTCDateTime

from obspy.signal.trigger import trigger_onset
from scipy.signal import hilbert

from locevdet.stations import Station

from locevdet.waveform_processing import trim_trace
from locevdet.utils import rolling_max, kurtosis_norm, starttimes_trigger_by_kurtosis
from locevdet.descriptors import envelope_fct, snr_calculation_fct

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
        self.kurtosis_params = kwargs.get('kurtosis_params')
        self.kurtosis_data = kwargs.get('kurtosis_data')

        # if start picked manually
        self.start_specific_manual = kwargs.get('start_specific_manually')

        # Envelope, SNR and trainwave's end detection
        self.snr = kwargs.get('snr')
        self.noise_level = kwargs.get('noise_level')
        self.all_endtimes_delta = kwargs.get('all_endtimes_delta')
        self.end_specific = kwargs.get('end_specific')
        self.form_ratio = kwargs.get('form_ratio')
        self.duration = kwargs.get('duration')
        self.picks = kwargs.get('picks')

        # Matlab other variables
        self.matlab_data = kwargs.get('matlab_data')

    def kurtosis(self, window, threshold_on, threshold_off=0.25):
        """ Calculate the kurtosis matrix, the start of the event (and all others potential starts)
        for a given trainwave.
        
        Args:
            window : slidding time window in seconds for the kurtosis
            threshold_on : threshold on to detect the start of the event by kurtosis
            threshold_off : threshold off to detect the end of the event by kurtosis

        Returns:
            The dictionary kurtosis_data containing:
                The array of the kurtosis for the given trace (calculated by the function kurt_norm)
                All potential starts of the event from kurtosis matrix in seconds from the 
                    trace's start
                The start of the event 'start_specific' in UTCDateTime, which is the closest of 
                    start_global 
        """
        kurt_norm = kurtosis_norm(self.trace_filtered, window)
        kurtosis_data = None  # To reset
        if len(kurt_norm) != 0:
            all_starttimes = starttimes_trigger_by_kurtosis(kurt_norm, threshold_on, threshold_off)
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
        """ Calculated the trainwave's envelope for a given type of trace (trace_type)
        
        Args:
            trace_type : Type of trace which envelope is calculated
                    either "trimmed_filtered" or "trace_filtered"
            rolling_max_window : slidding window for the calculation of the trace's envelope
        
        Returns:
            The array of the trainwave's envelope
        """
        if trace_type == 'trimmed_filtered':
            trace_filtered = self.trace_filtered
            envelope = envelope_fct(
                trace_type, rolling_max_window=rolling_max_window, trace_filtered=trace_filtered)

        elif trace_type == 'trace_filtered':
            trace = self.trace.copy()
            freqmin = self.freqmin_interest
            freqmax = self.freqmax_interest
            envelope = envelope_fct(
                trace_type, rolling_max_window=rolling_max_window, 
                trace=trace, freqmin=freqmin, freqmax=freqmax)   

        return envelope
    
    def snr_calculation(self, rolling_max_window:float=0):
        """ TODO """
        envelope_trim_filt = self.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)
        envelope_filt = self.envelope(trace_type='trace_filtered', rolling_max_window=rolling_max_window)

        noise_level, snr = snr_calculation_fct(
            envelope_trim_filt=envelope_trim_filt, envelope_filt=envelope_filt)

        self.noise_level = noise_level
        self.snr = snr
        print('snr_calculation : snr ', self.snr)

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

        print("thr_snr_purcent :", thr_snr_purcent)
        print("self.noise_level :", self.noise_level)
        threshold_snr_end = thr_snr_purcent * self.noise_level
        print("threshold_snr_end :", threshold_snr_end)
        triggersnr_samples_detection = trigger_onset(envelope, thrsedhold_on, threshold_snr_end)

        all_endtimes = triggersnr_samples_detection[:,1]
        all_endtimes_delta = all_endtimes * self.trace_filtered.stats.delta
        if self.kurtosis_data is not None:
            all_endtimes_utc = [
                self.trace_filtered.stats.starttime + ends
                for ends in all_endtimes_delta
                if self.trace_filtered.stats.starttime + ends > self.start_global + time_restricted and \
                    self.trace_filtered.stats.starttime + ends > self.kurtosis_data['start_specific'] + time_restricted
            ]
        else:
            all_endtimes_utc = [
                self.trace_filtered.stats.starttime + ends
                for ends in all_endtimes_delta
                if self.trace_filtered.stats.starttime + ends > self.start_global + time_restricted
            ]
        self.all_endtimes_delta = all_endtimes_delta

        if len(all_endtimes_utc) != 0:
            end_specific = UTCDateTime(all_endtimes_utc[0])
            self.end_specific = end_specific
        else :
            end_specific = None
        return all_endtimes_utc, end_specific
    
    def form_ratio_and_duration(self, rolling_max_window):
        """ TODO """
        if self.end_specific is None :
            _, end_specific = self.endtime_detection(rolling_max_window=rolling_max_window)
            self.end_specific = end_specific
        
        envelope = self.envelope(trace_type='trimmed_filtered', rolling_max_window=rolling_max_window)
        max_level = np.max(envelope)
        max_level_index = np.where(envelope == np.max(envelope))
        times = self.trace_trimmed.times()
        max_level_time = []

        for index in max_level_index[0]:
            utc_max_level_time = UTCDateTime(self.trace_filtered.stats.starttime + times[index])
            if utc_max_level_time > self.start_global \
                and utc_max_level_time < self.end_specific :
                max_level_time = [utc_max_level_time]

        if len(max_level_time) != 0:
            duration = self.end_specific - self.start_global
            start_to_max = self.end_specific - max_level_time[0]
            self.duration = duration
            self.form_ratio = start_to_max / duration
            return start_to_max / duration, duration

    def number_picks(self, rolling_max_window):
        """ TODO """
        envelope_trace = self.envelope(trace_type='trace_filtered', rolling_max_window=rolling_max_window)

        max_level = np.max(envelope_trace)
        threshold = max_level * 0.9
        
        number_picks = trigger_onset(envelope_trace, threshold, threshold)

        self.picks = len(number_picks[:,0])

    def __repr__(self):
        return f"Trainwave{self.trace}"
    
