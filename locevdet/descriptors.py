""" Module of useful functions to determine descriptors of seismograms """

import numpy as np

from scipy.signal import hilbert

from locevdet.waveform_processing import trim_trace
from locevdet.utils import rolling_max

def envelope_fct(trace_type:str, rolling_max_window:float=0, **kwargs):
    """ Calculate the envelope of the given trace
    specifying if the trace is trimmed and filtered or only filtered
    
    Args:
        trace_type : Specify the type of trace that will have an envelope
        rolling_max_window : length of the moving window (in seconds) to calculate the envelope
        freqmin : minimum frequency of the trace's filter (Hz)
        freqmax : maximum frequency of the trace's filter (Hz)
        kwargs : other useful arguments as
                - trace: Trace of one trainwave, component z. ( See Trace from obspy.core.trace)
                - trace_filtered : Trace of one trainwave already filtered, component z.
                        (See Trace from obspy.core.trace)
                - start_global : global start of the event calculated by STA-LTA method (UTCDateTime)
                - pre_trigger : duration before the start_global to cut the trace, in seconds
                        (see trim_trace function from waveform_processing.py )
                - post_trigger : duration after the start_global to cut the trace, in seconds
                        (see trim_trace function from waveform_processing.py )
    Returns :
        The envelope array of the given trace
    """
    if trace_type == 'trimmed_filtered':
        trace_filtered = kwargs.get('trace_filtered', None)

        if trace_filtered is None:
            trace = kwargs.get('trace', None)
            freqmin_filter = kwargs.get('freqmin', 0.5)
            freqmax_filter = kwargs.get('freqmax', 50)
            pre_trigger = kwargs.get('pre_trigger', 10)
            post_trigger = kwargs.get('post_trigger', 50)
            start_global = kwargs.get('start_global', None)
            try : 
                trace_trimmed = trim_trace(trace.copy(), start_global, pre_trigger, post_trigger)
                trace_filtered = trace_trimmed.filter(
                    'bandpass', freqmin=freqmin_filter, freqmax=freqmax_filter)
                
            except TypeError:
                print(" no trace_filtered calculated")
                raise
        if trace_filtered is not None:
            analytic_signal = hilbert(trace_filtered)

    elif trace_type == 'trace_filtered':
        trace = kwargs.get('trace', None)
        freqmin_filter = kwargs.get('freqmin', 0.5)
        freqmax_filter = kwargs.get('freqmax', 50)

        try:
            trace_filtered = trace.filter(
                'bandpass', freqmin=freqmin_filter, freqmax=freqmax_filter)
        except TypeError:
            print(" no trace_filtered calculated")
            raise

        if trace is not None:
            analytic_signal = hilbert(trace_filtered)
    
    envelope = np.abs(analytic_signal)
    if rolling_max_window > 0:
        envelope = rolling_max(envelope, rolling_max_window)
    return envelope

def snr_calculation_fct(envelope_trim_filt:np.ndarray, envelope_filt:np.ndarray):
    """ Calculate the signal to noise ratio from the given envelopes of the filtered trace and
    the trimmed and filtered trace, by calculating the noise and signal level.

    snr = {signal_level} / {noise_level}

    The noise level is estimated as the mean of the whole filtered trace's envelope
    The signal level is estimated as the maximum value of the trimmed and filterd trace's envelope
    
    Args:
        envelope_trim_filt : Array of the envelope of the trimmed and filtered trace
        envelope_filt : Array of the envelope of the filtered trace
        rolling_max_window : length of the moving window (in seconds) to calculate the envelope
    Returns:
        the noise level and the signal to noise ratio (snr)
    """
    noise_level = np.quantile(envelope_filt, 0.5)
    signal_level = np.max(envelope_trim_filt)
    snr = signal_level / noise_level

    return noise_level, snr
