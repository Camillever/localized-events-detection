""" Module for waveform processing functions """
from obspy import UTCDateTime
from obspy.core import Trace

def standard_process_waveforms(waveforms, apodisation:float, physical_quantity:str, pre_filt):
    """ Process the given waveforms (Stream) by applying :
    Detrend (remove the trend):
        https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.detrend.html
    Apodization on borders :
        https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.taper.html
    Remove instrumental response and Bandpass filter :
        https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html

    Args:
        waveforms : Stream of multiple traces (see Obspy)
        apodisation : Decimal percentage of taper at one end (ranging from 0. to 0.5)
            
        physical_quantity : Specify the output unit of waveforms
                "DISP" : Displacement
                "VEL" : Velocity
                "ACC" : Acceleration
        pre_filt: frequency f1,f2,f3,f4 to characterize the band-pass filter.
            (see https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html)
        
    Returns:
        Waveforms (Stream) processed
    """
    waveforms.detrend()
    waveforms.taper(apodisation)
    waveforms.remove_response(output=physical_quantity, pre_filt=pre_filt)
    return waveforms

def trim_trace(trace:Trace, center_value:UTCDateTime, pre_offset:float=0, post_offset:float=0):
    """  Cut in time the given trace before and after a given central value

    Args:
        trace: Trace (see https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.html)
        center_value : Central value in UTCDateTime
            (Be careful, this value must be integrated into the time band of the trace)
        pre_offset = Time intervall before the center_value which be cut (in seconds)
        post_offset = Time intervall after the center_value which be cut (in seconds)
    
    Returns:
        Trace cut in time
    """
    trace_copy = trace.copy()
    trace_copy.trim(\
        center_value - pre_offset,
        center_value + post_offset,
        nearest_sample=True
    )
    return trace_copy
