""" Module for waveform processing functions """

def standard_process_waveforms(waveforms, apodisation, physical_quantity, pre_filt):
    """
    pre_filt: frequency f1,f2,f3,f4 to characterize the band-pass filter.
        (see https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html)

    """
    waveforms.detrend()
    waveforms.taper(apodisation)
    waveforms.remove_response(output=physical_quantity, pre_filt=pre_filt)
    return waveforms
