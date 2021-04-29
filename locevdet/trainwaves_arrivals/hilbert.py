""" Module to calculate Hilbert envelope and signal to noise ratio from a given signal"""

import os
import scipy.io
import matplotlib.pyplot as plt
import numpy as np

from obspy.core import read
import scipy.signal as signal
from obspy.signal.filter import envelope

from locevdet.utils import ti_matlab_results, ti_to_utcdatetime


def hilbert(x, N=None, axis=-1):
    """
    Compute the analytic signal, using the Hilbert transform.

    The transformation is done along the last axis by default.

    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    N : int, optional
        Number of Fourier components.  Default: ``x.shape[axis]``
    axis : int, optional
        Axis along which to do the transformation.  Default: -1.

    Returns
    -------
    xa : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`

    Notes
    -----
    The analytic signal ``x_a(t)`` of signal ``x(t)`` is:

    .. math:: x_a = F^{-1}(F(x) 2U) = x + i y

    where `F` is the Fourier transform, `U` the unit step function,
    and `y` the Hilbert transform of `x`. [1]_

    In other words, the negative half of the frequency spectrum is zeroed
    out, turning the real-valued signal into a complex signal.  The Hilbert
    transformed signal can be obtained from ``np.imag(hilbert(x))``, and the
    original signal from ``np.real(hilbert(x))``.

    Examples
    ---------
    In this example we use the Hilbert transform to determine the amplitude
    envelope and instantaneous frequency of an amplitude-modulated signal.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import hilbert, chirp

    >>> duration = 1.0
    >>> fs = 400.0
    >>> samples = int(fs*duration)
    >>> t = np.arange(samples) / fs

    We create a chirp of which the frequency increases from 20 Hz to 100 Hz and
    apply an amplitude modulation.

    >>> signal = chirp(t, 20.0, t[-1], 100.0)
    >>> signal *= (1.0 + 0.5 * np.sin(2.0*np.pi*3.0*t) )

    The amplitude envelope is given by magnitude of the analytic signal. The
    instantaneous frequency can be obtained by differentiating the
    instantaneous phase in respect to time. The instantaneous phase corresponds
    to the phase angle of the analytic signal.

    >>> analytic_signal = hilbert(signal)
    >>> amplitude_envelope = np.abs(analytic_signal)
    >>> instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    >>> instantaneous_frequency = (np.diff(instantaneous_phase) /
    ...                            (2.0*np.pi) * fs)

    >>> fig = plt.figure()
    >>> ax0 = fig.add_subplot(211)
    >>> ax0.plot(t, signal, label='signal')
    >>> ax0.plot(t, amplitude_envelope, label='envelope')
    >>> ax0.set_xlabel("time in seconds")
    >>> ax0.legend()
    >>> ax1 = fig.add_subplot(212)
    >>> ax1.plot(t[1:], instantaneous_frequency)
    >>> ax1.set_xlabel("time in seconds")
    >>> ax1.set_ylim(0.0, 120.0)

    References
    ----------
    .. [1] Wikipedia, "Analytic signal".
           https://en.wikipedia.org/wiki/Analytic_signal
    .. [2] Leon Cohen, "Time-Frequency Analysis", 1995. Chapter 2.
    .. [3] Alan V. Oppenheim, Ronald W. Schafer. Discrete-Time Signal
           Processing, Third Edition, 2009. Chapter 12.
           ISBN 13: 978-1292-02572-8

    """
    from scipy import linalg, fft as sp_fft

    x = np.asarray(x)
    if np.iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = sp_fft.fft(x, N, axis=axis)
    h = np.zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    if x.ndim > 1:
        ind = [np.newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[tuple(ind)]
    x = sp_fft.ifft(Xf * h, axis=axis)
    return x

def rooling_max(signal, win_lenght=None):
    if win_lenght is None:
        win_lenght = len(signal) // 20
    
    max_signal = -np.inf * np.ones_like(signal)
    for i in range(len(signal)):
        if len(signal) - i < win_lenght:
            max_signal[i] = np.max(signal[i:])
        else:
            max_signal[i] = np.max(signal[i:i+win_lenght])

    return max_signal


def envelope_by_hilbert(trace_data:np.array):
    """ TODO
    
    See : https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html
    """

    analytic_signal = hilbert(trace_data)
    envelope = np.abs(analytic_signal)
    return envelope

def get_max_envelope(trace_data, threshold=3):
    envelope_hilbert = envelope_by_hilbert(trace_data)
    envelope_max = rooling_max(envelope_hilbert, 50)
    minmax_env = np.min(envelope_max)  # To determine differently ()
    envelope = envelope_max * (rooling_max(envelope_max, 50) >= threshold * minmax_env)
    return envelope


def similarity_ti(matlab_folder_path:str, mseeds_path:str, stations:list):
    """TODO
    """
    all_data_mat = [ 
        mat for mat in os.listdir(matlab_folder_path)
        if mat.endswith('.mat')
        ]
        
    all_seismogram = os.listdir(mseeds_path)
    
    all_ti_from_mat = []

    for matname in all_data_mat:
        print(matname)
        ti = ti_matlab_results(matlab_folder_path, matname)
        print("ti :", ti)
        ti_utc_datetime = ti_to_utcdatetime(matname, ti)
        all_ti_from_mat.append(ti_utc_datetime)
    
    return all_ti_from_mat


