""" Module for utilities functions """

from typing import Tuple, List
from obspy import UTCDateTime

def clean_utc_str(utc_datetime:UTCDateTime) -> str:
    """
    To reduce the length of date (in UTCDateTime)
    (Useful for nomenclature/name of file)

    Args:
        utc_datetime: date in UTCDateTime (YYYY-MM-DDThh:mm:ss.xxxx)

    Returns:
        A string of the date as "YYYY-MM-DDThh-mm-ss"

    """
    return str(utc_datetime).split('.')[0].replace(':','-')


def ref_duration_to_time_window(time_reference:UTCDateTime, duration:float,
        time_offset:Tuple[int]=(0, 0)) -> Tuple[float]:
    """
    TODO

    Args:
        time_reference: TODO
        duration: TODO
        time_offset: TODO

    Returns:
        TODO
    """
    start_time = time_reference - time_offset[0]
    end_time = time_reference + duration + time_offset[1]
    return start_time, end_time


def get_period(list_filenames: List[str])-> List[str]:
    """
    Compile all periods of events.

    Args:
        list_filenames: all filenames of MSEED seismograms,
                        with nomenclature : 'network_station_starttime_endtime'.

    Returns:
        List of periods of each events, written : 'starttime_endtime'

    """
    periods_list = []
    for filename in  list_filenames:
        period = '_'.join([filename.split('_')[2],filename.split('_')[3]])
        if period not in periods_list:
            periods_list.append(period)
    return periods_list


def get_starttime_trainwave(filename_trainwave):
    """
    Give the start time string of the trainwave dictionary name

    Args:
        filename_trainwave: file name of the trainwave dictionary

    Returns:
        String of the start time of this trainwave

    """
    return filename_trainwave.split('_')[-1]
