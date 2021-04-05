""" Module Kurtosis to detect specific and more precise start times of trainwaves, per station """

import os
import json

import numpy as np

from obspy import UTCDateTime
from obspy.core import read
from obspy.realtime.signal import kurtosis

from locevdet.utils import get_starttime_trainwave, get_period


def kurtosis_for_all_seismograms(seismograms_path:str, trainwaves_path:str, 
    pre_trigger:float, post_trigger:float):
    """ TODO
    
    """

    all_seismogram = [
        filename for filename in os.listdir(seismograms_path)
    ]

    all_trainwaves = [
        trainwave_filename for trainwave_filename in os.listdir(trainwaves_path)
    ]

    for filename in all_seismogram:
        filepath = os.path.join(seismograms_path, filename)
        seismogram = read(filepath)

        for _,seismo in enumerate(seismogram):
            if seismo.stats.component == 'Z':
                trace = seismo
        period = get_period([filename])
        trace_start = UTCDateTime(period[0].split('_')[0])
        trace_end = UTCDateTime(period[0].split('_')[1])
        print(f"trace start: {trace_start} and end: {trace_end}")

        for trainwave in all_trainwaves:
            utc_starttime_trainwave = get_starttime_trainwave(trainwave)
            if utc_starttime_trainwave > trace_start and utc_starttime_trainwave < trace_end:
                trainwave_path = os.path.join(trainwaves_path, trainwave)
                with open(trainwave_path) as json_file:
                    trainwave_data = json.load(json_file)
                    print("Type:", type(trainwave_data) )
                    print("trainwave_data :", trainwave_data)

                    start_global = UTCDateTime(trainwave_data["starttime_global"])
                    print("start_global :", start_global, 'type :', type(start_global))
                    print("trace.stats.starttime :", trace.stats.starttime, 'type :', type(trace.stats.starttime))
                    
                    offset_start_trainwave = np.abs(trace.stats.starttime - start_global)
                    print("offset_start_trainwave :", offset_start_trainwave)

                    pre_offset = offset_start_trainwave - pre_trigger
                    post_offset = offset_start_trainwave + post_trigger
                    print(f"Offsets : {pre_offset} and {post_offset}")
                    trace_time = trace.times().copy()
                    pre_index_time = np.where(trace_time==int(pre_offset))
                    post_index_time = np.where(trace_time==int(post_offset))
                    print(f"index_time : {pre_index_time} and {post_index_time}")

                    trace_window = trace.data[int(pre_index_time[0]) : int(post_index_time[0])]
                    print("trace_window :", trace_window)

        
