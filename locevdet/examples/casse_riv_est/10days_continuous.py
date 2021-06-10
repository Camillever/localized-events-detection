""" Script to identify events in continuous seismograms """

import os

from obspy import UTCDateTime

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.utils import get_info_from_mseedname
##################################################################################################
################################# 01/02/2020 - 11/02/2020 ########################################
mseeds_path_continuous = os.path.join("seismograms", "01-02-2020_11-02-2020", "continuous")

start_time_continuous = UTCDateTime("2020-02-01T00:00:00.00")
end_time_continuous = UTCDateTime("2020-02-11T00:00:00.00")

################################## STA-LTA #######################################################

config_stalta = {
    # Filter
    "freqmin": 2,
    "freqmax": 10,

    # Sta-lta parameters
    "nsta_time": 1, # in seconds
    "nlta_time": 30, # in seconds
    "thr_on": 2.9,
    "thr_off": 1,
    "thr_coincidence_sum": 2,

    # Remove useless events (duplicate, false)
    "minimum_time" : 2
}

before_07fev_seismograms =[ 
    mseed for mseed in os.listdir(mseeds_path_continuous)
    if get_info_from_mseedname(mseed)['starttime'] < start_time_continuous
    if not mseed.startswith('PF_NSR') or not mseed.startswith('G_RER')
    ]
after_07fev_seismograms =[ 
    mseed for mseed in os.listdir(mseeds_path_continuous)
    if get_info_from_mseedname(mseed)['starttime'] > start_time_continuous
    if not mseed.startswith('G_RER')
    ]
continuous_seismograms = [*before_07fev_seismograms, *after_07fev_seismograms]
print('seismograms :', continuous_seismograms)

events_list_continuous, duplicate_ev_removed, false_ev_removed = stalta_detect_events(
    folder_in=mseeds_path_continuous,
    all_seismogram=continuous_seismograms,
    **config_stalta
)

print(f"Events ({len(events_list_continuous)}): {events_list_continuous}")
print(f"{len(duplicate_ev_removed)} dictionnaires ont été supprimés, dont : {duplicate_ev_removed}")
print(f"{len(false_ev_removed)} dictionnaires ont été supprimés, dont : {false_ev_removed}")


# Plots
## Per event
save_stalta_fig_path = os.path.join("captures", "01-02-2020_11-02-2020", "continuous")
for event in events_list_continuous:
    event.plot(type_graph='stalta', seismograms_type='continuous', save_fig_path=save_stalta_fig_path)
