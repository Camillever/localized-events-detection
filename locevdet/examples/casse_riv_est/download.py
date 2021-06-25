""" Download TODO: finish docstring """

import os
from obspy import UTCDateTime

from locevdet.waveform_processing import standard_process_waveforms
from locevdet.download import read_events_times_from_csv, download_from_stations
from locevdet.download import get_events_times_from_some_mseeds
from locevdet.utils import ref_duration_to_time_window

CLIENT = 'RESIF'

pre_filt = (0.1, 0.5, 49, 50)
apodisation = 0.05
my_process_waveforms = \
    lambda waveforms: standard_process_waveforms(waveforms, apodisation, "VEL", pre_filt)


##################################################################################################
################################  01-02-2020_11-02-2020  #########################################
# ROCKFALL - DISCONTINUATED

## All magnitudes 
# mseeds_all_magnit_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "rockfall", "all_magnitude")
# HIM_csv_events = os.path.join("events_csv", "HIM_01-02-2020_11-02-2020_MC3_dump_bulletin.csv")
# early_feb_events_times = read_events_times_from_csv(HIM_csv_events, time_offset=(35, 35))
# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "NSR", "TTR", "PER"],
#     events_times=early_feb_events_times,
#     save_path=mseeds_all_magnit_path,
#     process_waveforms=my_process_waveforms
# )

# download_from_stations(
#     CLIENT,
#     network="G",
#     stations=["RER"],
#     events_times=early_feb_events_times,
#     save_path=mseeds_all_magnit_path,
#     process_waveforms=my_process_waveforms,
#     channel="BH*"
# )

# ## Main event (2020-02-06)
# main_event_start_time = UTCDateTime("2020-02-06T18:33:14.88")
# main_event_end_time = UTCDateTime("2020-02-6T18:34:40.00")
# mseeds_main_event_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "rockfall", "main_event")

# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "TTR", "PER"],
#     events_times=[(main_event_start_time, main_event_end_time)],
#     save_path=mseeds_main_event_path,
#     process_waveforms=my_process_waveforms
# )

# download_from_stations(
#     CLIENT,
#     network="G",
#     stations=["RER"],
#     events_times=[(main_event_start_time, main_event_end_time)],
#     save_path=mseeds_main_event_path,
#     process_waveforms=my_process_waveforms
# )

## Average and Strong magnitudes
## Download others mseeds from initial mseeds you have (on same periods)
old_mseeds = os.path.join('seismograms','01-02-2020_11-02-2020','rockfall','average_strong')
events_times = get_events_times_from_some_mseeds(old_mseeds)
save_path_renew = old_mseeds

# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "TTR", "PER", "NSR"],
#     events_times=events_times,
#     save_path=save_path_renew,
#     process_waveforms=my_process_waveforms
# )

# download_from_stations(
#     CLIENT,
#     network="G",
#     stations=["RER"],
#     events_times=events_times,
#     save_path=save_path_renew,
#     process_waveforms=my_process_waveforms,
#     channel="BH*"
# )

# CONTINUOUS
start_time_continuous = UTCDateTime("2020-02-01T00:00:00.00")
end_time_continuous = UTCDateTime("2020-02-11T00:00:00.00")
# save_path_continuous = os.path.join("seismograms", "01-02-2020_11-02-2020", "continuous")
save_path_continuous_usb = os.path.join("D:","seismograms","01-02-2020_11-02-2020","continuous")

events_times_continuous = [ref_duration_to_time_window(
    time_reference=start_time_continuous, 
    duration=3600, 
    time_offset=(0, 100))] 
end_check = events_times_continuous[0][1]
while end_check < end_time_continuous :
    events_times_continuous.append(ref_duration_to_time_window(
    time_reference=events_times_continuous[-1][0] + 3600, 
    duration=3600, 
    time_offset=(200, 100)))
    end_check = events_times_continuous[-1][1]

last_tuple = events_times_continuous[-1]
events_times_continuous.pop(-1)
events_times_continuous.append((last_tuple[0], end_time_continuous))

download_from_stations(
    CLIENT,
    network="PF",
    stations=["FRE", "HIM", "PER"],
    events_times=events_times_continuous,
    save_path=save_path_continuous_usb,
    process_waveforms=my_process_waveforms, 
    pause=2
)

download_from_stations(
    CLIENT,
    network="G",
    stations=["RER"],
    events_times=events_times_continuous,
    save_path=save_path_continuous_usb,
    process_waveforms=my_process_waveforms, 
    channel="BH*",
    pause=2
)

nsr_starttime = UTCDateTime("2020-02-07T10:00:00.00")
event_times_nsr_continous = [
    dates for dates in events_times_continuous
    if dates[0] > nsr_starttime
]

download_from_stations(
    CLIENT,
    network="PF",
    stations=["NSR"],
    events_times=event_times_nsr_continous,
    save_path=save_path_continuous_usb,
    process_waveforms=my_process_waveforms, 
    pause=2
)


# VOLCANO-TECTONICS Events
# csv_sommital_events =  \
#     os.path.join("events_csv", "volcano_tectonic_01-02-2020_11-02-2020_MC3_dump_bulletin.csv")
# volcanotectonic_events_times = read_events_times_from_csv(csv_sommital_events, time_offset=(35, 65))
# mseeds_volcanotectonic_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "volcanotectonic")
# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "NSR", "TTR", "PER"],
#     events_times=volcanotectonic_events_times,
#     save_path=mseeds_volcanotectonic_path,
#     process_waveforms=my_process_waveforms
# )

# download_from_stations(
#     CLIENT,
#     network="G",
#     stations=["RER"],
#     events_times=volcanotectonic_events_times,
#     save_path=mseeds_volcanotectonic_path,
#     process_waveforms=my_process_waveforms,
#     channel="BH*"
# )


# # LOCAL SEISMS events
# csv_locals_events = os.path.join("events_csv", "locals_01-02-2020_11-02-2020_MC3_dump_bulletin.csv")
# locals_events_times = read_events_times_from_csv(csv_locals_events, time_offset=(35, 65))
# mseeds_locals_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "local_seisms")
# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "NSR", "TTR", "PER"],
#     events_times=locals_events_times,
#     save_path=mseeds_locals_path,
#     process_waveforms=my_process_waveforms
# )
