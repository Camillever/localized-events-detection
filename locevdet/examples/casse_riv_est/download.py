""" Download TODO: finish docstring """

import os
from obspy import UTCDateTime

from locevdet.waveform_processing import standard_process_waveforms
from locevdet.download import read_events_times_from_csv, download_from_stations
from locevdet.download import get_events_times_from_some_mseeds

CLIENT = 'RESIF'

pre_filt = (0.1, 0.5, 49, 50)
my_process_waveforms = \
    lambda waveforms: standard_process_waveforms(waveforms, 0.05, "VEL", pre_filt)



################################  01-02-2020_11-02-2020  #########################################
# ROCKFALL

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

# ## Average and Strong magnitudes
# ## Download others mseeds from initial mseeds you have (on same periods)
# old_mseeds = os.path.join('mseed_RESIF', 'seismograms_PF-HIM_PF-FRE_G-RER')
# events_times = get_events_times_from_some_mseeds(old_mseeds)
# save_path_renew = os.path.join('mseed_RESIF', 'seismograms_all_stations_average_and_strong')

# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "TTR", "PER", "NSR"],
#     events_times=events_times,
#     save_path=save_path_renew,
#     process_waveforms=my_process_waveforms
# )

# VOLCANO-TECTONICS Events
csv_sommital_events =  \
    os.path.join("events_csv", "volcano_tectonic_01-02-2020_11-02-2020_MC3_dump_bulletin.csv")
volcanotectonic_events_times = read_events_times_from_csv(csv_sommital_events, time_offset=(35, 65))
mseeds_volcanotectonic_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "volcanotectonic")
# download_from_stations(
#     CLIENT,
#     network="PF",
#     stations=["FRE", "HIM", "NSR", "TTR", "PER"],
#     events_times=volcanotectonic_events_times,
#     save_path=mseeds_volcanotectonic_path,
#     process_waveforms=my_process_waveforms
# )

download_from_stations(
    CLIENT,
    network="G",
    stations=["RER"],
    events_times=volcanotectonic_events_times,
    save_path=mseeds_volcanotectonic_path,
    process_waveforms=my_process_waveforms,
    channel="BH*"
)


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
