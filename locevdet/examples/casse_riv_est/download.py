""" Download TODO: finish docstring """

import os
from obspy import UTCDateTime

from locevdet.waveform_processing import standard_process_waveforms
from locevdet.download import read_events_times_from_csv, download_from_stations

CLIENT = 'RESIF'

pre_filt = (0.1, 0.5, 49, 50)
my_process_waveforms = \
    lambda waveforms: standard_process_waveforms(waveforms, 0.05, "VEL", pre_filt)

HIM_csv_events = os.path.join("events_csv", "HIM_01-02-2020_11-02-2020_MC3_dump_bulletin.csv")
mseeds_client_path = os.path.join("mseeds", CLIENT)

# Catalog events (2020-02-01 to 2020-02-11)
early_feb_events_times = read_events_times_from_csv(HIM_csv_events, time_offset=(110, 35))
download_from_stations(
    CLIENT,
    network="PF",
    stations=["FRE", "HIM", "NSR", "TTR", "PER"],
    events_times=early_feb_events_times,
    save_path=mseeds_client_path,
    process_waveforms=my_process_waveforms
)

download_from_stations(
    CLIENT,
    network="G",
    stations=["RER"],
    events_times=early_feb_events_times,
    save_path=mseeds_client_path,
    process_waveforms=my_process_waveforms,
    channel="BH*"
)

# Main event (2020-02-06)
main_event_start_time = UTCDateTime("2020-02-06T18:33:14.88")
main_event_end_time = UTCDateTime("2020-02-6T18:34:40.00")

mseeds_main_event_path = os.path.join("mseeds", "main_event")

download_from_stations(
    CLIENT,
    network="PF",
    stations=["FRE", "HIM", "TTR", "PER"],
    events_times=[(main_event_start_time, main_event_end_time)],
    save_path=mseeds_main_event_path,
    process_waveforms=my_process_waveforms
)

download_from_stations(
    CLIENT,
    network="G",
    stations=["RER"],
    events_times=[(main_event_start_time, main_event_end_time)],
    save_path=mseeds_main_event_path,
    process_waveforms=my_process_waveforms
)
