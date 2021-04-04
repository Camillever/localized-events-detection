""" To apply STA-LTA"""

import os

from locevdet.trainwaves_arrivals.sta_lta import stalta_per_event_coincidence_trigger
from locevdet.trainwaves_arrivals.sta_lta import trainwaves_too_close_remove

#paths
mseeds_path = os.path.join("mseeds", "RESIF")
sta_lta_save_path = os.path.join("sta_lta")

config = {
    # Filter
    "freqmin": 2,
    "freqmax": 10,

    # Sta-lta parameters
    "nsta_time": 1, # in seconds
    "nlta_time": 30, # in seconds
    "thr_on": 2.9,
    "thr_off": 1,
    "thr_coincidence_sum": 2,
}

stalta_per_event_coincidence_trigger(
    folder_in=mseeds_path,
    save_path=sta_lta_save_path,
    **config
)

# To remove useless trainwaves, too close of each others
json_sta_lta_save_path = os.path.join(sta_lta_save_path, "JSON")
intervall = 2 # in seconds

trainwaves_too_close_remove(json_sta_lta_save_path, intervall)