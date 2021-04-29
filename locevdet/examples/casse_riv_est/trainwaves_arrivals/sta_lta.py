""" Application of STA-LTA : 
    Determination of global start of trainwaves triggered by simultaned stations """

import os

from locevdet.trainwaves_arrivals.sta_lta import stalta_per_event_coincidence_trigger
from locevdet.trainwaves_arrivals.sta_lta import trainwaves_too_close_remove
from locevdet.trainwaves_arrivals.sta_lta import false_trainwaves

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
    "override":True
}

number_trainwaves = stalta_per_event_coincidence_trigger(
    folder_in=mseeds_path,
    save_path=sta_lta_save_path,
    **config
)
print(f"{number_trainwaves} dictionnaires ont été importés")

# To remove useless trainwaves, too close of each others
json_sta_lta_save_path = os.path.join(sta_lta_save_path, "JSON")
INTERVALL = 2 # in seconds

dict_removed = trainwaves_too_close_remove(json_sta_lta_save_path, INTERVALL)
print(f"{len(dict_removed)} dictionnaires ont été supprimés, dont : {dict_removed}")

# To remove false trainwaves
OFFSET = 2 # in seconds
dict_removed = false_trainwaves(mseeds_path,
    json_sta_lta_save_path, config["nlta_time"], time_offset=OFFSET)
print(f"{len(dict_removed)} dictionnaires ont été supprimés, dont : {dict_removed}")