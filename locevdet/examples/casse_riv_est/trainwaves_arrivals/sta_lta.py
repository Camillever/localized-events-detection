""" Application of STA-LTA : 
    Determination of global start of trainwaves triggered by simultaned stations """

import os


from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.trainwaves_arrivals.sta_lta import remove_too_close_trainwaves
from locevdet.trainwaves_arrivals.sta_lta import remove_border_stalta_false_trainwaves
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms

#Paths
mseeds_path = os.path.join("mseed_RESIF", "seismograms_all_stations_average_and_strong")
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


all_seismogram_reduced = list_seismograms(mseeds_path)

event_list = stalta_detect_events(
    folder_in=mseeds_path,
    all_seismogram=all_seismogram_reduced,
    **config
)
print(f"Events ({len(event_list)}): {event_list}")

# To remove useless trainwaves, too close of each others
dict_removed = remove_too_close_trainwaves(event_list, minimum_time=2)
print(f"{len(dict_removed)} dictionnaires ont été supprimés, dont : {dict_removed}")

# To remove false trainwaves
events_removed = remove_border_stalta_false_trainwaves(event_list, config["nlta_time"])
print(f"{len(events_removed)} dictionnaires ont été supprimés, dont : {events_removed}")

# Save EventList
eventlist_pickle_filepath = os.path.join('data', 'eventlist', 'eventlist_01-02-2020_11-02-2020.pkl')
event_list.save(eventlist_pickle_filepath, override=True)

### To visualize this save
# import pickle
# with open(eventlist_pickle_filepath, "rb") as content:
#     evenlist = pickle.load(content)
# # print(evenlist)
# print(len(evenlist))
