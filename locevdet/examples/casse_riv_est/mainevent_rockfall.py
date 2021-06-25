""" Script to pick and locate arrivals of the main event """
import os

import numpy as np

from obspy import UTCDateTime
from scipy.io import savemat

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.utils import get_info_from_mseedname

# Paths of seismograms
mseeds_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "rockfall", "average_strong")
start_periodtime = UTCDateTime("2020-02-06T17:00:00")
end_periodtime = UTCDateTime("2020-02-06T19:00:00")

main_event_seismograms = [ 
    mseed for mseed in os.listdir(mseeds_path)
    if get_info_from_mseedname(mseed)['starttime'] > start_periodtime and \
        get_info_from_mseedname(mseed)['endtime'] < end_periodtime
    if mseed.startswith('PF_HIM') or mseed.startswith('PF_FRE') or mseed.startswith('PF_PER')
    ]
print(main_event_seismograms)

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

mainevents_list, duplicate_ev_removed, false_ev_removed = stalta_detect_events(
    folder_in=mseeds_path,
    all_seismogram=main_event_seismograms,
    **config_stalta
)
print(f"{len(mainevents_list)} events have been detected by STA-LTA  on the main event : {mainevents_list}")
print(f"{len(duplicate_ev_removed)} dictionaries have been deleted, among : {duplicate_ev_removed}")
print(f"{len(false_ev_removed)} dictionaries have been deleted, among : {false_ev_removed}")

############################## ADD RER for each event ###########################################

for event in mainevents_list:
    event.add_station_trainwaves(mseeds_folder=mseeds_path, network='G', station='RER')

# Save freqmin and freqmax :
for event in mainevents_list:
    for _, trainwave in event.trainwaves.items():
        trainwave.freqmin_interest = config_stalta["freqmin"]
        trainwave.freqmax_interest = config_stalta["freqmax"]

############################## ADD matlab's results to Evenlist ###################################

folder_matlab_path = os.path.join('results_matlab', 'data_matlab')
mainevents_list.set_matlab_data(folder_matlab_path, max_time_difference=2)

####################################### Arrivals detected by kurtosis ############################
config_kurtosis = {
    "window" : 200,                   # in seconds
    "threshold_on" : 0.35
}

for event in mainevents_list:
    for _, trainwave in event.trainwaves.items():
        kurtosis_data = trainwave.kurtosis(**config_kurtosis)

################################### PICK MANUALLY ARRIVALS ########################################
# for event in mainevents_list:
#     event.pick_arrivals_manually()

###########################  LOCALISATION : Convertion of arrivals .mat ############################
order_stations = ['HIM', 'FRE', 'NSR', 'PER', 'TTR', 'RER'] 

starts_to_save = [
    # 'start_manual', 
    'ti', 
    'td', 
    'start_specific']
save_folderpath = os.path.join("dictionaries_data", "main_event", "matlab")

from locevdet.eventlist import EventList
def starts_save_to_mat(eventlist:EventList, type_start, folder_save:str):
    number_events = len(mainevents_list)
    mainevents = np.zeros((number_events, len(order_stations)), dtype=object)
    for n_ev, event in enumerate(eventlist):
        event_row = ['None', 'None', 'None', 'None', 'None', 'None']
        for number, station in enumerate(order_stations) :
            for i, trainwave in event.trainwaves.items():
                
                if trainwave.station.name == station:
                    start = None
                    if type_start == 'start_manual' :
                        if trainwave.start_specific_manual is not None:
                            start = trainwave.start_specific_manual
                    elif type_start == 'ti':
                        if trainwave.matlab_data is not None :
                            if trainwave.matlab_data['trainwave'] is not None:
                                start = trainwave.matlab_data['trainwave']['initial_time']
                    elif type_start == 'td':
                        if trainwave.matlab_data is not None :
                            if trainwave.matlab_data['trainwave'] is not None:
                                start = trainwave.matlab_data['trainwave']['central_time']
                    elif type_start == 'start_specific':
                        if trainwave.kurtosis_data['start_specific'] is not None:
                            start = trainwave.kurtosis_data['start_specific']
                    event_row[number] = str(start)
                elif trainwave.station.name != station:
                    continue          
        mainevents[n_ev,:] = event_row
    print(mainevents)
    save_filename = f"mainevent_{type_start}.mat"
    savemat(os.path.join(folder_save, save_filename), mdict={'mainevents': mainevents})

for start in starts_to_save:
    starts_save_to_mat(mainevents_list, start, save_folderpath)
