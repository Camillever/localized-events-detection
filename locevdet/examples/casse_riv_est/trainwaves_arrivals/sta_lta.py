""" Application of STA-LTA and Kurtosis on data from 01/02/2020 to 11/02/2020 :
    - Determination of global start of events triggered by simultaned stations
    - Determination of specific starts of trainwaves """

import os

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.trainwaves_arrivals.kurtosis_test import kurtosistest_for_all_seismograms
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms
from locevdet.utils import freq_band_interest
from locevdet.visualisation.plots import plot_stalta_per_event, plot_kurtosis_per_event

# Paths
mseeds_path = os.path.join("mseed_RESIF", "seismograms_all_stations_average_and_strong")
sta_lta_save_path = os.path.join("sta_lta")



# Seismograms concerned
## (see locevdet.examples.casse_riv_est.all_seismo python file to modify conditions)
all_seismogram_reduced = list_seismograms(mseeds_path)


# Plots all seismograms
## Per event


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

event_list, duplicate_ev_removed, false_ev_removed = stalta_detect_events(
    folder_in=mseeds_path,
    all_seismogram=all_seismogram_reduced,
    **config_stalta
)
print(f"Events ({len(event_list)}): {event_list}")
print(f"{len(duplicate_ev_removed)} dictionnaires ont été supprimés, dont : {duplicate_ev_removed}")
print(f"{len(false_ev_removed)} dictionnaires ont été supprimés, dont : {false_ev_removed}")

# Plots
## Per event
save_stalta_fig_path = os.path.join("sta_lta", "PNG_01-02-2020_11-02-2020")
# plot_stalta_per_event(event_list, save_fig_path=save_stalta_fig_path, show=False)

############################## ADD matlab's results to Evenlist ###################################

folder_matlab_path = os.path.join('results_matlab', 'data_matlab')
event_list.set_matlab_data(folder_matlab_path)


# for event in event_list:
#     for _, trainwave in event.trainwaves.items():
#         if trainwave.matlab_data is not None : 
#             print(trainwave.matlab_data['trainwave']['initial_time'])
#             print(trainwave.matlab_data['trainwave']['fmin'])



##################################### KURTOSIS #####################################################
# The tool to determine parameters for kurtosis

# kurtosistest_for_all_seismograms(
#     mseeds_path, 
#     sta_lta_save_path, 
#     PRE_TRIGGER, 
#     POST_TRIGGER,  
#     win=win
# )


# Extraction of specific starts of trainwaves
save_hist_path = os.path.join('captures', 'histograms')
###

config_kurtosis = {

    # Filter
    "freqmin": freq_band_interest(event_list, save_fig_path=save_hist_path, show=False)[0],
    "freqmax": freq_band_interest(event_list, save_fig_path=save_hist_path, show=False)[1],

    # kurtosis parameters
    "window" : 200, 
    "threshold_on" : 0.35
}
# for event in event_list:
#     for _, trainwave in event.trainwaves.items():
#         kurtosis_data = trainwave.kurtosis(**config_kurtosis)


# Plots
## Per event
save_kurt_fig_path = os.path.join('kurtosis', 'PNG_01-02-2020_11-02-2020')
# plot_kurtosis_per_event(event_list, save_fig_path=save_kurt_fig_path, show=False, **config_kurtosis)

# ################################### ENVELOPE AND SNR #############################################
# # Add Hilbert envelop of trainwave's windows to EventList
save_env_fig_path = os.path.join('envelope', 'PNG_01-02-2020_11-02-2020')
for event in event_list:
    event.plot_envelope(save_fig_path=save_env_fig_path, show=False)

# Add SNR of trainwave's windows to EventList

# for event in event_list:
#     for trainwave in event.trainwave:
#         trainwave.nsr(rolling_max_window=) # To define in Event class


# Detection of end of trainwaves + Add to EventList

# Plots
## Per event


# ################################## SAVE EVENLIST ##################################################
# # Save to pickle format
# eventlist_pickle_filepath = os.path.join('data', 'eventlist', 'eventlist_01-02-2020_11-02-2020.pkl')
# event_list.save(eventlist_pickle_filepath, override=True)

## To visualize this save
# import pickle
# with open(eventlist_pickle_filepath, "rb") as content:
#     evenlist = pickle.load(content)
# # print(evenlist)
# print(len(evenlist))


######################### COMPARAISON specific starts of trainwaves ############################
# Plots 
# save_compare_fig_path = os.path.join("captures", "comparaison_ti")
# plots_comparaison_starts_HIM_FRE(eventlist, save_compare_fig_path, show=True)
