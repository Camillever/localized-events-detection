""" Application of STA-LTA and Kurtosis on data from 01/02/2020 to 11/02/2020 :
    - Determination of global start of events triggered by simultaned stations
    - Determination of specific starts of trainwaves """

import os
import pickle

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.trainwaves_arrivals.kurtosis_test import kurtosistest_for_all_seismograms
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms
from locevdet.trainwaves_arrivals.band_freq_interest import freq_band_interest
# from locevdet.visualisation.plots import plot_stalta_per_event, plot_kurtosis_per_event

# Paths of seismograms
mseeds_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "rockfall", "average_strong")

# Seismograms concerned
## (see locevdet.examples.casse_riv_est.all_seismo python file to modify conditions)
all_seismogram_reduced = list_seismograms(mseeds_path)


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
save_stalta_fig_path = os.path.join("captures", "01-02-2020_11-02-2020", "rockfall", "stalta")
# for event in event_list:
#     event.plot(type_graph='stalta', save_fig_path=save_stalta_fig_path, show=False)

############################## ADD matlab's results to Evenlist ###################################

folder_matlab_path = os.path.join('results_matlab', 'data_matlab')
event_list.set_matlab_data(folder_matlab_path, max_time_difference=2)

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
save_hist_path = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall')
###

config_kurtosis = {
    "window" : 200,                   # in seconds
    "threshold_on" : 0.35
}

# Filter on each trimmed traces
freqmin = freq_band_interest(event_list, save_fig_path=save_hist_path, show=False)[0]
freqmax = freq_band_interest(event_list, save_fig_path=save_hist_path, show=False)[1]

for event in event_list:
    for _, trainwave in event.trainwaves.items():
        trainwave.freqmin_interest = freqmin
        trainwave.freqmax_interest = freqmax
        kurtosis_data = trainwave.kurtosis(**config_kurtosis)

# Plots
## Per event
save_kurt_fig_path = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall', 'kurtosis')
for event in event_list:
    event.plot(type_graph='kurtosis', save_fig_path=save_kurt_fig_path, show=False, **config_kurtosis)

# ################################### ENVELOPE AND SNR #############################################
config_envelope={
    "rolling_max_window" : 1,          # in seconds
    "thr_snr_purcent" : 1.1,           # 1 = 100%
    "time_intervall_inspect" : [1,10], # in seconds
    "window_inspect" : 5               # in seconds
}
# # Calculation of the snr
# for event in event_list:
#     for _, trainwave in event.trainwaves.items():
#         trainwave.snr_calculation(rolling_max_window=config_envelope['rolling_max_window'],
#             time_intervall_inspect=config_envelope['time_intervall_inspect'],
#             window_inspect=config_envelope['window_inspect'])
        
# # Determination of ends of trainwaves
# for event in event_list:
#     for _, trainwave in event.trainwaves.items():
#         trainwave.endtime_detection(rolling_max_window=config_envelope['rolling_max_window'],
#             time_intervall_inspect=config_envelope['time_intervall_inspect'],
#             thr_snr_purcent=config_envelope['thr_snr_purcent'])

# Plots
save_envelope_fig_path = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall', 'envelope')
# for event in event_list:
#     event.plot(type_graph='envelope', save_fig_path=save_envelope_fig_path, show=False, **config_envelope)

# ################################## SAVE EVENLIST ##################################################
# Save to pickle format
# dictname = 'eventlist_01-02-2020_11-02-2020.pkl'
# eventlist_pickle_filepath = \
#     os.path.join('dictionaries_data', '01-02-2020_11-02-2020', 'rockfall', 'eventlist', dictname)
# event_list.save(eventlist_pickle_filepath, override=True)

# # To visualize this save

# with open(eventlist_pickle_filepath, "rb") as content:
#     evenlist = pickle.load(content)
# print(evenlist)
# print(len(evenlist))


######################### COMPARAISON specific starts of trainwaves ############################
# Plots 
save_compare_fig_path = \
    os.path.join("captures", "01-02-2020_11-02-2020", "rockfall", "comparaison_ti")

# "ti_with_snr"
event_list.plots_time_compare_HIM_FRE(
    type_graph="ti_with_snr", 
    save_fig_path=save_compare_fig_path, 
    show=False
)

# # " diff_ti_in_fct_snr"


# "hist_diff_ti"
event_list.plots_time_compare_HIM_FRE(
    type_graph="hist_diff_ti", 
    save_fig_path=save_compare_fig_path, 
    show=False
)

# "dt"
event_list.plots_time_compare_HIM_FRE(
    type_graph="dt", 
    save_fig_path=save_compare_fig_path, 
    show=False
)