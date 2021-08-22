""" Application of STA-LTA and Kurtosis on data from 01/02/2020 to 11/02/2020 :
    - Determination of global start of events triggered by simultaned stations
    - Determination of specific starts of trainwaves """

import os
import pickle

from tqdm import tqdm

from obspy import UTCDateTime

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms
from locevdet.trainwaves_arrivals.band_freq_interest import freq_band_interest
from locevdet.utils import clean_utc_str
# from locevdet.visualisation.plots import plot_stalta_per_event, plot_kurtosis_per_event

start_period = UTCDateTime("2020-02-01T00:00:00")
end_period = UTCDateTime("2020-02-11T00:00:00")
#Save to pickle format
period_name = '_'.join((clean_utc_str(start_period), clean_utc_str(end_period)))
dictname = f"eventlist_{period_name}.pkl"
eventlist_pickle_filepath = \
    os.path.join('dictionaries_data', '01-02-2020_11-02-2020', 'rockfall', 'eventlist', dictname)
# Paths of seismograms
mseeds_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "rockfall", "average_strong")


# Seismograms concerned
## (see locevdet.examples.casse_riv_est.all_seismo python file to modify conditions)
all_seismogram_reduced = list_seismograms(mseeds_path)

stalta_detection_events = True
add_RER = False
add_matlab_results = False
kurtosis_tool_params = False
kurtosis = True
envelope_enddetection = False
descriptors = False
azimut_distribution = False

analysis_plots = True

check_pickle_content = False
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
        "minimum_time" : 2,

        "rolling_max_window" : 50  # in seconds
    }
if stalta_detection_events:
    event_list, duplicate_ev_removed, false_ev_removed, _ = stalta_detect_events(
        folder_in=mseeds_path,
        all_seismogram=all_seismogram_reduced,
        low_snr_remove=False,
        **config_stalta
    )

    print(f"{len(event_list)} Events")
    print(f"{len(duplicate_ev_removed)} dictionnaires ont été supprimés, dont : \
        {duplicate_ev_removed}")
    print(f"{len(false_ev_removed)} dictionnaires ont été supprimés, dont : {false_ev_removed}")

    # Plots
    ## Per event
    save_stalta_fig_path = os.path.join("captures", "01-02-2020_11-02-2020", "rockfall", "stalta")
    print("STA-LTA captures in process")
    for event in tqdm(event_list):
        event.plot(type_graph='stalta', save_fig_path=save_stalta_fig_path, show=False)
    
    #Save to pickle format
    event_list.save(eventlist_pickle_filepath, override=True)


############################## ADD RER for each event ###########################################
if add_RER:
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

            for event in eventlist:
                event.add_station_trainwaves(mseeds_folder=mseeds_path, network='G', station='RER')

                eventlist.save(eventlist_pickle_filepath, override=True)

############################## ADD matlab's results to Evenlist ###################################
if add_matlab_results :
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

            folder_matlab_path = os.path.join('results_matlab', 'data_matlab')
            eventlist.set_matlab_data(folder_matlab_path, max_time_difference=2)

            # Filter on each trimmed traces
            save_hist_path = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall')
            freqmin = freq_band_interest(eventlist, save_fig_path=save_hist_path, show=False)[0]
            freqmax = freq_band_interest(eventlist, save_fig_path=save_hist_path, show=False)[1]

            save_clustering = os.path.join("captures", "01-02-2020_11-02-2020", \
                "rockfall", "clustering")
            for event in eventlist:
                event.fc_td_clustering(save_clustering)
            
            eventlist.save(eventlist_pickle_filepath, override=True)

##################################### KURTOSIS #####################################################
# The tool to determine parameters for kurtosis
config_kurtosis_tool = {
        "thr_on_init" : 0.75,
        "thr_off_init" : 0.25,
        "win_init" : 100                  # in seconds
    }
if kurtosis_tool_params :
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

            for event in eventlist:
                event.kurt_param_sliders(**config_kurtosis_tool)

config_kurtosis = {
        "window" : 200,                   # in seconds
        "threshold_on" : 0.35
    }
if kurtosis :
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

            for event in eventlist:
                for _, trainwave in event.trainwaves.items():
                    try :
                        trainwave.freqmin_interest = freqmin
                        trainwave.freqmax_interest = freqmax
                    except NameError:
                        freqmin, freqmax = 2, 10
                        trainwave.freqmin_interest = freqmin
                        trainwave.freqmax_interest = freqmax
                    
                    kurtosis_data = trainwave.kurtosis(**config_kurtosis)

            # Plots
            ## Per event
            # save_kurt_fig_path = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall', 'kurtosis')
            save_kurt_fig_path = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall', \
                'kurtosis_rapport')
            for event in eventlist:
                event.plot(type_graph='kurtosis', save_fig_path=save_kurt_fig_path, \
                    show=False, **config_kurtosis)
        
        #Save to pickle format
        eventlist.save(eventlist_pickle_filepath, override=True)

# ################################### ENVELOPE AND End event ######################################
config_envelope={
    "rolling_max_window" : config_stalta["rolling_max_window"],  # number of points 
    "thr_snr_purcent" : 1.1,           # 1 = 100%
    "time_restricted" : 5,             # in seconds
    "time_inspect_startglobal" : 1               # in seconds
}
if envelope_enddetection:
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)
          
            save_envelope_fig_path = os.path.join('captures', '01-02-2020_11-02-2020', \
                'rockfall', 'envelope_rapport')
            for event in eventlist:
                for _, trainwave in event.trainwaves.items():
                    if trainwave.station.name != 'RER':
                        trainwave.endtime_detection(\
                            rolling_max_window=config_envelope['rolling_max_window'],
                            time_restricted=config_envelope['time_restricted'],
                            time_inspect_startglobal=config_envelope['time_inspect_startglobal'],
                            thr_snr_purcent=config_envelope['thr_snr_purcent'])
            print("Envelope captures")
            for event in tqdm(eventlist):
                event.plot(type_graph='envelope', save_fig_path=save_envelope_fig_path, \
                    show=False, **config_envelope)
        
            #Save to pickle format
            eventlist.save(eventlist_pickle_filepath, override=True)

################################ DESCRIPTORS ######################################################
if descriptors:
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

            for event in eventlist:
                for _, trainwave in event.trainwaves.items():
                    if trainwave.station.name != 'RER':
                        print('snr and noise level:')
                        print(f"snr {trainwave.snr} and noise {trainwave.noise_level}")

            # Ratio form and duration descriptors
            for event in eventlist:
                for _, trainwave in event.trainwaves.items():
                    if trainwave.station.name != 'RER':
                        trainwave.endtime_detection(\
                            rolling_max_window=config_envelope['rolling_max_window'],
                            time_restricted=config_envelope['time_restricted'],
                            time_inspect_startglobal=config_envelope['time_inspect_startglobal'],
                            thr_snr_purcent=config_envelope['thr_snr_purcent'])

                        trainwave.form_ratio_and_duration(\
                            rolling_max_window=config_envelope['rolling_max_window'])

            # Number of picks (descriptor)
            for event in eventlist:
                for _, trainwave in event.trainwaves.items():
                    if trainwave.station.name != 'RER':
                        trainwave.number_picks(\
                            rolling_max_window=config_envelope['rolling_max_window'])
        
            #Save to pickle format
            eventlist.save(eventlist_pickle_filepath, override=True)

###################################  AZIMUTS per station ##########################################
if azimut_distribution:
    stations_with_azimut = ['HIM', 'FRE', 'RER']
    save_fig_azimut = os.path.join('captures', '01-02-2020_11-02-2020', 'rockfall', 'azimut_dist')

    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

            eventlist.azimut_hist(
                list_stations=stations_with_azimut, 
                save_fig_path=save_fig_azimut
                )

            eventlist.save(eventlist_pickle_filepath, override=True)

######################### COMPARAISON specific starts of trainwaves ############################
# Plots 
save_compare_fig_path = \
        os.path.join("captures", "01-02-2020_11-02-2020", "rockfall", "comparaison_ti")
if analysis_plots :
    all_plots_type = [\
        # "ti_with_snr", 
        "snr_in_fct_diff_ti", 
        # "hist_diff_ti", 
        # "dt", 
        # "temporal_snr_diff_ti"
    ]

    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)
            for type_plot in all_plots_type :
                eventlist.plots_time_compare_HIM_FRE(
                    type_graph=type_plot, 
                    save_fig_path=save_compare_fig_path, 
                    show=False
                )

################################### TO CHECK EVENTLIST PICKLE ####################################
if check_pickle_content :
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'rockfall', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)
            for event in eventlist:
                for _, trainwave in event.trainwaves.items():
                    print("trainwave.snr :", trainwave.snr)
                    print("trainwave.noise_level :", trainwave.noise_level)
                    print("trainwave.duration :", trainwave.duration)