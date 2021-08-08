""" Script to identify events in continuous seismograms """

import shutil, os

import pickle
from tqdm import tqdm

from obspy import UTCDateTime
from obspy.core import read

from scipy.io import savemat

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.utils import ref_duration_to_time_window, clean_utc_str
from locevdet.utils import get_info_from_mseedname, localisation
from locevdet.waveform_processing import trim_trace
##################################################################################################
################################# 01/02/2020 - 11/02/2020 ########################################

mseeds_path_continuous = os.path.join("seismograms","01-02-2020_11-02-2020","continuous")
start_time_continuous = UTCDateTime("2020-02-01T00:00:00.00")
end_time_continuous = UTCDateTime("2020-02-11T00:00:00.00")


date_07fev = UTCDateTime("2020-02-07T00:00:00.00")

before_07fev_seismograms =[
    mseed for mseed in os.listdir(mseeds_path_continuous)
    if get_info_from_mseedname(mseed)['starttime'] < date_07fev
    if not mseed.startswith('PF_NSR') or not mseed.startswith('G_RER')
    ]
after_07fev_seismograms =[
    mseed for mseed in os.listdir(mseeds_path_continuous)
    if get_info_from_mseedname(mseed)['starttime'] > date_07fev
    if not mseed.startswith('G_RER')
    ]
continuous_seismograms = [*before_07fev_seismograms, *after_07fev_seismograms]


download_eventlist = False
descriptors_eventlist = False
convertion_mat = True

check_pickle_content = False


################################## EVENTS DETECTION BY STA-LTA ####################################

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

    # To remove low snr events
    "rolling_max_window" : 50,  # number of points (100pts = 1s if sampling_rate = 100Hz)
    "general_thr" : 3.6,
    "pre_trigger" : 10,        # in seconds
    "post_trigger" : 50        # in seconds
}

## Divide period to avoid MemoryError
periods_dates = [ref_duration_to_time_window(
    time_reference=start_time_continuous,
    duration=86400/2,    # per day
    time_offset=(0, 0))] 
check_date = periods_dates[0][1]
while check_date < end_time_continuous :
    periods_dates.append(ref_duration_to_time_window(
    time_reference=periods_dates[-1][1], 
    duration=86400/2, 
    time_offset=(0,0)))
    check_date = periods_dates[-1][1]
print("periods_dates :", periods_dates)

if download_eventlist :
    for (start_period, end_period) in periods_dates:
        print(f"{start_period} - {end_period}")
        seismograms_in_period = [
            mseed for mseed in continuous_seismograms
            if get_info_from_mseedname(mseed)['starttime'] >= start_period and \
                get_info_from_mseedname(mseed)['endtime'] <= end_period
        ]
        events_list_continuous, duplicate_ev, false_ev, lowsnr_removed = stalta_detect_events(
            folder_in=mseeds_path_continuous,
            all_seismogram=seismograms_in_period,
            low_snr_remove=False,
            **config_stalta
        )

        print(f"{len(duplicate_ev)} événements ont été supprimés")
        print(f"{len(false_ev)} événements ont été supprimés")
        print(f"{len(lowsnr_removed)} événements ont été supprimés")
        print(f"Events ({len(events_list_continuous)})")

        # Add band of frequence of interest
        print('Add band of frequence')
        for event in tqdm(events_list_continuous):
            for _, trainwave in event.trainwaves.items():
                trainwave.freqmin_interest = config_stalta["freqmin"]
                trainwave.freqmax_interest = config_stalta["freqmax"]

        # Plots
        ## Per event
        print('Save plots of stalta detection')
        save_stalta_fig_path = os.path.join("captures", "01-02-2020_11-02-2020",\
            "continuous", "stalta")
        for event in tqdm(events_list_continuous):
            event.plot(
                type_graph='stalta',
                seismograms_type='continuous',
                save_fig_path=save_stalta_fig_path,
                freqmin=config_stalta["freqmin"],
                freqmax=config_stalta["freqmax"],
                thr_on=config_stalta["thr_on"],
                pre_trigger=config_stalta["pre_trigger"], # in seconds (trace_trimmed)
                post_trigger=config_stalta["post_trigger"],  # in seconds (trace_trimmed),
                rolling_max_window=config_stalta["rolling_max_window"]
            )

        # Remove all low snr events (if low_snr_remove is False in stalta_detect_events fct)
        print('Remove all low snr events')
        low_snr_events = events_list_continuous.delete_low_snr(
            general_thr=config_stalta["general_thr"], 
            rolling_max_window=config_stalta["rolling_max_window"])
        print(f"{len(low_snr_events)} événements ont été supprimés")

        print(f"Events ({len(events_list_continuous)})")

        # To remove events captures
        all_start_global = []
        for event in events_list_continuous:
            all_start_global.append(clean_utc_str(event.start_global))
        print('all_start_global len :', len(all_start_global))

        all_stalta_captures = [ 
            fig for fig in os.listdir(save_stalta_fig_path)
            if fig.endswith('.png')
        ]
        for capture in all_stalta_captures:
            title_png = capture.split('_')
            start_global_event = clean_utc_str(UTCDateTime(title_png[1].split('.')[0]))
            if start_global_event >= start_period and start_global_event <= end_period \
                and start_global_event not in all_start_global :
                # os.remove(os.path.join(save_stalta_fig_path, capture))
                # Déplacer un fichier du répertoire rep1 vers rep2
                shutil.move(\
                    os.path.join(save_stalta_fig_path, capture), 
                    os.path.join("captures", "01-02-2020_11-02-2020",\
            "continuous","stalta_removed", capture))

    ################################## SAVE EVENLIST ##############################################
        #Save to pickle format
        period_name = '_'.join((clean_utc_str(start_period), clean_utc_str(end_period)))
        dictname = f"eventlist_continuous_{period_name}.pkl"
        eventlist_pickle_filepath = \
            os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
                'continuous', 'eventlist', dictname)
        events_list_continuous.save(eventlist_pickle_filepath, override=True)

###################################################################################################
################################### DESCRIPTORS ###################################################
# Load
if descriptors_eventlist:
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist'))
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        start_period = eventlist_file.split('_')[2]
        end_period = (eventlist_file.split('_')[3]).split('.')[0]
        print(f"start and end period : {start_period} -{end_period}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)

        ################################### ENVELOPE AND SNR #######################################
        config_envelope={
            "rolling_max_window" : 4,          # in seconds
            "thr_snr_purcent" : 1.1,           # 1 = 100%
            "time_restricted" : 5,             # in seconds
            "time_inspect_startglobal" : 1               # in seconds
        }
        
                
        # Determination of ends of trainwaves
        print('Descriptors calculation')
        for event in tqdm(eventlist):
            for _, trainwave in event.trainwaves.items():
                if trainwave.station.name != 'RER':
                    trainwave.endtime_detection(rolling_max_window=config_envelope['rolling_max_window'],
                        time_restricted=config_envelope['time_restricted'],
                        time_inspect_startglobal=config_envelope['time_inspect_startglobal'],
                        thr_snr_purcent=config_envelope['thr_snr_purcent'])
                    
                    trainwave.form_ratio_and_duration(rolling_max_window=config_envelope['rolling_max_window'])

                    trainwave.number_picks(rolling_max_window=config_envelope['rolling_max_window'])
            # Plots
        save_envelope_fig_path = os.path.join('captures', '01-02-2020_11-02-2020', 'continuous', 'envelope')
        print('Save envelope plots')
        for event in tqdm(eventlist):
            event.plot(type_graph='envelope', save_fig_path=save_envelope_fig_path, \
                show=False, **config_envelope)

    ################################## SAVE EVENLIST ##############################################
        #Save to pickle format
        period_name = '_'.join((start_period, end_period))
        dictname = f"eventlist_continuous_{period_name}.pkl"
        eventlist_pickle_filepath = \
            os.path.join('dictionaries_data', '01-02-2020_11-02-2020', 'continuous', 'eventlist', dictname)
        eventlist.save(eventlist_pickle_filepath, override=True)



################################################################################################
######################################### CONVERTION TO .MAT ##################################
all_seismogram = [
        filename for filename in os.listdir(mseeds_path_continuous)
        if filename.startswith('PF_FRE') 
        or filename.startswith('PF_HIM') 
        or filename.startswith('G_RER')
    ]

stations_authorized = ['PF_HIM', 'PF_FRE', 'G_RER']


save_mat_continuous_events = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
    'continuous', 'seismograms_matlab')

if convertion_mat :
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
        'continuous', 'eventlist'))

    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)
            
            print('Convertion of (trimmed) events into mat file')
            for event in tqdm(eventlist):
                start_global = event.start_global
                pre_trigger=config_stalta["pre_trigger"] # in seconds (trace_trimmed)
                post_trigger=config_stalta["post_trigger"]  # in seconds (trace_trimmed)
                freqmin=config_stalta["freqmin"]
                freqmax=config_stalta["freqmax"]

                for station in stations_authorized :
                    all_seismograms_station = [
                        mseed for mseed in all_seismogram
                        if mseed.startswith(station)
                    ]

                    station_ckecked = None
                    for seismogram_file in all_seismograms_station:
                        network_seismo = get_info_from_mseedname(seismogram_file)['network']
                        station_seismo = get_info_from_mseedname(seismogram_file)['station']
                        starttime_seismo = get_info_from_mseedname(seismogram_file)['starttime']
                        endtime_seismo = get_info_from_mseedname(seismogram_file)['endtime']

                        if start_global > starttime_seismo and start_global < endtime_seismo \
                            and station_ckecked is None:
                            filepath = os.path.join(mseeds_path_continuous, seismogram_file)
                            seismogram = read(filepath)
                            filemat = {"E":{},"N":{},"Z":{}}
                            station_ckecked = 1

                            for _,seismo in enumerate(seismogram):
                                if seismo.stats.npts <= 50 :
                                    raise ('Shape Array incorrect')

                                seismogram_trim = trim_trace(\
                                    seismo, start_global, pre_trigger, post_trigger)
                                try :
                                    seismogram_trim_filt = seismogram_trim.copy().filter(\
                                        'bandpass',freqmin=freqmin, freqmax=freqmax)
                                except ValueError:
                                    seismogram_trim_filt = seismogram_trim

                                component = str(seismo.stats.component)
                                filemat[component] = {
                                        "network" : seismo.stats.network, 
                                        "station" : seismo.stats.station, 
                                        "channel" : seismo.stats.channel, 
                                        "starttime" : str(seismogram_trim_filt.stats.starttime), 
                                        "sampling_rate" : seismo.stats.sampling_rate, 
                                        "delta": seismo.stats.delta,
                                        "npts": seismogram_trim_filt.stats.npts, 
                                        "t": seismogram_trim_filt.times(), 
                                        "data": seismogram_trim_filt.data, 
                                        "localisation":    {
                                            "rgr92": {
                                                "x":localisation(seismogram_file)["x"], 
                                                "y":localisation(seismogram_file)["y"], 
                                                "elevation":localisation(seismogram_file)["z"]           
                                                        }, 
                                            "wgs94": {
                                                "latitude":localisation(seismogram_file)["latitude"], 
                                                "longitude":localisation(seismogram_file)["longitude"], 
                                                "elevation":localisation(seismogram_file)["z"]
                                                            }
                                                        }
                                        }
                            filename = '_'.join((network_seismo, station_seismo, clean_utc_str(start_global)))
                            print(filename)
                            savemat(os.path.join(save_mat_continuous_events,filename+".mat"), filemat)
                        continue


################################### TO CHECK EVENTLIST PICKLE ####################################
if check_pickle_content :
    all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist'))
    number_events = 0
    for eventlist_file in all_eventlists:
        print(f"Eventlist file in process :{eventlist_file}")
        eventlist_filepath = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist', eventlist_file)
        with open(eventlist_filepath, "rb") as content:
            eventlist = pickle.load(content)
            number_events += len(eventlist)

            print(f"{len(eventlist)} events for {eventlist_file}")
            print("")
    print(f"Total : {number_events} events" )
            # for event in eventlist:
            #     for _, trainwave in event.trainwaves.items():
            #         print(trainwave.start_global)