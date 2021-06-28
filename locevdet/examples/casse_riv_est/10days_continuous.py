""" Script to identify events in continuous seismograms """

import os
from tqdm import tqdm

from obspy import UTCDateTime

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.utils import get_info_from_mseedname, ref_duration_to_time_window, clean_utc_str
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

    # To remove low snr events
    "rolling_max_window" : 50,  # in seconds
    "general_thr" : 3,
    "pre_trigger" : 10,        # in seconds
    "post_trigger" : 50        # in seconds
}

## Divide period to avoid MemoryError
periods_dates = [ref_duration_to_time_window(
    time_reference=start_time_continuous, 
    duration=86400,    # per day
    time_offset=(0, 0))] 
check_date = periods_dates[0][1]
while check_date < end_time_continuous :
    periods_dates.append(ref_duration_to_time_window(
    time_reference=periods_dates[-1][1], 
    duration=86400, 
    time_offset=(0,0)))
    check_date = periods_dates[-1][1]

for (start_period, end_period) in periods_dates:
    seismograms_in_period = [
        mseed for mseed in continuous_seismograms
        if get_info_from_mseedname(mseed)['starttime'] >= start_period and \
            get_info_from_mseedname(mseed)['endtime'] <= end_period
    ]
    events_list_continuous, duplicate_ev, false_ev, lowsnr_removed = stalta_detect_events(
        folder_in=mseeds_path_continuous,
        all_seismogram=seismograms_in_period,
        low_snr_remove=True,
        **config_stalta
    )

    print(f"Events ({len(events_list_continuous)})")
    print(f"{len(duplicate_ev)} événements ont été supprimés")
    print(f"{len(false_ev)} événements ont été supprimés")
    print(f"{len(lowsnr_removed)} événements ont été supprimés")

    # # Remove all low snr events
    # low_snr_events = events_list_continuous.delete_low_snr(
    #     general_thr=3, rolling_max_window=4)
    # print(f"{len(low_snr_events)} événements ont été supprimés, dont : {low_snr_events}")

    # # Add band of frequence of interest and Calculation of the snr and addition to trainwaves
    # for event in events_list_continuous:
    #     for _, trainwave in event.trainwaves.items():
    #         trainwave.freqmin_interest = config_stalta["freqmin"]
    #         trainwave.freqmax_interest = config_stalta["freqmax"]
    #         if trainwave.station.name != 'RER':
    #             trainwave.snr_calculation(rolling_max_window=config_stalta['rolling_max_window'])

    # Plots
    ## Per event
    save_stalta_fig_path = os.path.join("captures","01-02-2020_11-02-2020","continuous","stalta")
    for event in tqdm(events_list_continuous):
        event.plot(
            type_graph='stalta',
            seismograms_type='continuous',
            save_fig_path=save_stalta_fig_path,
            freqmin=config_stalta["freqmin"],
            freqmax=config_stalta["freqmax"],
            thr_on=config_stalta["thr_on"],
            pre_trigger=config_stalta["pre_trigger"], # in seconds (trace_trimmed)
            post_trigger=config_stalta["post_trigger"]  # in seconds (trace_trimmed)
        )


################################## SAVE EVENLIST ################################################
    #Save to pickle format
    period_name = '_'.join((clean_utc_str(start_period), clean_utc_str(end_period)))
    dictname = f"eventlist_continuous_{period_name}.pkl"
    eventlist_pickle_filepath = \
        os.path.join('dictionaries_data', '01-02-2020_11-02-2020', 'continuous', 'eventlist', dictname)
    events_list_continuous.save(eventlist_pickle_filepath, override=True)

    # # To visualize this save

    # with open(eventlist_pickle_filepath, "rb") as content:
    #     evenlist = pickle.load(content)
    # print(len(evenlist))
   
################################### DESCRIPTORS ###################################################
