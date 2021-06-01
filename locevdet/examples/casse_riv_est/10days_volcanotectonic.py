""" Application of STA-LTA from 01/02/2020 to 11/02/2020 :
    - Determination of global start of events triggered by simultaned stations
    - Convertion of pertinent volcanotectonic events """

import os

import pandas as pd

from obspy import UTCDateTime
from obspy.core import read

from scipy.io import savemat

from locevdet.trainwaves_arrivals.sta_lta import stalta_detect_events
from locevdet.utils import get_info_from_mseedname, localisation

# Paths of seismograms
mseeds_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "volcanotectonic")

all_seismogram = [
        filename for filename in os.listdir(mseeds_path)
        if filename.startswith('PF_FRE') 
        or filename.startswith('PF_HIM') 
        or filename.startswith('G_RER')
    ]
print("all_seismogram :", all_seismogram)

####################################### STA-LTA ##################################################

config_stalta = {
    # Filter
    "freqmin": 1,
    "freqmax": 50,

    # Sta-lta parameters
    "nsta_time": 1, # in seconds
    "nlta_time": 30, # in seconds
    "thr_on": 2.9,
    "thr_off": 1,
    "thr_coincidence_sum": 2,

    # Remove useless events (duplicate, false)
    "minimum_time" : 2
}

## Are nsta and nlta parameters pertinents ??

event_list, duplicate_ev_removed, false_ev_removed = stalta_detect_events(
    folder_in=mseeds_path,
    all_seismogram=all_seismogram,
    **config_stalta
)
print(f"Events ({len(event_list)}): {event_list}")
print(f"{len(duplicate_ev_removed)} dictionnaires ont été supprimés, dont : {duplicate_ev_removed}")
print(f"{len(false_ev_removed)} dictionnaires ont été supprimés, dont : {false_ev_removed}")

# Plots
## Per event
save_stalta_path = os.path.join("captures", "01-02-2020_11-02-2020", "volcanotectonic", "stalta")
# for event in event_list:
#     event.plot(type_graph='stalta', save_fig_path=save_stalta_path, show=False)


##################################### Convert into .mat ##########################################
time_events_detected = [
    (event.split('.')[0]).split('_')[1]
    for event in os.listdir(save_stalta_path)
    ]
seismograms_detected = []
for event in time_events_detected:
    seismogram_associated = []
    for seismogram in all_seismogram:
        starttime_seismo = get_info_from_mseedname(seismogram)['starttime']
        endtime_seismo = get_info_from_mseedname(seismogram)['endtime']
        station_seismo = get_info_from_mseedname(seismogram)['station']
        if UTCDateTime(event) > starttime_seismo and UTCDateTime(event) < endtime_seismo \
            and station_seismo not in seismogram_associated:
            seismograms_detected.append(seismogram)
            seismogram_associated.append(station_seismo)

# To be sure to not have same event on several seismograms
seismograms_detected = list(set(seismograms_detected)) 

save_mat_volcano_events = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
    'volcanotectonic', 'seismograms_matlab')

for seismogram_detected in seismograms_detected:
    filepath = os.path.join(mseeds_path, seismogram_detected)
    seismogram = read(filepath)
    # filemat ={}
    filemat = {"E":{},"N":{},"Z":{}}
    for _,seismo in enumerate(seismogram):
        component = str(seismo.stats.component)
        filemat[component] = {
                "network" : seismo.stats.network, 
                "station" : seismo.stats.station, 
                "channel" : seismo.stats.channel, 
                "starttime" : str(seismo.stats.starttime), 
                "sampling_rate" : seismo.stats.sampling_rate, 
                "delta": seismo.stats.delta,
                "npts": seismo.stats.npts, 
                "t": seismo.times(), 
                "data": seismo.data, 
                "localisation":    {
                    "rgr92": {
                        "x":localisation(seismogram_detected)["x"], 
                        "y":localisation(seismogram_detected)["y"], 
                        "elevation":localisation(seismogram_detected)["z"]           
                                }, 
                    "wgs94": {
                        "latitude":localisation(seismogram_detected)["latitude"], 
                        "longitude":localisation(seismogram_detected)["longitude"], 
                        "elevation":localisation(seismogram_detected)["z"]
                                    }
                                }
                }
    savemat(os.path.join(save_mat_volcano_events,seismogram_detected+".mat"), filemat)
