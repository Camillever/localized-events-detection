""" Kurtosis application : 
with a tool to finf wisely paramaters to have starts of trainwaves"""

import os

from obspy import UTCDateTime

from locevdet.utils import freq_band_interest
from locevdet.trainwaves_arrivals.kurtosis_test import kurtosistest_for_all_seismograms
from locevdet.trainwaves_arrivals.kurtosis_test import kurtosis_starts_extraction

mseeds_path = os.path.join("mseeds", "RESIF") # Seismograms
sta_lta_save_path = os.path.join("sta_lta", "JSON") # Trainwaves dictionaries

PRE_TRIGGER = 5 # in seconds
POST_TRIGGER = 5 # in seconds

#Parameters
win = 200 # in seconds
THR_ON = 0.35

folder_matlab = os.path.join('results_matlab', 'data_matlab')
FREQMIN, FREQMAX = freq_band_interest(folder_matlab)

# To save figures and dictionnary
save_path = os.path.join("kurtosis")


# # TOOL to find parameters 
# kurtosistest_for_all_seismograms(
#     mseeds_path, 
#     sta_lta_save_path, 
#     PRE_TRIGGER, 
#     POST_TRIGGER,  
#     win=win
# )

# Extraction of starts when parameters were determined
## MSEEDS

all_seismogram = [
        filename for filename in os.listdir(mseeds_path)
        if not filename.startswith('G_RER') 
        if not filename.startswith('PF_TTR')  
    ]
    
    # Remove specific station on some period
    filenames = []
    all_seismogram_reduced = []
    for filename in all_seismogram:
        starttime = UTCDateTime(get_info_from_mseedname(filename, 'starttime'))
        station = get_info_from_mseedname(filename, 'station')
        if filename.startswith('PF_NSR') and starttime < UTCDateTime("2020-02-07T10:32:10"):
            filenames.append(filename)
        else:
            all_seismogram_reduced.append(filename)


kurtosis_starts_extraction(
    mseeds_path, 
    sta_lta_save_path, 
    PRE_TRIGGER, 
    POST_TRIGGER, 
    thr_on=THR_ON, 
    freqmin=FREQMIN, 
    freqmax=FREQMAX, 
    win=win, 
    save_path=save_path
    )
