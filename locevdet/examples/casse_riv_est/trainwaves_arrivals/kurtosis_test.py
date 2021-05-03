""" Kurtosis application : 
with a tool to finf wisely paramaters to have starts of trainwaves"""

import os

from obspy import UTCDateTime

from locevdet.utils import freq_band_interest
from locevdet.trainwaves_arrivals.kurtosis_test import kurtosistest_for_all_seismograms
from locevdet.trainwaves_arrivals.kurtosis_test import kurtosis_starts_extraction
from locevdet.trainwaves_arrivals.kurtosis_test import specific_starts_extraction_by_kurtosis
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms

mseeds_path = os.path.join("mseed_RESIF", "seismograms_all_stations_average_and_strong") # Seismograms
sta_lta_save_path = os.path.join("sta_lta", "JSON") # Trainwaves dictionaries

PRE_TRIGGER = 5 # in seconds
POST_TRIGGER = 20 # in seconds

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

#  all_seismogram_reduced = list_seismograms(mseeds_path)

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
    