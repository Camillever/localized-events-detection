""" Kurtosis test - application"""

import os

from locevdet.trainwaves_arrivals.kurtosis_test import kurtosis_starts_extraction

mseeds_path = os.path.join("mseeds", "RESIF") # Seismograms
sta_lta_save_path = os.path.join("sta_lta", "JSON") # Trainwaves dictionaries

PRE_TRIGGER = 5 # in seconds
POST_TRIGGER = 5 # in seconds

win = 200 # in seconds

#Filter parameters
kurtosis_type = 'PANDAS_ROLLING' # 'OBPSY' or 'PANDAS_ROLLING'
THR_ON = 0.60

FREQMIN = 2
FREQMAX = 13

kurtosis_starts_extraction(
    mseeds_path, 
    sta_lta_save_path, 
    PRE_TRIGGER, 
    POST_TRIGGER, 
    thr_on=THR_ON, 
    freqmin=FREQMIN, 
    freqmax=FREQMAX, 
    win=win
    )

