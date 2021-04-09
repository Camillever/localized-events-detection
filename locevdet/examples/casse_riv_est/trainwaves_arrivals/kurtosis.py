""" Kurtosis to specify more precisely start time of trainwaves"""

import os

from locevdet.trainwaves_arrivals.kurtosis import kurtosis_for_all_seismograms

mseeds_path = os.path.join("mseeds", "RESIF") # Seismograms
sta_lta_save_path = os.path.join("sta_lta", "JSON") # Trainwaves dictionaries

PRE_TRIGGER = 5 # in seconds
POST_TRIGGER = 5 # in seconds
THR_ON = 0.18 # in purcentage of the maximum value of kurtosis
THR_OFF = 0.1 # in purcentage of the maximum value of kurtosis

win = 0.03 # in seconds

#Filter parameters
FREQMIN=3
FREQMAX=8

kurtosis_for_all_seismograms(mseeds_path, sta_lta_save_path, PRE_TRIGGER, POST_TRIGGER, THR_ON, THR_OFF, FREQMIN, FREQMAX, win)

