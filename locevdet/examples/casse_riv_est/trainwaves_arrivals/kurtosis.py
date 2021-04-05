""" Kurtosis to specify more precisely start time of trainwaves"""

import os

from locevdet.trainwaves_arrivals.kurtosis import kurtosis_for_all_seismograms

mseeds_path = os.path.join("mseeds", "RESIF") # Seismograms
sta_lta_save_path = os.path.join("sta_lta", "JSON") # Trainwaves dictionaries

PRE_TRIGGER = 5 # in seconds
POST_TRIGGER = 5 # in seconds

kurtosis_for_all_seismograms(mseeds_path, sta_lta_save_path, PRE_TRIGGER, POST_TRIGGER)

