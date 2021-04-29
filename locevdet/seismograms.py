import os
import json
from typing import List

from obspy.core import Stream

from locevdet.stations import Station

class Trace():  
    def __init__(self, station:Station, **kwargs):

        self.starttime =  # in UTCDatetime 
        self.endtime =    # in UTCDatetime
        self.station = station
        self.signal =    # 
        self.radial_signal = kwargs.get('radial', None)  # Extracted from Kristel's results (Default None)
        

    def __repr__(self):
        return f"{self.startime} - {self.endtime} ({self.station})"


# All seismograms compilation
class Seismograms():  
    def __init__(self, )