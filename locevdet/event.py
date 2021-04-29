import os
import json
from typing import List

from obspy.core import Stream

from locevdet.stations import Station
from locevdet.seismograms import Trace
from locevdet.utils import clean_utc_str

class Trainwave():

    def __init__(self, signal:Trace, station:Station, **kwargs):
        self.signal_associated = signal
        self.station = station

        self.start_specific = kwargs.get('start_specific', None)
        self.kurtosis_params = kwargs.get('kurtosis_params', None)

    def __repr__(self):
        return f"{self.signal.shape} ({self.station})"

class Event():

    def __init__(self, start_global, stations:List[Station]):
        self.start_global = start_global
        self.stations = stations
        self.trainwaves = {}

    def add_trainwaves(self, stream:Stream):
        stream_stations = [trace.stats.station for trace in stream]
        for station in self.stations:
            station_index = stream_stations.index(station.name)
            signal = stream[station_index].data
            self.trainwaves[station] = Trainwave(signal, station)

    def save(self, filepath, format_save, override=False):
        if override or not os.path.isfile(filepath):
            if format_save == 'JSON':
                jsonable_event = {
                    'start_global': self.start_global,
                    'stations': self.stations
                }
                with open(filepath, 'w') as content:
                    json.dump(jsonable_event, content,  indent=1)
            
            # elif format_save == 'PICKLE':

            #     picklable_event = 


            #     with open(filepath, 'w') as content:
            #         pickle.dump(picklable_event, content)

    def __repr__(self):
        return f"Start global: {clean_utc_str(self.start_global)} {self.trainwaves}"

class EventList():

    def __init__(self, events:List[Event]=None):
        if events is None:
            self.events = []
        else:
            for event in events:
                if not isinstance(event, Event):
                    raise TypeError('EventList must be composed of Event objects')
            self.events = events
    
    def append(self, event:Event):
        self.events.append(event) 
    
    def save(self, filepath, format_save='PICKLE', override=False):
        self.save = Event.save(self, filepath, format_save, override)  # Save each event , doesn't it ? How save the List of events? 



