""" Classes Trainwave / Event / Eventlist"""

import os

import json
import pickle

from typing import List

from obspy.core import Stream

from locevdet.stations import Station
from locevdet.utils import clean_utc_str

class Trainwave():

    def __init__(self, trace, station:Station, **kwargs):
        self.trace = trace
        self.station = station

        self.start_specific = kwargs.get('start_specific', None)
        self.kurtosis_params = kwargs.get('kurtosis_params', None)

    def get_data_from_matlab_file(self, matlab_filepath):
        raise NotImplementedError

    def __repr__(self):
        return f"Trainwave{self.trace}"

class Event():

    def __init__(self, start_global, stations:List[Station]):
        self.start_global = start_global
        self.stations = stations
        self.trainwaves = {}

    def add_trainwaves(self, stream:Stream):
        stream_stations = [trace.stats.station for trace in stream]
        for station in self.stations:
            station_index = stream_stations.index(station.name)
            trace = stream[station_index]
            self.trainwaves[station] = Trainwave(trace, station)
    

    def to_json(self):
        return {
            'start_global': str(self.start_global),
            'stations': str(self.stations)
        }

    def save(self, filepath, format_save, override=False):
        if override or not os.path.isfile(filepath):
            if format_save == 'JSON':
                with open(filepath, 'w') as content:
                    json.dump(self.to_json(), content,  indent=1)

    def __eq__(self, other):
        if isinstance(other, Event):
            return self.start_global == other.start_global
        else:
            raise TypeError("Only two Events can be compared")

    def __str__(self):
        return clean_utc_str(self.start_global)

    def __repr__(self):
        return f"Event({str(self)}, {self.trainwaves})"

class EventList(list):

    def __init__(self, events:List[Event]=None):
        events = events if events is not None else []
        super().__init__(events)
        if events is None:
            self.events = []
        else:
            for event in events:
                if not isinstance(event, Event):
                    raise TypeError('EventList must be composed of Event objects')
            self.events = events

    def to_json(self):
        return {str(event): event.to_json() for event in self.events}

    def save(self, filepath, format_save='PICKLE', override=False):
        if override or not os.path.isfile(filepath):
            if format_save == 'JSON':
                with open(filepath, 'w') as content:
                    json.dump(self.to_json(), content,  indent=1)
            elif format_save == 'PICKLE':
                with open(filepath, 'wb') as content:
                    pickle.dump(self, content)
