""" Module for waveforme preprocessing and downloading form a distant client. """

import os
from typing import List, Tuple, Callable

import pandas as pd
from tqdm import tqdm
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

from locevdet.utils import clean_utc_str, ref_duration_to_time_window

def download_from_stations(client:Client, network:str, stations:List[str],
        events_times:List[Tuple[float]], process_waveforms:Callable=lambda x:x,
        location:str='*', channel:str='*', save_path:str=None, override:bool=False):
    """
    Download MSEED seismograms from the given client.

    Args:
        client: Name of url Client (see: URL_MAPPINGS from obspy.clients.fdsn.header)
        network: Name of the seismic network.
        stations: List of the seismic station's names.
        events_times: List of (start, end) times of events in UTCDateTime format.
        process_waveforms: Function to process waveforms
        location: Location identifiers
        channel: Name of the seismic channel associated to the previous station.
        save_path: Directory path of downloaded and processed MSEED files
        override: False per default to avoid override existing folder containing MSEED files

    Returns:
        Download and save processed MSEED files
        The saved files have the nomenclature : 'network_station_start-time_end-time'

    """

    if not isinstance(client, Client):
        client = Client(client)

    if save_path is None:
        save_path = os.getcwd()
    os.makedirs(save_path, exist_ok=True)

    for i, station in enumerate(stations):
        network_txt = network + " " * (3 - len(network))
        desc_msg = f"Network {network_txt} - Station {i+1}/{len(stations)} ({station})"
        for start_time, end_time in tqdm(events_times, desc=desc_msg, total=len(events_times)):
            mseed_filename = '_'.join((
                network, station,
                clean_utc_str(start_time),
                clean_utc_str(end_time)
            ))
            mseed_save_path = os.path.join(save_path, mseed_filename)

            if override or not os.path.isfile(mseed_save_path):
                waveforms = client.get_waveforms(
                    network, station, location, channel,
                    start_time, end_time, attach_response=True
                )
                waveforms = process_waveforms(waveforms)
                waveforms.write(mseed_save_path, format='MSEED', encoding='FLOAT64')

def read_events_times_from_csv(csv_events_path:str, csv_delimiter:str=';',
        time_offset:Tuple[int]=(0, 0)) -> List[Tuple[float]]:
    """
    Read time of each detected events from a given csv catalog
    and return the start and end time of window

    Args:
        csv_events_path: Dictionary path of the catalog of events in format csv.
        csv_delimiter: Delimiter type of columns in the given csv file.
        time_offset: Offset before and after the time window in the catalog.

    Returns:
        List of (start, end) times in UTCDateTime per events extracted from the catalog

    """
    csv_events = pd.read_csv(csv_events_path, delimiter=csv_delimiter)
    events_times = []
    for _, row in csv_events.iterrows():
        time_window = ref_duration_to_time_window(
            time_reference=UTCDateTime(row[0]),
            duration=float(row[2]),
            time_offset=time_offset
        )
        events_times.append(time_window)
    return events_times

if __name__ == "__main__":
    pass
