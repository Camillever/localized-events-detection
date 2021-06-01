""" Module to visualize content of trace or stream of traces"""

import os

from typing import List

from obspy import read

def content_traces(all_seismogram:List[str], folder_in:str):
    """
    To visualize the content of all seismograms in the given folder
    Give few details of traces contained into these seismograms.
    
    Args:
        all_seismogram : List of seismogram's filenames 
        folder_in: Directory path of MSEED seismograms
    Returns :
        print all the traces contained into each seismograms of the given folder
        and return the number of seismograms in this folder
    """
    for filename in all_seismogram:
        filepath = os.path.join(folder_in, filename)
        seismogram = read(filepath)
        print(seismogram)
    return (f'Nombre de sismogrammes : {len(all_seismogram)}')
