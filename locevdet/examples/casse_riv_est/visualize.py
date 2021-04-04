""" Visualisation : content of traces, plot of seismograms """

import os

from locevdet.visualisation.content_trace import content_traces


CLIENT = 'RESIF'
mseeds_client_path = os.path.join("mseeds", CLIENT)

# To see content of all seismograms MSEED
content_traces(folder_in=mseeds_client_path)

# Plot
## Per station

## Per event

## Per event with spectrograms