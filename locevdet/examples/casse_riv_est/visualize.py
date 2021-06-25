""" Visualisation : content of traces, plot of seismograms """

import os

from obspy import UTCDateTime

from locevdet.utils import get_info_from_mseedname
from locevdet.visualisation.content_trace import content_traces
from locevdet.visualisation.plots import captures_per_event
from locevdet.examples.casse_riv_est.all_seismo import list_seismograms

CLIENT = 'RESIF'
# ###################################### 01-02-2020_11-02-2020 #####################################
# # ROCKFALL 
# mseeds_rockfall_path = \
#     os.path.join("seismograms", "01-02-2020_11-02-2020", "rockfall", "average_strong")
# all_seismograms = list_seismograms(mseeds_rockfall_path)

# ## To see content of all seismograms MSEED
# content_traces(all_seismogram=all_seismograms, folder_in=mseeds_rockfall_path)

# ## Plot
# ### Per station

# ### Per event
# save_fig_rockfall_event_path = \
#     os.path.join("captures", "01-02-2020_11-02-2020", "rockfall", "seismograms", "per_event")
# captures_per_event(mseeds_rockfall_path, save_fig_path=save_fig_rockfall_event_path)

# ### Per event with spectrograms


# # VOLCANO-TECTONIC
# mseeds_volcanotectonic_path = \
#     os.path.join("seismograms", "01-02-2020_11-02-2020", "volcanotectonic")
# ## Content of seismograms MSEED
# mseeds_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "volcanotectonic")
# all_volcanotecto = [
#         filename for filename in os.listdir(mseeds_path)
#         if filename.startswith('PF_FRE') 
#         or filename.startswith('PF_HIM') 
#         or filename.startswith('G_RER')
# ]
# print(all_volcanotecto)
# content_traces(all_seismogram=all_volcanotecto, folder_in=mseeds_volcanotectonic_path)
# Plot
## Per Station

# ## Per event
# save_fig_volcanotectonic_event_path = \
#     os.path.join("captures", "01-02-2020_11-02-2020", "volcanotectonic", "seismograms", "per_event")
# captures_per_event(mseeds_volcanotectonic_path, save_fig_path=save_fig_volcanotectonic_event_path)

### Per event with spectograms


# # LOCALS SEISMS
# mseeds_locals_path = os.path.join("seismograms", "01-02-2020_11-02-2020", "local_seisms")

# ## Content of seismograms MSEED
# all_locals = os.listdir(mseeds_locals_path)
# content_traces(all_seismogram=all_locals, folder_in=mseeds_locals_path)

# ## Plot
# ### Per Station

# ### Per event
# save_fig_locals_event_path = \
#     os.path.join("captures", "01-02-2020_11-02-2020", "local_seisms", "seismograms", "per_event")
# captures_per_event(mseeds_locals_path, save_fig_path=save_fig_locals_event_path)
# ### Per event with spectograms


# CONTINUOUS
mseeds_path_continuous = os.path.join("seismograms","01-02-2020_11-02-2020","continuous")
## Content of seismograms MSEED

date_07fev = UTCDateTime("2020-02-07T00:00:00.00")
before_07fev_seismograms =[ 
    mseed for mseed in os.listdir(mseeds_path_continuous)
    if get_info_from_mseedname(mseed)['starttime'] < date_07fev
    if not mseed.startswith('PF_NSR') or not mseed.startswith('G_RER')
    ]
after_07fev_seismograms =[ 
    mseed for mseed in os.listdir(mseeds_path_continuous)
    if get_info_from_mseedname(mseed)['starttime'] > date_07fev
    if not mseed.startswith('G_RER')
    ]
continuous_seismograms = [*before_07fev_seismograms, *after_07fev_seismograms]

content_traces(all_seismogram=continuous_seismograms, folder_in=mseeds_path_continuous)