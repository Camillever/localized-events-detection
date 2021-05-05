""" Class EventList """
import os
import json
import pickle

from typing import List
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as dates

from scipy import stats

from obspy import UTCDateTime

from mat4py import loadmat

from locevdet.event import Event

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


    def set_matlab_data(self, matlab_folder_path, max_time_difference=5):
        matlab_filepaths = [
            os.path.join(matlab_folder_path, filename)
            for filename in os.listdir(matlab_folder_path)
            if filename.endswith('.mat')
        ]

        for matlab_filepath in matlab_filepaths:
            
            title_data = matlab_filepath.split('_')
            
            matlab_data = {
                'network': title_data[3],
                'station': title_data[4],
                'starttime': UTCDateTime(title_data[5]),
                'endtime': UTCDateTime(title_data[6].split('.')[0]),
            }
            matlab_data['periodtime'] = "_".join((str(matlab_data['starttime']), str(matlab_data['endtime'])))
            matlab_starttime = UTCDateTime(matlab_data['starttime'])

            for event in self:
                for _, trainwave in event.trainwaves.items():
                    starttime = UTCDateTime(trainwave.trace.stats.starttime)
                    same_starttime = abs(starttime - matlab_starttime) 
                    if same_starttime < max_time_difference and trainwave.station.name == matlab_data['station']:
                        mat = loadmat(matlab_filepath)
                        matlab_data['trainwave'] = {
                            'radial_signal': mat["wavetrains"]['r'][0],
                            'initial_time': matlab_data['starttime'] + mat["wavetrains"]['Tiwp'][0],
                            'fmin': mat["fmin"],
                            'centralfrequency': mat["wavetrains"]['centralfrequency'][0],
                            'duration': mat["wavetrains"]['duration'][0]
                        }
                        trainwave.matlab_data = matlab_data

    def plots_time_compare_HIM_FRE(self, type_graph:str, 
            save_fig_path:str=None, show:bool=True):
        """ TODO
    
    Args:
        type_graph : Type of graphique to plot as     
            - "ti_with_snr" : Initials times (with snr values) of trainwaves determined by
                    kurtosis in function of initials times determined by Stockwell's spectograms
            - " diff_ti_in_fct_snr" : Difference of initials times determined by kurtosis
                    and Stockwell's spectograms in function of snr values (for ti by kurtosis)
            - "hist_diff_ti" : Histograms of the difference of initials times determined by
                    kurtosis and Stockwell's spectograms
             - "dt" : delta time between initials times of trainwaves between the two stations.
                    Comparaison between dt found by Stockwell's spectograms and kurtosis.

        save_fig_path : Directory path of the folder where to save the figure.
                If None, the figure will be not saved.

        show : True, show the figure in a matplotlib window. False, the figure will be not shown.

    Returns:
        TODO
        """
        
        
        plt.close("all")

        list_HIM_starts_specific =[]
        list_HIM_ti = []
        list_HIM_nsr = []

        list_FRE_starts_specific =[]
        list_FRE_ti = []
        list_FRE_nsr = []

        for event in self:
            for _, trainwave in event.trainwaves.items():
                if trainwave.matlab_data is not None :

                    if trainwave.station.name == 'HIM':
                        list_HIM_ti.append(np.datetime64(trainwave.matlab_data['trainwave']['initial_time']))
                        list_HIM_starts_specific.append(np.datetime64(trainwave.kurtosis_data["start_specific"]))
                        # list_HIM_nsr.append(trainwave.snr)

                    elif trainwave.station.name == 'FRE':
                        list_FRE_ti.append(np.datetime64(trainwave.matlab_data['trainwave']['initial_time']))
                        list_FRE_starts_specific.append(np.datetime64(trainwave.kurtosis_data["start_specific"]))
                        # list_FRE_nsr.append(trainwave.snr)
        
        fig_comp = plt.figure('01-02-2020_11-02-2020')
        dateformat = dates.DateFormatter('%m %d %H:%M')

        if save_fig_path is not None:
            fig_comp.set_size_inches((20, 10), forward=False)

        ax1 = fig_comp.add_subplot(121)
        ax1.set_title('HIM', fontsize=15, fontweight='bold')
        him_sc = ax1.scatter(list_HIM_ti, list_HIM_starts_specific)

        slope, intercept, r_value, p_value, std_err = stats.linregress(list_HIM_ti, list_HIM_starts_specific)
        ax1.plt(list_HIM_ti, slope * list_HIM_ti + intercept)

        ax1.set_xlabel('Débuts de trains d\'ondes par méthode de Stockwell')
        ax1.set_ylabel('Débuts de trains d\'ondes par Kurtosis')
        ax1.xaxis.set_major_formatter(dateformat)
        ax1.yaxis.set_major_formatter(dateformat)

        # plt.colorbar(him_sc)

        ax2 = fig_comp.add_subplot(122)
        ax2.set_title('FRE', fontsize=15, fontweight='bold')
        fre_sc = ax2.scatter(list_FRE_ti, list_FRE_starts_specific)
        ax2.set_xlabel('Débuts de trains d\'ondes par méthode de Stockwell')
        ax2.set_ylabel('Débuts de trains d\'ondes par Kurtosis')
        ax2.xaxis.set_major_formatter(dateformat)
        ax2.yaxis.set_major_formatter(dateformat)

        # fig_comp.colorbar(fre_sc)
        plt.gcf().autofmt_xdate()
        plt.tight_layout()

        if show is True:
            plt.show()

        if save_fig_path is not None:
            figname = "compare_ti_01-02-2020_11-02-2020.png"
            fig_save_path = os.path.join(save_fig_path, figname)
            fig_comp.savefig(fig_save_path, bbox_inches='tight')
