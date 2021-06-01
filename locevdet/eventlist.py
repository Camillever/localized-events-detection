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


    def set_matlab_data(self, matlab_folder_path, max_time_difference:float=2):
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

                    
                    if same_starttime < max_time_difference and \
                        trainwave.station.name == matlab_data['station'] :
                        mat = loadmat(matlab_filepath)
                        ti = matlab_data['starttime'] + mat["wavetrains"]['Tiwp'][0]
                        
                        matlab_data['trainwave'] = {
                            'radial_signal': mat["wavetrains"]['r'][0],
                            'initial_time': ti,
                            'fmin': mat["fmin"],
                            'centralfrequency': mat["wavetrains"]['centralfrequency'][0],
                            'duration': mat["wavetrains"]['duration'][0]
                        }
                        trainwave.matlab_data = matlab_data
                    else:
                        continue

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
        import matplotlib.dates as dates
        from sklearn.metrics import mean_squared_error, r2_score
        from locevdet.utils import linear_regression_on_dates
        from locevdet.visualisation.plots import demo_con_style

        plt.close("all")

        list_HIM_starts_specific, list_HIM_ti, diff_ti_HIM = [], [], []
        list_HIM_nsr = []

        list_FRE_starts_specific, list_FRE_ti,  diff_ti_FRE = [], [], []
        list_FRE_nsr = []

        fig_comp = plt.figure('01-02-2020_11-02-2020')
        if save_fig_path is not None:
            fig_comp.set_size_inches((20, 10), forward=False)

        for event in self:
            for _, trainwave in event.trainwaves.items():
                
                if trainwave.matlab_data is not None and trainwave.kurtosis_data is not None:
                    ti = trainwave.matlab_data['trainwave']['initial_time']
                    start_specific = trainwave.kurtosis_data["start_specific"]
                    if trainwave.station.name == 'HIM':
                        list_HIM_ti.append(ti)
                        list_HIM_starts_specific.append(UTCDateTime(start_specific)) # ELSE : Nan 
                        diff_ti_HIM.append(
                            np.abs(ti - start_specific))
                        # list_HIM_nsr.append(trainwave.snr)

                    elif trainwave.station.name == 'FRE':
                        list_FRE_ti.append(ti)
                        list_FRE_starts_specific.append(UTCDateTime(start_specific))  ## ELSE : Nan
                        diff_ti_FRE.append(
                            np.abs(ti - start_specific))
                        # list_FRE_nsr.append(trainwave.snr)
        print("list_HIM_ti :", list_HIM_ti)
        print("list_HIM_starts_specific :", list_HIM_starts_specific)
        print("diff_ti_HIM :", diff_ti_HIM)
        print("diff HIM shape :", len(diff_ti_HIM))
        
        list_HIM_ti_date64 = [np.datetime64(date) for date in list_HIM_ti]
        list_HIM_starts_specific_date64 = [np.datetime64(date) for date in list_HIM_starts_specific]
        list_FRE_ti_date64 = [np.datetime64(date) for date in list_FRE_ti]
        list_FRE_starts_specific_date64 = [np.datetime64(date) for date in list_FRE_starts_specific]
        

        myFmt = dates.DateFormatter('%m %d %H:%M')
        secondsFmt = dates.DateFormatter('%S')

        
        if type_graph != "dt":
            ax1 = fig_comp.add_subplot(121)
            ax1.set_title('HIM', fontsize=15, fontweight='bold')
        else:
            ax = plt.subplot(111)
            ax.set_title('Intervalle de temps d\'arrivées entre HIM et FRE', fontsize=15, fontweight='bold')
        
        if type_graph == "ti_with_snr":
            
            him_sc = ax1.scatter(list_HIM_ti_date64, list_HIM_starts_specific_date64)

            # Linear Regression
            X_HIM, Y_HIM, Y_HIM_pred, modeleRegHIM = linear_regression_on_dates(list_HIM_ti_date64, list_HIM_starts_specific_date64)
            ax1.plot(X_HIM, Y_HIM_pred, c='red')
            slopeHIM = format(modeleRegHIM.coef_[0][0], '.8f')
            interceptHIM = format(modeleRegHIM.intercept_[0], '.2f')
            mean_sq_errHIM_seconds = mean_squared_error(Y_HIM, Y_HIM_pred)*24*60*60  # in seconds
            mean_sq_errHIM = format(mean_sq_errHIM_seconds, '.2e')
            
            demo_con_style(
                ax1, 
                f"Coefficients (a ; b): ({slopeHIM} ; {interceptHIM}), Erreur quadratique moyenne: {mean_sq_errHIM} s", 
                colorborder='red'
            )
            # The coefficient of determination: 1 is perfect prediction
            print('Coefficient of determination: %.2f'
                % r2_score(Y_HIM, Y_HIM_pred))

            ax1.set_xlabel('Débuts de trains d\'ondes par méthode de Stockwell')
            ax1.set_ylabel('Débuts de trains d\'ondes par Kurtosis')
            ax1.xaxis.set_major_formatter(myFmt)
            ax1.yaxis.set_major_formatter(myFmt)

            # plt.colorbar(him_sc)
        
        elif type_graph == "hist_diff_ti":
            bins = np.arange(0, np.max(diff_ti_HIM), 0.5)
            ax1.hist(diff_ti_HIM, bins=bins, alpha=0.5)
            ax1.set_xlabel('abs( ti - start_specific)     (secondes)')
            ax1.set_ylabel('Nombre de trains d\'ondes')
            ax1.grid(True)
        
        elif type_graph == "dt":
            dt_stockwell = [
                list_HIM_ti[i] - list_FRE_ti[i]
                for i in range(len(list_HIM_ti))
            ]
            dt_kurtosis= [
                list_HIM_starts_specific[i] - list_FRE_starts_specific[i]
                for i in range(len(list_HIM_ti))]
            print("dt_stockwell", dt_stockwell)
            stockwell_label = 'Intervalle de temps entre HIM et FRE selon la méthode des spectogrammes de Stockwell'
            stem_stockwell = ax.stem(list_HIM_ti_date64, dt_stockwell, label=stockwell_label)
            kurtosis_label = 'Intervalle de temps entre HIM et FRE selon la méthode de Kurtosis'
            stem_kurtosis = ax.stem(list_HIM_ti_date64, dt_kurtosis, label=kurtosis_label)
            plt.setp(stem_stockwell, color=(1, 0, 0, 0.3))
            plt.setp(stem_kurtosis, color=(0, 0, 1, 0.3))
            # ax.set_yscale("log")
            ax.set_ylabel('dt entre HIM et FRE (secondes)')
            ax.xaxis.set_major_formatter(myFmt)
            ax.grid(True)
            plt.legend()

        if type_graph != "dt":
            ax2 = fig_comp.add_subplot(122)
            ax2.set_title('FRE', fontsize=15, fontweight='bold')

        if type_graph == "ti_with_snr":
            fre_sc = ax2.scatter(list_FRE_ti_date64, list_FRE_starts_specific_date64)

            # Linear Regression
            X_FRE, Y_FRE, Y_FRE_pred, modeleRegFRE = linear_regression_on_dates(list_FRE_ti_date64, list_FRE_starts_specific_date64)
            ax2.plot(X_FRE, Y_FRE_pred, c='red')
            slopeFRE = format(modeleRegFRE.coef_[0][0], '.8f')
            interceptFRE = format(modeleRegFRE.intercept_[0], '.2f')
            mean_sq_errFRE_seconds = mean_squared_error(Y_FRE, Y_FRE_pred)*24*60*60  # in seconds
            mean_sq_errFRE = format(mean_sq_errFRE_seconds, '.2e')
            demo_con_style(
                ax2, 
                f"Coefficients (a ; b): ({slopeFRE} ; {interceptFRE}), Erreur quadratique moyenne: {mean_sq_errFRE} s", 
                colorborder='red'
            )
            # The coefficient of determination: 1 is perfect prediction
            print('Coefficient of determination: %.2f'
                % r2_score(Y_FRE, Y_FRE_pred))

            ax2.set_xlabel('Débuts de trains d\'ondes par méthode de Stockwell')
            ax2.set_ylabel('Débuts de trains d\'ondes par Kurtosis')
            ax2.xaxis.set_major_formatter(myFmt)
            ax2.yaxis.set_major_formatter(myFmt)
        
        elif type_graph == "hist_diff_ti":
            bins = np.arange(0, np.max(diff_ti_FRE), 0.5)
            ax2.hist(diff_ti_FRE, bins=bins, alpha=0.5)
            ax2.set_xlabel('abs( ti - start_specific)     (secondes)')
            ax2.set_ylabel('Nombre de trains d\'ondes')
            ax2.grid(True)
            # fig_comp.colorbar(fre_sc)

        plt.gcf().autofmt_xdate()
        plt.tight_layout()

        if show is True:
            plt.show()

        if save_fig_path is not None:
            if type_graph == "ti_with_snr":
                figname = "compare_ti_01-02-2020_11-02-2020.png"
            elif type_graph == "diff_ti_in_fct_snr":
                figname = "diff_ti_fct_snr_01-02-2020_11-02-2020.png"
            elif type_graph == "hist_diff_ti":
                figname = "hist_diff_ti_01-02-2020_11-02-2020.png"
            elif type_graph == "dt":
                figname = "dt_01-02-2020_11-02-2020.png"
            fig_save_path = os.path.join(save_fig_path, figname)
            fig_comp.savefig(fig_save_path, bbox_inches='tight')
