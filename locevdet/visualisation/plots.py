""" Module for plotting waveforms """

import os
import matplotlib.pyplot as plt

from obspy import read
from tqdm import tqdm



def capture_one_event(folder_in, folder_out, start_time):
    """
    Returns a PNG image of one event started from the given start_time from all seismograms (MSEED) of this event into the given folder
    
    Args :
        folder_in : directory path of downloaded and processed MSEED files (str)
        folder_out : directory path of folder where the png image will be saved
        start_time : UTC start-time of the record of the event (cf nomenclatures of MSEED seismograms)
    Returns :
        A PNG image with all seismograms of this event with the nomenclature : start-time.png
        
    """
    all_seismogram = [
        filename for filename in os.listdir(folder_in)
        if filename.startswith("PF_") or filename.startswith("G_RER")
    ]
    start_time_name = str(start_time)[:-8]
    
    fig = plt.figure(start_time_name, figsize = [25, 15])
    filenames_event = []
    for filename in all_seismogram :
        
        fig.suptitle(filename.split('_')[2]+' --- '+filename.split('_')[3])
        
        if start_time_name.replace(':', '-') == filename.split('_')[2]:
            filepath = os.path.join(folder_in, filename)
            filenames_event.append(filename)
            seismogram = read(filepath)
            
            if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_PER"):
                trace = seismogram[0]
                
            elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
                trace = seismogram[0]
                trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
                
            else:
                trace = seismogram[2]
            
            
            ax = fig.add_subplot(len(stations), 1, len(filenames_event))
            ax.plot(trace.times(), trace.data)
            ax.set_title(filename.split('_')[0]+'_'+filename.split('_')[1], loc='right', fontweight ="bold")
            
            plt.tight_layout()
    
    fig.savefig(os.path.join(folder_out, start_time_name.replace(':', '-')+'.png'))
    
    plt.show()
    
    return (print('Figure saved into :', folder_out, 'with the name :', start_time_name.replace(':', '-')+'.png'))




# Par station
def captures_per_station(folder_in:str, folder_out:str, all_seismogram, networks_stations, n_cols, n_lines):
    """
    
    
    """

    for network_station in networks_stations :
        i=0
        for seismogram_number, filename in enumerate(all_seismogram):
            if filename.startswith(network_station):
                i+=1

                filepath = os.path.join(folder, filename)
                seismogram = read(filepath)


                if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_PER"):
                    trace = seismogram[0] # Une seule composante et elle est verticale
                elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
                    trace = seismogram[0]
                    trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
                else:
                    trace = seismogram[2] # Composant verticale

                dt = trace.stats.starttime.timestamp
                t = np.linspace(trace.stats.starttime.timestamp - dt,
                                trace.stats.endtime.timestamp - dt,
                                trace.stats.npts)

                fig = plt.figure(network_station, figsize = [30, 20])
                plt.title(network_station)
                plt.subplot(n_lines, n_cols, i)
                plt.plot(trace.times(), trace.data)
                plt.tight_layout()

        fig.savefig(os.path.join(folder_out, network_station+'_processed.png'))


n_cols = 8  # for 72 events
n_lines = 9

networks_stations = ['PF_FRE', 'PF_HIM', 'PF_NSR', 'PF_PER', 'PF_TTR', 'G_RER']
folder_path = os.path.join(path, 'mseed_RESIF')
folder_out_per_station = os.path.join(folder_path, 'captures','per_station')

all_seismogram = [
        filename for filename in os.listdir(folder_path)
        if filename.startswith("PF_") or filename.startswith("G_RER")
    ]

# captures_per_station(folder_path, folder_out_per_station, all_seismogram, networks_stations, n_cols, n_lines)


#Par événement 

def get_starttime(filename: str)-> str:
    post_name = filename.split('_')[2]
    return post_name



def captures_per_event(folder, folder_out, all_seismogram, networks_stations):
    starttime_filename = []

    HIM_filenames = [
        filename for filename in os.listdir(folder)
        if filename.startswith('PF_HIM') 
    ]

    for starttime in HIM_filenames:
        starttime_filename.append(get_starttime(starttime))


    for start_name in tqdm(starttime_filename):
        filenames_event = []
        fig = plt.figure(start_name, figsize = [25, 15])

        for filename in all_seismogram :

            if start_name == filename.split('_')[2]:
                fig.suptitle(filename.split('_')[2]+' --- '+filename.split('_')[3]) 
                filepath = os.path.join(folder, filename)
                filenames_event.append(filename)

                seismogram = read(filepath)


                if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_TTR") or filename.startswith("PF_PER"):
                    trace = seismogram[0]
                elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
                    trace = seismogram[0]
                    trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
                else:
                    trace = seismogram[2]

                dt = trace.stats.starttime.timestamp
                t = np.linspace(trace.stats.starttime.timestamp - dt,
                                        trace.stats.endtime.timestamp - dt,
                                        trace.stats.npts)

                ax = fig.add_subplot(len(networks_stations), 1, len(filenames_event))

                ax.plot(trace.times(), trace.data)
                ax.set_title(filename.split('_')[0]+'_'+filename.split('_')[1], loc='right', fontweight ="bold")
                plt.tight_layout()

        fig.savefig(os.path.join(folder_out, start_name+'_processed.png'))
#         plt.show()

networks_stations = ['PF_FRE', 'PF_HIM', 'PF_NSR', 'PF_PER', 'PF_TTR', 'G_RER']
folder_path = os.path.join(path, 'mseed_RESIF')
folder_out_per_event = os.path.join(folder_path, 'captures', 'per_event')

all_seismogram = [
        filename for filename in os.listdir(folder_path)
        if filename.startswith("PF_") or filename.startswith("G_RER")
    ]

# captures_per_event(folder_path, folder_out_per_event, all_seismogram, networks_stations)


# Par événement (et avec spectogrammes)
import math
from scipy import signal

def nextpow2(i):
    """
    Find 2^n that is equal to or greater than.
    """
    
    n = 1
    while n < i:
        n *= 2
        
    return n

def captures_per_event_spectogramms(folder, folder_out, all_seismogram, networks_stations):
    
    starttime_filename = []
    HIM_filenames = [
        filename for filename in os.listdir(folder)
        if filename.startswith('PF_HIM') 
    ]

    for starttime in HIM_filenames:
        starttime_filename.append(get_starttime(starttime))
    
    
    for start_name in tqdm(starttime_filename):
        filenames_event = []
        fig = plt.figure(start_name, figsize = [25, 15])

        for filename in all_seismogram :

            if start_name == filename.split('_')[2]:
                fig.suptitle(filename.split('_')[2]+' --- '+filename.split('_')[3]) 

                filepath = os.path.join(folder, filename)
                filenames_event.append(filename)

                seismogram = read(filepath)

                if filename.startswith("PF_RER") or filename.startswith("PF_NSR") or filename.startswith("PF_TTR") or filename.startswith("PF_PER"):
                    trace = seismogram[0]
                elif filename.startswith("G_RER") or filename.startswith("PF_TTR"):
                    trace = seismogram[0]
                    trace = trace.filter('bandpass', freqmin=1.5, freqmax=20)  # High-pass filter for these stations
                else:
                    trace = seismogram[2]

                dt = trace.stats.starttime.timestamp
                t = np.linspace(trace.stats.starttime.timestamp - dt,
                                        trace.stats.endtime.timestamp - dt,
                                        trace.stats.npts)

                # Seismograms
                ax1 = fig.add_subplot(len(networks_stations), 2, len(filenames_event))


                ax1.plot(trace.times(), trace.data, 'black')

                ax1.set_title(filename.split('_')[0]+'_'+filename.split('_')[1], loc='right', fontweight ="bold")
                ax1.set_xlabel('Temps (s)')
                ax1.set_ylabel('Vitesse (m/s)')

                # Spectograms
                window_nb_ech = int(5*trace.stats.sampling_rate)
                noverlap = int(math.floor(window_nb_ech*90/100))
                nfft = int(nextpow2(window_nb_ech))
                signal_windowed = np.ones(int(nfft))
                window = signal_windowed *signal.tukey(len(signal_windowed), 0.05)
                number_pts_padded = int(nfft)

                filenames_event.append('spectro')
                print(len(filenames_event))
                ax2 = fig.add_subplot(len(networks_stations), 2, len(filenames_event))

                f, t, sxx = signal.spectrogram(trace.data, trace.stats.sampling_rate, window=window, noverlap=noverlap, nfft=nfft, mode='psd')
                sxx_log10_normalized = (np.log10(sxx) - np.min(np.log10(sxx)))/(np.max(np.log10(sxx))-np.min(np.log10(sxx)))
                spectro = plt.pcolor(t, f, sxx_log10_normalized>0.9 , cmap='magma')


                ax2.set_xlabel('Temps (s)')
                ax2.set_ylabel('Fréquence (Hz)')
                if filename.startswith("G_RER") :
                    ax2.set_ylim([0.5, 10])
                else : 
                    ax2.set_ylim([0.5, 25])

                fig.colorbar(spectro, ax=ax2)
                plt.tight_layout()


        fig.savefig(os.path.join(folder_out, start_name+'_processed.png'))
        plt.show() 


folder_path = os.path.join(path, 'mseed_RESIF')
folder_out_per_event_spectograms = os.path.join(folder_path, 'captures', 'per_event_with_spectrogram', 'max_energy')
networks_stations = ['PF_FRE', 'PF_HIM', 'PF_NSR', 'PF_PER', 'PF_TTR', 'G_RER']

all_seismogram = [
        filename for filename in os.listdir(folder_path)
        if filename.startswith("PF_") or filename.startswith("G_RER")
    ]

# captures_per_event_spectogramms(folder_path, folder_out_per_event_spectograms, all_seismogram, networks_stations)