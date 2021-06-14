""" Module to define the class Trajectories
To analyse simulated trajectories simulated on Pierre98
"""
import os

import matplotlib.pyplot as plt

import numpy as np

from locevdet.utils import skipper

class Trajectories():
    def __init__(self, volume_boulder, **kwargs):
        self.volume_boulder = volume_boulder
        self.number_runs = kwargs.get('number_runs', 5)
        self.number_simul_per_run = kwargs.get('number_simul_per_run', 1000)
        self.runs = {}
    
    def add_runs(self, folder_path:str):
        """ Add all runs (and simulations) into arrays in the dictionary 'Trajectories.runs'
        
        Args:
            folder_path : Directory path containing the X runs folders
        
        Returns:
            Save all simulations of trajectories into Trajectories.runs
            
        """
        runs = ['run' +str(i) for i in range(1,self.number_runs+1)]
        for num, run in enumerate(runs):
            folder_path_run = os.path.join(folder_path, run)
            all_traj = [
                traj for traj in os.listdir(folder_path_run)
                if traj.endswith(".txt") and traj.startswith("Traj")
            ]
            all_traj.remove('Traj.txt')
            all_traj.remove('Traj_material.txt')
            all_data = []
            for traj in all_traj:
                data_onefile = np.loadtxt(skipper(os.path.join(folder_path, run, traj), header=True))
                all_data.append(data_onefile)
            self.runs[str(num)] = all_data
        
    def energy_fct_time(self, save_folder:str):
        """ Save the plot of seismic energy of simulated trajectories in function of time
        
        Args :
            save_folder : Directory path where plot will be saved
        Returns:
            Save the plot in the given save_folder
        """
        plt.close("all")

        fig = plt.figure("energy_in_fct_time")

        if save_folder is not None:
            fig.set_size_inches((20, 10), forward=False)
        time = []
        energy = []
        for run in range(self.number_runs):
            run_dict = self.runs[str(run)]
            for traj in range(self.number_simul_per_run):
                trajectory_dict = run_dict[traj]
                time += list(trajectory_dict[:,0])
                energy += list(trajectory_dict[:,6])


        plt.scatter(time, energy)
        plt.xlabel('Temps (s)')
        plt.ylabel('Energie (J)')

        title = (
            f" Volume : {self.volume_boulder} m$^{3}$ "
        )
        fig.suptitle(title, fontsize=18)
        plt.tight_layout()
        if save_folder is not None:
            figname = f"traj_{self.volume_boulder}.png"
            fig_save_path = os.path.join(save_folder, figname)
            fig.savefig(fig_save_path, bbox_inches='tight')
 
    def __repr__(self):
        return f"Trajectories - Volume {self.volume_boulder} m$^{3}$ \
            runs : {self.runs}"