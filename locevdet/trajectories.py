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
        self.runs_with_bounces = {}
    
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
            all_data_reduced = []
            for traj in all_traj:
                data_onefile_with_bounces = np.loadtxt(skipper(os.path.join(folder_path, run, traj), header=True))

                ## Remove row with bounces 
                y_ground = data_onefile_with_bounces[:, 2]
                y_bounce = data_onefile_with_bounces[:, 3]
                data_onefile = data_onefile_with_bounces[~(y_ground == y_bounce)]
                
                # for row in range(np.size(data_onefile,0)):
                #     y_ground = data_onefile[row, 2]
                #     y_bounce = data_onefile[row, 3]
                #     if y_ground == y_bounce:

                
                all_data.append(data_onefile_with_bounces)
                all_data_reduced.append(data_onefile)
            self.runs[str(num)] = all_data_reduced
            self.runs_with_bounces[str(num)] = all_data
        
    def energy_fct_time(self, save_folder:str):
        """ Save the plot of seismic energy of simulated trajectories in function of time
        
        Args :
            save_folder : Directory path where plot will be saved
        Returns:
            Save the plot in the given save_folder
        """
        plt.close("all")

        fig = plt.figure("energy_in_fct_time")

        if self.volume_boulder == 5 :
            masse = 17000
        elif self.volume_boulder == 10:
            masse = 33500

        g_pesant = 9.81

        if save_folder is not None:
            fig.set_size_inches((20, 10), forward=False)
        time = []
        energy = []
        y_ground = []
        for run in range(self.number_runs):
            run_dict = self.runs[str(run)]
            for traj in range(self.number_simul_per_run):
                trajectory_dict = run_dict[traj]
                y_ground += list(trajectory_dict[:,2])
                time += list(trajectory_dict[:,0])

        for _, y in enumerate(y_ground):
            energy.append(masse*g_pesant*y)


        # plt.bar(time, energy, width=1, alpha=0.5)
        plt.plot(time, energy, '.')
        plt.yscale("log")
        plt.xlabel('Temps (s)')
        plt.ylabel('Energie (J)')

        plt.grid(True)

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