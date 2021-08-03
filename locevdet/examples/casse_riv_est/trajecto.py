""" Script to analyse rockfall's trajectories simulations """
import os

from locevdet.trajectories import Trajectories

number_runs = 5
simulations = 1000

volume_type_foldername = [("Traj_5m3", 5), ("Traj_10m3_Distrib_uniform_15t-50t", 10)]

for (volume_folder, volume) in volume_type_foldername:
    trajecto_folder = os.path.join('documents', 'trajecto', 'Trajecto_Pierre98', volume_folder)
    traj_5m3 = Trajectories(volume_boulder=volume, number_runs=number_runs, number_simul_per_run=simulations)
    traj_5m3.add_runs(folder_path=trajecto_folder)

    fig_5m3_save_folder = os.path.join("captures", "trajectories", "energy_fct_time")
    traj_5m3.energy_fct_time(save_folder=fig_5m3_save_folder)
