"""Module to compute the localisation map"""

import os
import gdal
import numpy as np

def tif_to_ascii(inputname, outputname):
    gdal_translate -of "ascii" inputname.tif outputname.ascii




tif_path = os.path.join('topo', 'mnt')
tifname = 'MNT_Grille_Riviere_de_lEst'
topo_tif = os.path.join(tif_path, tifname)

tif_to_ascii(topo_tif, tifname)