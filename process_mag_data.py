import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re

def plot_mag_vs_timesteps(filename):
    data = pd.read_csv(filename, sep=r"[ \t]+") 
    #data.info()

    d = data['Fe_m']
    df = pd.DataFrame(data=d)
    arr = pd.to_numpy()

#    fig, ax = plt.subplots()

#    df.plot('Time', 'm_z', 'scatter', ax)
#    plt.savefig('mag_vs_timesteps.png')

def _plot_mag_vs_timesteps(filename):
    with open(filename) as f:
        lines = np.array([line for line in f])
        i = re.findall(lines[0], "Fe_m")[0].start()

def plot_mag_vs_temperature(lattice_type, layers):
    for layer in layers:
        print(layer)
        _dir = lattice_type + "/" + str(layer)

        temps_layer = []
        for path, dirnames, filenames in os.walk(_dir):
            if bool(re.search(r"^.*K.*$", path)):
                for m in re.finditer(r"[+-]?([0-9]*[.])?[0-9]+", path):
                    temp = path[m.start():m.end()]
                    temps_layer.append(float(temp))

filename = "fcc_111/3/50.0K/jams_mag.tsv"

_plot_mag_vs_timesteps(filename)
#plot_mag_vs_temperature("fcc_111", [2,3,4,5,6])
