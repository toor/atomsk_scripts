import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt

def plot_mag_vs_timesteps(filename):
    data = pd.read_csv(filename, sep=r"[ \t]+") 
    data.info()

    d = {'Time': data['time'], 'm_z': data['Fe_mz']}
    df = pd.DataFrame(data=d)

    fig, ax = plt.subplots()

    df.plot('Time', 'm_z', 'scatter', ax)
    plt.savefig('mag_vs_timesteps.png')

def plot_mag_vs_temperature(lattice_type):
    # list of paths to layers
    (path, layers, files) = os.walk(lattice_type)
    for layer in layers:
        print(layer)
        (path, temperatures, files) = os.walk(layer)
        for T in temperatures:
            # strip Kelvin suffix
            temp = float(T.strip("K"))

            filename = lattice_type + "/" + layer + "/" + T + "/jams_mag.tsv"
            
            # TODO slice data after equilibration
            data = pd.read_csv(filename, sep=r"[ \t]+")['Fe_mz']

            return

plot_mag_vs_temperature("fcc_111")
