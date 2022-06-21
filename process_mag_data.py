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
    for (path, layers, files) in os.walk(lattice_type):
        print("path = " + path)
        for layer in layers:
            print("layer =" + layer)
            for (path, temperatures, files) in os.walk(layer):
                for temperature in temperatures:
                    T = float(temperature.split("K"))

                    filename = lattice_type

                    return

plot_mag_vs_temperature("fcc_111")
