import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re

def plot_mag_vs_timesteps(filename):
    data = pd.read_csv(filename, sep=r"[ \t]+") 
    data.info()

    d = {'Time': data['time'], 'm_z': data['Fe_m']}
    df = pd.DataFrame(data=d)

    fig, ax = plt.subplots()

    df.plot('Time', 'm_z', 'scatter', ax)
    plt.savefig('mag_vs_timesteps.png')

def plot_mag_vs_temperature(lattice_type, layers):
    for layer in layers:
        print(layer)
        _dir = lattice_type + "/" + str(layer)
            
        print("Walking " + _dir)

        for path, dirnames, filenames in os.walk(_dir):
            print(path)
            true = bool(re.search(r"^.*K.*$", path))
            if true:
                print("i have entered this if statement, hallelujah")

plot_mag_vs_temperature("fcc_111", [2,3,4,5,6])
