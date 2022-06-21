import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt

def plot_mag_vs_timesteps(filename):
    data = pd.read_csv(filename)

    d = {'Time': data['Time'], 'm_z': data['Fe_mz']}
    df = pd.DataFrame(data=d)

    fig, ax = plt.subplots()

    df.plot('Time', 'm_z', 'scatter', ax)


filename = 'bcc_100/6/480.0K/jams_mag.tsv'

plot_mag_vs_timesteps(filename)
