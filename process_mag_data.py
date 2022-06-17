import numpy
import scipy
import pandas as pd
import matplotlib

def plot_mag_vs_timesteps(filename):
    data = pd.read_csv(filename, sep='')
