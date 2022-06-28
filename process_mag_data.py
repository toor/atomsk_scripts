import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re

# use 4 space delimiter
sep = '    '

#output_file = snakemake.output
lattice = snakemake.params[0]
layers = snakemake.params[1]
temp = snakemake.params[2]

output_file = lattice + "/" + str(layers) + "mag_vs_temp.dat"

with open(output_file, 'a') as f:
    mag_data = np.genfromtxt(input_file, skip_header=1, usecols=8)
    mag = np.mean(mag_data[200:(mag_data.size - 1)])

    w_str = str(temp) + sep + str(mag)

    f.write(w_str)
    f.close()
