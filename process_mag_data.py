import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re

# use 4 space delimiter
sep = '    '

input_file = snakemake.input
output_file = snakemake.output
temperature = snakemake.params.get("temp")

mag_data = np.genfromtxt(input_file, skip_header=1, usecols=8)
mag_data = mag_data[200:mag_data.size]

mag = np.mean(mag_data)

with open(output_file, 'a') as f:
    w_str = str(temperature) + sep + str(mag)
    f.write(w_str)
    f.close()
