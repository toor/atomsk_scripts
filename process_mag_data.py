import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re

TEMP_REGEX = r'[+-]?([0-9]*[.])[0-9]+'

# use 4 space delimiter

sep = '\t'

in_files = snakemake.input
out_file = snakemake.output

with open(str(out_file), 'a') as f:
    for in_file in in_files:
        # Skip the first 200 steps as equilibration. TODO: Use pymbar.
        mag_data = np.genfromtxt(in_file, skip_header=201, usecols=7)
        mag = np.mean(mag_data)

        m = re.search(TEMP_REGEX, in_file)
        temp = float(in_file[m.start():m.end()])

        w = str(temp) + sep + str(round(mag, 3)) + "\n"
        f.write(w)
    f.close()
