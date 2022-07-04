import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re
import pymbar

TEMP_REGEX = r'[+-]?([0-9]*[.])[0]+'

# use 4 space delimiter
sep = '\t'

in_files = snakemake.input
out_file = str(snakemake.output)

with open(out_file, 'a') as f:
    for in_file in in_files: 
        print(in_file)
        # TODO maybe avoid hardcoding the column number for elegance
        mag_data = np.genfromtxt(in_file, skip_header=1, usecols=7) # m_z
        mag = np.mean(mag_data)

        m = re.search(TEMP_REGEX, in_file)
        temp = float(in_file[m.start():m.end()])
        print(str(temp))

        w = str(temp) + sep + str(round(mag, 3)) + "\n"
        f.write(w)
    f.close()
