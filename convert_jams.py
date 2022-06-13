# This converts an XSF file into a JAMS config file. Options should be added
# to specify the exchange properties and various other options, ideally at a top level
# i.e. a Snakefile to allow easy modification of these options.

import numpy as np
import libconf


def read_unitcell(filename):
    separator = '        '

    # this function should read the correct parameters for a unit cell
    with open(filename) as f:
        lines = [line for line in f]
        for j in range(0, len(lines)):
            if lines[j].strip() == "CRYSTAL":
                # start of file; primitive vectors appear immediately
                # after the CRYSTAL declaration
                a1 = lines[j+2].strip().split(separator)
                a2 = lines[j+3].strip().split(separator)
                a3 = lines[j+4].strip().split(separator)
 
                a1 = [float(a1[0]), float(a1[1]), float(a1[2])]
                a2 = [float(a2[0]), float(a2[1]), float(a2[2])]
                a3 = [float(a3[0]), float(a3[1]), float(a3[2])] 

    return 0
