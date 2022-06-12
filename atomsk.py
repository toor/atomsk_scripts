import os
import subprocess
import numpy as np
import sys
import math

element = "Fe"

def truncate(num, digits) -> float:
    dec = len(str(num).split('.')[1])
    if dec <= digits:
        return num
    step = 10.0 ** digits
    return math.trunc(step*num) / step

# This function takes in a parameter which specifies the distance between neighbouring
# atoms, and returns a dictionary of allowed unit cell types and relevant parameters
def generate_cells(d_nbr):
    cells = {
        "sc_100" : {
            "orient"  : "[100] [010] [001]",
            "a"       : d_nbr,
            "d_plane" : d_nbr
        },
        "bcc_100": {
            "orient"  : "[100] [010] [001]",
            "a"       : d_nbr * (2/np.sqrt(3)),
            "d_plane" : d_nbr/np.sqrt(3)
        },
        "bcc_110": {
            "orient"  : "[1-10] [001] [110]",
            "a"       : d_nbr * (2/np.sqrt(3)),
            "d_plane" : d_nbr * (np.sqrt(2)/np.sqrt(3))
        },
        "fcc_100": {
            "orient"  : "[100] [010] [001]",
            "a"       : d_nbr * np.sqrt(2),
            "d_plane" : d_nbr * (np.sqrt(2)/2)
        },
        "fcc_111": {
            "orient"  : "[11-2] [1-10] [111]",
            "a"       : d_nbr * np.sqrt(2),
            "d_plane" : d_nbr * (np.sqrt(2)/np.sqrt(3))
        },
    }

    return cells

# This generates the correct unit cell for a specified element,
# atomic spacing and lattice type/orientation.
# Returns: the filename of the relevant XSF file produced as well as the
def generate_unit_cell(input_file, d_nbr, lattice_type):
    cells = generate_cells(d_nbr)
    
    # Extract information about the relevant cell
    orient = "orient " + cells[lattice_type]["orient"]
    a = truncate(cells[lattice_type]["a"], 3)
    d_plane = cells[lattice_type]["d_plane"]

    lattice = lattice_type.split('_')[0]

    args = ["atomsk",
            "--create",
            lattice,
            str(a),
            element,
            orient,
            "-fractional",
            input_file,
            "cfg",
            "-ow"]
    cmd = ' '.join(args)
    print("atomsk: Generating unit cell")
    print(cmd)
    # TODO: subprocess.run() seems to perhaps pass the wrong intput to atomsk, in that
    # atomsk reads the orientation argument as the filename and thus names the files
    # wrongly, causing atomsk to pick up filenames incorrectly. Also, I wonder whether
    # snakemake is capturing the unit cell files, in the format:
    # <lattice>_<orientation>_<element>.cfg.
    os.system(cmd)
    #subprocess.run(args)

    return d_plane

def generate_supercell(input_file, output_file, d_plane, layers):   
    args = ["atomsk",
            input_file,
            "-duplicate",
            "1 1 " + str(layers),
            "-fractional",
            output_file,
            "cfg",
            "-ow"]
    cmd = ' '.join(args)
    print("atomsk: Generating supercell.")
    print(cmd)

    #subprocess.run(args)
    os.system(cmd)

    return

# Argument format: Smallest number of monolayers, largest number of monolayers,
# lattice type/orientaion (specified in the format <lattice_type>_orientation)
def main(lattice_type, a, layers, output_file):   
    input_file = lattice_type + "_" + element + ".cfg"

    d_plane = generate_unit_cell(input_file, float(a), lattice_type)
    generate_supercell(input_file, output_file, d_plane, int(layers))
    
    return

main(snakemake.params[0], snakemake.params[1], snakemake.params[2], snakemake.output[0])
