import os
import subprocess
import numpy as np
import sys
import math

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
def generate_unit_cell(input_file, d_nbr, element, lattice_type):
    cells = generate_cells(d_nbr)
    
    # Extract information about the relevant cell
    orient = "orient " + cells[lattice_type]["orient"]
    a = truncate(cells[lattice_type]["a"], 3)
    d_plane = cells[lattice_type]["d_plane"]

    lattice = lattice_type.split('_')[0]

    args = ["atomsk",
            " --create ",
            lattice,
            str(a),
            element,
            orient,
            " -fractional ",
            input_file,
            "-nthreads",
            "1"]
    cmd = ' '.join(args)
    print("atomsk: Generating unit cell")
    print(cmd)
    subprocess.run(args)

    return d_plane

def generate_supercell(input_file, output_file, d_plane, layers):   
    args = ["atomsk ",
            input_file,
            " -duplicate 1 1 ",
            str(layers),
            " -fractional ",
            output_file]
    cmd = ' '.join(args)
    print("atomsk: Generating supercell.")
    print(cmd)

    #subprocess.run(args)

    return

# Argument format: Smallest number of monolayers, largest number of monolayers,
# lattice type/orientaion (specified in the format <lattice_type>_orientation)
def main():
    args = sys.argv[1:]

    if len(args) < 5:
        print("Expected 5 arguments, but only " + str(len(args)) + " were given.")
        sys.exit(0)
    thickness1 = int(args[0])
    thickness2 = int(args[1])
    
    lattice_type = args[2]
    d_nbr = float(args[3])
    element = args[4]
    
    input_file = lattice_type + "_" + element + ".cfg"

    d_plane = generate_unit_cell(input_file, d_nbr, element, lattice_type)
    
    for j in range(thickness1, thickness2 + 1):
        output_file = input_file.split('_')[0] + "_" + str(j) + ".cfg"
        generate_supercell(input_file, output_file, d_plane, j)

    return

main()
