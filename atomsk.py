import os
import subprocess
import numpy as np
import sys

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
def generate_unit_cell(d_nbr, element, lattice_type):
    cells = generate_cells(d_nbr)
    
    # Extract information about the relevant cell
    orient = "orient " + cells[lattice_type]["orient"]
    a = cells[lattice_type]["a"]
    d_plane = cells[lattice_type]["d_plane"]

    lattice = lattice_type.split('_')[0]

    filename = lattice_type + "_" + element + ".xsf"
    rm = "rm " + filename
    os.system(rm)

    args = ["atomsk --create",
            lattice,
            str(a),
            element,
            orient,
            filename]
    cmd = ' '.join(args)
    print(cmd)

    os.system(cmd)

    return filename, d_plane

def generate_supercell(filename, d_plane, layer):
    cut_plane =  (float(layer) - 1)*d_plane + 0.01
    cell_z = float(layer)*d_plane
    
    filename2 = filename.split('.')[0] + "_" + str(layer) + ".xsf"
    # remove annoying duplicates
    rm = "rm " + filename2
    os.system(rm)

    args = ["atomsk " + filename,
            "-duplicate 1 1",
            str(layer),
            #"-cut above ",
            #str(cut_plane),
            #"z",
            #"-cell set ",
            #str(d_plane),
            #"H3",
            #"-fractional",
            filename2]
    cmd = ' '.join(args)
    print(cmd)

    os.system(cmd)

    return

# Argument format: Smallest number of monolayers, largest number of monolayers,
# lattice type/orientaion (specified in the format <lattice_type>_orientation)
def main():
    if len(sys.argv[1:]) < 5:
        print("Expected 5 arguments, but only " + str(len(sys.argv) - 1) + " were given.")
        sys.exit(0)
    args = sys.argv[1:] 
    
    thickness1 = int(args[0])
    thickness2 = int(args[1])
    
    lattice_type = args[2]
    d_nbr = float(args[3])
    element = args[4]
    
    filename, d_plane = generate_unit_cell(d_nbr, element, lattice_type)
    
    for j in range(thickness1, thickness2 + 1):
        generate_supercell(filename, d_plane, j)

    return

main()
