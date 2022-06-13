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

    orient = cells[lattice_type]["orient"].split(' ')
    print(orient)
    a = truncate(cells[lattice_type]["a"], 3)
    d_plane = cells[lattice_type]["d_plane"]

    lattice = lattice_type.split('_')[0]
    
    make_unit_cell = ["atomsk",
            "--create",
            lattice,
            str(a),
            element,
            "orient",
            orient[0],
            orient[1],
            orient[2],
            "-ow",
            "-fractional",
            input_file,
            "cfg"]
    # convert to XCrysDen XSF format for easy conversion
    # to JAMS
    convert_unit_cell_xsf = ["atomsk",
            input_file,
            "-ow",
            "xsf"]

    subprocess.run(make_unit_cell)
    subprocess.run(convert_unit_cell_xsf)

    return d_plane

def generate_supercell(input_file, output_file, d_plane, layers):   
    cut_plane = str((layers - 1)*d_plane + 0.01)
    cell_z = str(layers*d_plane)

    make_supercell = ["atomsk",
            input_file,
            "-duplicate",
            "1",
            "1",
            str(layers),
            "-cut",
            "above",
            cut_plane,
            "z",
            "-cell",
            "set",
            cell_z,
            "H3",
            "-fractional",
            "-ow",
            output_file,
            "cfg"]

    convert_supercell_xsf = ["atomsk",
            output_file,
            "-ow",
            "xsf"]

    subprocess.run(make_supercell)
    # convert to an XSF file for visualisation and easy conversion
    # to the relevant file format
    subprocess.run(convert_supercell_xsf)
    return

# Argument format: Smallest number of monolayers, largest number of monolayers,
# lattice type/orientaion (specified in the format <lattice_type>_orientation)
def main(lattice_type, a, layers, output_file):   
    input_file = lattice_type + "_" + element + ".cfg"

    d_plane = generate_unit_cell(input_file, float(a), lattice_type)
    generate_supercell(input_file, output_file, d_plane, int(layers))
    
    return

main(snakemake.params[0], snakemake.params[1], snakemake.params[2], snakemake.output[0])
