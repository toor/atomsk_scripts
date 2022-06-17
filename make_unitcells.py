import subprocess
import numpy as np
import sys
import math
import re

def truncate(num, digits) -> float:
    dec = len(str(num).split('.')[1])
    if dec <= digits:
        return num
    step = 10.0 ** digits
    return math.trunc(step*num) / step

# This function takes in a parameter which specifies the distance between neighbouring
# atoms, and returns a dictionary of allowed unit cell types and relevant parameters
def generate_cells(a):
    cells = {
        "sc_100" : {
            "orient"  : "[100] [010] [001]",
            "a"       : a,
            "d_plane" : a
        },
        "bcc_100": {
            "orient"  : "[100] [010] [001]",
            "a"       : a,
            "d_plane" : a/2
        },
        "bcc_110": {
            "orient"  : "[1-10] [001] [110]",
            "a"       : a,
            "d_plane" : a/np.sqrt(2)
        },
        "fcc_100": {
            "orient"  : "[100] [010] [001]",
            "a"       : a,
            "d_plane" : a/2
        },
        "fcc_111": {
            "orient"  : "[11-2] [1-10] [111]",
            "a"       : a,
            "d_plane" : a /np.sqrt(3)
        },
    }

    return cells

# This generates the correct unit cell for a specified element,
# atomic spacing and lattice type/orientation.
def generate_unit_cell(input_file, a, lattice_type, element):
    cells = generate_cells(a)

    orient = cells[lattice_type]["orient"].split(' ')
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
            "-fractional",
            "-ow",
            input_file,
            "cfg"]
    # convert to XCrysDen XSF format for easy conversion
    # to JAMS
    convert_unit_cell_xsf = ["atomsk",
            input_file,
            "-ow",
            "-fractional",
            "xsf"]

    subprocess.run(make_unit_cell,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT)
    #subprocess.run(convert_unit_cell_xsf)

    return d_plane

def generate_supercell(input_file, output_file, d_plane, layers):   
    print("Generating a supercell with " + str(layers) + " layers from input " + input_file)
    print("Supercell name: " + output_file)
    cut_plane = str((layers - 1)*d_plane + 0.01)
    cell_z = str((layers - 1)*d_plane)

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
            "-fractional",
            "-ow",
            output_file,
            "cfg"]

    cut_layers = ["atomsk",
            output_file,
            "-cut",
            "above",
            cut_plane,
            "z",
            "-ow",
            output_file, 
            "cfg"]
    print("before supercell")
    subprocess.run(make_supercell,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT)
    print("ran supercell")
    #subprocess.run(cut_layers,
            #stdout=subprocess.DEVNULL,
            #stderr=subprocess.STDOUT)
    print("ran cut_layers")
    return

def read_unit_cell(filename, element, constant):   
    a1 = []
    a2 = []
    a3 = []
    lattice_vectors = [a1, a2, a3]
    coords = []

    # this function should read the correct parameters for a unit cell
    with open(filename) as f:
        lines = [line for line in f]

        header_start = 0
        atom_pos_start = 0

        for j in range(0, len(lines)):
            # check if line contains equals sign and whether
            # header_start has been updated yet. This if statement
            # should only be entered if the latter is zero.
            # Lattice vectors are guaranteed to always appear directly after
            # the definition of the Angstrom
            cond1 = re.search("A = ", lines[j]) 
            #cond2 = re.search("#", lines[j])
            if bool(cond1) and (not header_start):
                # Lattice vectors appear directly after the declaration of the
                # Angstrom unit.
                header_start = j + 1
            # same spiel here - the second condition is not strictly
            # necessary
            if (lines[j].strip() == element) and not atom_pos_start:
                atom_pos_start = j + 1
                print("Positions of atoms start at line: " + str(atom_pos_start))

        a = lines[header_start:(header_start + 3)]
        b = lines[(header_start + 3):(header_start + 6)]
        c = lines[(header_start + 6):(header_start + 9)]
       
        positions = lines[atom_pos_start:len(lines)]
        for pos in positions:
            coords.append(pos.strip('\t\n\r').split('    '))

        for i, x in enumerate([a, b, c]):
            lattice_vectors[i].append(str((float(x[0].split('=')[1].strip()))/constant))
            lattice_vectors[i].append(str((float(x[1].split('=')[1].strip()))/constant))
            lattice_vectors[i].append(str((float(x[2].split('=')[1].strip()))/constant)) 

        a1 = ', '.join(a1)
        a2 = ', '.join(a2)
        a3 = ', '.join(a3)

        return a1, a2, a3, coords           

def convert_jams(element, a, input_file, output_file):
    A = 1e-10 # 1 angstrom
    a1, a2, a3, coords = read_unit_cell(input_file, element, a)
   
    # See https://github.com/stonerlab/jams.git for more examples
    # of JAMS' config file format - it uses the libconfig library
    with open(output_file, "w") as f:
        f.write("unitcell = {\n")
        f.write("  parameter = " + str(float(a)*A) + ";\n")
        f.write("\n")
        f.write("  basis = (\n")
        f.write("    [" + a1 + "],\n")
        f.write("    [" + a2 + "],\n")
        f.write("    [" + a3 + "]);\n")
        f.write("\n")
        f.write("  positions = (\n")
        
        atom_count = len(coords)
        for j in range(0, (atom_count)):
            coord = ', '.join(coords[j])
            f.write("    (\"" + element + "\", [" + str(coord) + "])")
            if j==(atom_count - 1):
                f.write(");\n")
            else:
                f.write(",\n")
        f.write("};\n")
        f.close()

        return

def main(element, lattice_type, a, layers, output_file1, output_file2):   
    input_file = lattice_type + "_" + element + ".cfg"

    d_plane = generate_unit_cell(input_file, float(a), lattice_type, element)
    generate_supercell(input_file, output_file1, d_plane, int(layers))

    convert_jams(element, a, output_file1, output_file2)

    return

main(snakemake.params[0],
    snakemake.params[1],
    snakemake.params[2],
    snakemake.params[3],
    snakemake.output[0],
    snakemake.output[1])
