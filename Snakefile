from datetime import date

today = date.today().strftime("%y_%m_%d")

const = 2.5
elem = "Fe"
sim_num = 1

output_name = "SIM_000" + str(sim_num) + "_WH_" + today


thickness_i = 3
thickness_f = 10
layers = range(thickness_i, thickness_f + 1)

T_i = 1
T_f = 100

temperatures = range(T_i, T_f)

CELLS = ["sc_100",
         "bcc_100",
         "bcc_110",
         "fcc_100",
         "fcc_111"]

rule all:
    input:
        expand("{cell}_{constant}_{layer}JAMS.cfg", cell=CELLS, constant=const, layer=layers)

# NOTE: The constant passed here is the nearest-neighbour distance
# between atomic sites - NOT the lattice constant/lattice parameter a.
# The conversion is done in `atomsk.py`
rule make_supercell:
    output:
        "{cell}_{constant}_{layer}JAMS.cfg"
    params:
        element=elem,
        cell=lambda wc: wc.cell,
        constant=lambda wc: wc.constant,
        layer=lambda wc: wc.layer
    script:
        "make_unitcells.py"

rule calculate_magnetisation:
    input:
        "{cell}_{constant}_{layer}JAMS.cfg"
    output:
        output_name + "/" + output_name + "_mag.tsv"
    params:
        folder=output[0].split('_')
    shell:
        "../jams {input} --output {output}"

