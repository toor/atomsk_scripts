const = 2.5
elem = "Fe"

layers = []
thickness_i = 3
thickness_f = 10

for i in range(thickness_i, thickness_f + 1):
    layers.append(i)

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
        "atomsk.py"
