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
        expand("{cell}_{constant}_{layer}.cfg", cell=CELLS, constant=2.5, layer=layers)

# NOTE: The constant passed here is the nearest-neighbour distance
# between atomic sites - NOT the lattice constant/lattice parameter a.
# The conversion is done in `atomsk.py`
rule make_supercell:
    output:
        "{cell}_{constant}_{layer}.cfg"
    params:
        cell=lambda wc: wc.cell,
        constant=lambda wc: wc.constant,
        layer=lambda wc: wc.layer
    script:
        "atomsk.py"

#rule cleanup:
    #shell:
        #"rm *.cfg"

#rule make_jams_cell:
#    input:
#        "{cell}_{constant}_{layer}.cfg"
#    output:
#        "{cell}_{constant}_{layer}_JAMS.cfg"
#    script:
#        "make_jams_cfg.py"
