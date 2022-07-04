import os
import numpy as np

localrules: gen_unitcell, render_cfg, analyse_magnetisation, plot_magnetisation

constant=2.5

d_nbrs = {
    "sc_100": constant,
    "bcc_100": constant*np.sqrt(3)/2,
    "bcc_110": constant*np.sqrt(3)/2,
    "fcc_100": constant/np.sqrt(2),
    "fcc_111": constant*np.sqrt(3)/np.sqrt(2)
}

temperatures = [f'{x:.1f}' for x in range(10, 410, 10)]

cells = {
    "sc_100",
    "bcc_100",
    "bcc_110",
    "fcc_100",
    "fcc_111"
}

layers = range(2, 6)

# applied field in Tesla
app_field = 0.1
# Number of times to repeat the unit cell in-plane
repeats=32


# generate the supercell
rule gen_unitcell:
    output:
        "{lattice}/{layer}/supercell.xsf",
        "{lattice}/{layer}/unitcell.cfg"
    params:
        element="Fe",
        lattice=lambda wc: wc.lattice,
        const=constant,
        layer=lambda wc: wc.layer,
    script:
        "make_unitcells.py" 

# produce the desired JAMS config from a jinja file
rule render_cfg:
    input:
        "JAMS_defaults.jinja2.cfg"
    output:
        "{lattice}/{layer}/{T}K/jams.cfg"
    params:
        d_nbr=lambda wc: d_nbrs[wc.lattice],
        field=app_field,
        rep=repeats
    template_engine:
        "jinja2"

rule calc_magnetisation:
    input:
        "{lattice}/{layer}/{T}K/jams.cfg",
        "{lattice}/{layer}/unitcell.cfg"
    output:
        "{lattice}/{layer}/{T}K/jams_mag.tsv"
    shell:
        "../jams --output=\"{wildcards.lattice}/{wildcards.layer}/{wildcards.T}K\" --name=\"jams\" {input}"

rule analyse_magnetisation:
    input:
        expand("{{lattice}}/{{layer}}/{T}K/jams_mag.tsv", T=temperatures)
    output:
        "{lattice}/{layer}/mag_vs_temp.dat"
    params:
        lattice=lambda wc: wc.lattice,
        layer=lambda wc: wc.layer,
        const=constant,
        field=app_field,
        rep=repeats
    script:
        "process_mag_data.py"

rule plot_magnetisation:
    input:
        "{lattice}/{layer}/mag_vs_temp.dat"
    output:
        "{lattice}/{layer}/mag_vs_temp.png"
    script:
        "plot_mag_data.py"
