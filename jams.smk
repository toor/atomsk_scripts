import os
import numpy as np

localrules: gen_unitcell, render_equi_cfg, render_mag_cfg, render_magnon_cfg, analyse_magnetisation, plot_magnetisation, plot_magnon_spectrum

constant=2.5

temperatures = [f'{x:.1f}' for x in range(10, 710, 10)]
 
d_nbrs = {
    "sc_100": constant,
    "bcc_100": (constant*np.sqrt(3))/2,
    "bcc_110": (constant*np.sqrt(3))/2,
    "fcc_100": constant/np.sqrt(2),
    "fcc_111": (constant*np.sqrt(3))/np.sqrt(2)
}

cells = {
    #"sc_100"
    "bcc_100",
    "bcc_110",
    "fcc_100",
    "fcc_111"
}

# applied field in Tesla
app_field = 0.1
# Number of times to repeat the unit cell in-plane
repeats=128

# generate the supercell
rule gen_unitcell:
    output:
        "results/{lattice}/{layer}/supercell.cfg",
        "results/{lattice}/{layer}/unitcell.cfg"
    params:
        element="Fe",
        lattice=lambda wc: wc.lattice,
        const=constant,
        layer=lambda wc: wc.layer,
    script:
        "make_unitcells.py" 

# produce the desired JAMS config from a jinja file
rule render_equi_cfg:
    input:
        "JAMS_equilibration.jinja2.cfg"
    output:
        "results/{lattice}/{layer}/{T2}K/jams_equilibration.cfg"
    params:
        r_cutoff=constant,
        field=app_field,
        rep=repeats,
        gilbert_damp=0.1,
        mu_s=3
    template_engine:
        "jinja2"

rule render_magnon_cfg:
    input:
        "JAMS_magnon.jinja2.cfg"
    output:
        "results/{lattice}/{layer}/{T2}K/jams_magnon.cfg"
    params:
        r_cutoff=constant,
        field=app_field,
        rep=repeats,
        gilbert_damp=0.01,
        mu_s=3
    template_engine:
        "jinja2"

# Render a different config that is used to calculate the magnon spectrum
rule render_mag_cfg:
    input:
        "JAMS_mag.jinja2.cfg"
    output:
        "results/{lattice}/{layer}/{T}K/jams_mag.cfg"
    params:
        r_cutoff=lambda wc: d_nbrs[wc.lattice],
        field=app_field,
        rep=repeats,
        gilbert_damp=0.1,
        mu_s=3
    template_engine:
        "jinja2"

# We need to allow the system to equilibrate first before calculating the magnon spectrum.
# This is done with a higher value of the gilbert damping parameter as specified
# in render_magnon_cfg above. After the system has been thermalised, the states of the spins
# in this thermalised state are passed to JAMS to determine the magnon spectrum. We use a more
# realistic value of the Gilbert damping parameter to determine the magnon spectrum, which also
# makes the spectrum less diffuse (WHY?)
rule equilibration:
    input:
        "results/{lattice}/{layer}/{T2}K/jams_equilibration.cfg",
        "results/{lattice}/{layer}/unitcell.cfg"
    output:
        "results/{lattice}/{layer}/{T2}K/jams_equilibrated_final.h5"
    shell:
        "../jams_latest --output=\"results/{wildcards.lattice}/{wildcards.layer}/{wildcards.T2}K\" --name=\"jams_equilibrated\" {input} > results/{wildcards.lattice}/{wildcards.layer}/{wildcards.T2}K/jams_equilibration.log"

rule calculate_thermodynamics:
    input:
        "results/{lattice}/{layer}/{T}K/jams_mag.cfg",
        "results/{lattice}/{layer}/unitcell.cfg"
        #"{lattice}/{layer}/{T}K/jams_equilibrated_final.h5"
    output:
        "results/{lattice}/{layer}/{T}K/jams_mag.tsv"
    shell:
        "../jams_latest --output=\"results/{wildcards.lattice}/{wildcards.layer}/{wildcards.T}K\" --name=\"jams\" {input[0]} {input[1]}  > results/{wildcards.lattice}/{wildcards.layer}/{wildcards.T}K/jams_magnetisation.log" 

# \'lattice: {{spins=\"{input[2]}\";}};\'

rule calculate_magnon_spectrum:
    input:
        "results/{lattice}/{layer}/{T2}K/jams_magnon.cfg",
        "results/{lattice}/{layer}/unitcell.cfg",
        "results/{lattice}/{layer}/{T2}K/jams_equilibrated_final.h5"
    output:
        "results/{lattice}/{layer}/{T2}K/jams_magnon_magnon_spectrum_path_0.tsv"
    shell:
        "../jams_latest --output=\"results/{wildcards.lattice}/{wildcards.layer}/{wildcards.T2}K\" --name=\"jams_magnon\" {input[0]} {input[1]} \'lattice: {{spins=\"{input[2]}\";}};\' > results/{wildcards.lattice}/{wildcards.layer}/{wildcards.T2}K/jams_magnon.log" 

rule analyse_magnetisation:
    input:
        expand("results/{{lattice}}/{{layer}}/{T}K/jams_mag.tsv", T=temperatures)
    output:
        "results/{lattice}/{layer}/mag_vs_temp.dat",
        "results/{lattice}/{layer}/sus_vs_temp.dat"
    params:
        lattice=lambda wc: wc.lattice,
        layer=lambda wc: wc.layer,
        const=constant,
        field=app_field,
        rep=repeats,
        mu_s=3
    script:
        "process_mag_data.py"

rule plot_magnetisation:
    input:
        "results/{lattice}/{layer}/mag_vs_temp.dat",
        "results/{lattice}/{layer}/sus_vs_temp.dat"
    output:
        "results/{lattice}/{layer}/mag_vs_temp.png",
        "results/{lattice}/{layer}/sus_vs_temp.png"
    params:
         lattice=lambda wc: wc.lattice,
         layer=lambda wc: wc.layer,
         const=constant,
         field=app_field,
         mu_s=3,
         rep=repeats
    script:
        "plot_mag_data.py"

rule plot_magnon_spectrum:
    input:
        "results/{lattice}/{layer}/{T2}K/jams_magnon_magnon_spectrum_path_0.tsv"
    output:
        "figures/{lattice}/{layer}/magnon_{T2}K/spectrum.pdf"
    params:
        #ylim=20.0,
        vmin=1e-12,
        temp=lambda wc: wc.T2
    script:
        "plot_magnon_spectrum.py"
