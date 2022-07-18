import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import numpy as np

def centered_meshgrid(x, y):
    '''
    Adds an extra row and column to X,Y and shifts all of the points so that
    when plotting with pcolormesh the data point is the center of the pixel.
    '''
    
    # we make new larger arrays and copy element by element because
    # np.resize or x_new.resize mangles the layout of the elements
    # because it resizes the linear memory without moving any elements
    # around
    x_new = np.zeros((x.shape[0]+1, x.shape[1]+1))
    y_new = np.zeros((y.shape[0]+1, y.shape[1]+1))
    
        
    # not strictly needed for this function, but pcolormesh will fail if they
    # are not the same shape
    assert(x_new.shape == y_new.shape)
    
    for i in range(0, x.shape[0]):
        for j in range(0, x.shape[1]):
            x_new[i, j] = x[i, j]
            
    for i in range(0, y.shape[0]):
        for j in range(0, y.shape[1]):
            y_new[i, j] = y[i, j]

    # we need to calculate the size of the pixels, i.e. the delta between points
    # NOTE: we assume points are evenly spaced
    
    # for maximum confusion, pcolormesh transposes things so the [0]
    # dimension shape is the y-axis and the [1] dimension is the x axis
    x_num_pixels = x.shape[1] - 1
    y_num_pixels = y.shape[0] - 1

    x_length = np.abs(x[0][-1] - x[0][1])
    y_length = np.abs(y[-1][0] - y[1][0])

    x_pixel_size = x_length / x_num_pixels
    y_pixel_size = y_length / y_num_pixels
    
    # populate the extra column
    x_new[:, -1] = x_new[:, -2] + x_pixel_size
    # copy the last row
    x_new[-1,:] = x_new[-2,:]
    
    # populate the extra row
    y_new[-1,:] = y_new[-2,:] + y_pixel_size
    # copy the last column
    y_new[:, -1] = y_new[:, -2]

    # shift all the data by half a pixel so the data points are 
    # now at the center of the pixel
    x_new = x_new - 0.5 * x_pixel_size
    y_new = y_new - 0.5 * y_pixel_size
    
    return x_new, y_new

def pandas_pcolor_plot(data, x_name, y_name, c_name, **args):
    '''
    Plot a pcolormesh color plot directly from a pandas dataFrame without
    using a pivot table hack.
    '''
    plot_data = data.copy()
    x_unique = data[x_name].unique()
    y_unique = data[y_name].unique()

    ncols=len(x_unique)
    nrows=len(y_unique)

    X = plot_data[x_name].values.reshape(nrows, ncols)
    Y = plot_data[y_name].values.reshape(nrows, ncols)
    C = plot_data[c_name].values.reshape(nrows, ncols)

    X,Y = centered_meshgrid(X, Y)

    plt.pcolormesh(X, Y, C, shading='flat', **args)

def plot_spectrum(path):
    '''
    Plot a magnon spectrum from a given path, encapsulating all the specific details
    of this simulation
    '''
    data = pd.read_table(path, sep="\s+")

    plt.ylabel('$E$ (meV)')

    if hasattr(snakemake.params, "ylim"):
        plt.ylim(0,snakemake.params["ylim"])
    else:
        plt.ylim(0,data['E_meV'].max())
    
    if hasattr(snakemake.params, "xlabel"):
        plt.xlabel(snakemake.params["xlabel"])
    else:
        plt.xlabel(r'$k$ ($2\pi / a$)')
        
    # set the aspect ratio (which is in the x/y axis units, not pixels or cm)
    if hasattr(snakemake.params, "aspect"):
        ax = plt.gca()
        ax.set_aspect(aspect=snakemake.params["aspect"])

    if hasattr(snakemake.params, "vmin"):
        colornorm = colors.LogNorm(vmin=snakemake.params["vmin"])
    else:
        colornorm = colors.LogNorm(vmin=1e-35)

    # make sure you use rasterized=True otherwise if you save the plot as a pdf you'll get nasty
    # results in some pdf viewers
    pandas_pcolor_plot(data, 'h', 'E_meV', 'Re_sqw_xx', cmap='jet', norm=colornorm, rasterized=True)
    
    # replace the xticks with the names of the symmetry points
    
plt.rcParams.update({'font.size': 14})

plot_spectrum(snakemake.input[0])

# plot phonon dispersion
v14 = 1.3     # nm/ps
a = 0.43337   # nm
hbar = 0.6582119569;  # meV ps
k = np.linspace(0,0.5)
E = v14*k*hbar*(2*np.pi / a)
plt.plot(k, E, "w--")
plt.plot(1.0-k, E, "w--")

plt.xlim(0,0.5)

ax = plt.gca()
#ax.set_aspect(aspect=0.3333)
#plt.xticks(np.arange(0.0, 1.0, 0.5))


if hasattr(snakemake.params, "title"):
    plt.title(snakemake.params["title"])
title = "Magnon spectrum at T = " + str(snakemake.params.get("temp"))

plt.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')
plt.close()

