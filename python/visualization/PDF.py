import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as stat
from quantities_dict import quantities_dict
from extra_quantities import compute_extra_quantity
from cycler import cycler

from parameters import Parameters


plt.rc( ('axes'),labelsize=8.5)
plt.rc( ('xtick','ytick'),labelsize=8)
plt.rc("text",usetex=True)
plt.rc(("legend"),fontsize=8)

plt.rc( ('axes'),labelsize=11.5)
plt.rc( ('xtick','ytick'),labelsize=11)
plt.rc("text",usetex=True)
plt.rc(("legend"),fontsize=11)
plt.rcParams['axes.prop_cycle'] = cycler('color',['#1f78b4','#a6cee3','#33a02c','#b2df8a',                                       '#e31a1c','#fb9a99','#ff7f00','#fdbf6f',
                                      '#6a3d9a','#cab2d6'])

plt.rcParams['lines.linewidth'] = 1.75



def prepare_PDF(quantity, h5_output, h5_input=None, redshift=0,
                vmin=None, vmax=None, position=1.0, pdf_type='normal',
                pos_type='relative', ax=None,
                mass_bins=np.array([1e7,1e8,1e9,1e10,1e11,1e12]),
                select_only='all', plot_histogram=False, **args):
    """
    Plots the Probability Distribution Function (PDF) of a given quantity for
    different galaxy mass bins at a given redshift.
    Optionally, plots the PDF of the log of the quantity.

    Input: quantity -> string containing the th
           h5_output -> hdf5 file containing Magnetizer output data.
           h5_input -> hdf5 file containing Magnetizer input data. If absent,
                       will assume that h5_input = h5_output.
           mass_bins -> array containing the edges of the mass bins.
                        Default: [1e7,1e8,1e9,1e10,1e11,1e12]
           redshift -> redshift at which the PDF should be computed. Default: 0
           vmin, vmax -> maximum and minimum values for the quantity to be
                         considered. If None, the maximum and minium values in
                         the model are used. Default: None
           pdf_type -> 'log' computes P(log(quantity)) while
                       'normal' computes P(quantity). Default: 'normal'
           position -> 'position', in number of half-mass radii, for the
                       calculation of the PDF (for the radial dependent quantities). Default: 1.0
           ax -> matplotlib axis which will contain the plot.

    Returns: the matplotlib figure object.
    """

    if ax is None:
        ax = plt.subplot(1,1,1)

    if h5_input is None:
      h5_input = h5_output

    i_target = 0 
    
    params = Parameters(h5_output)
    # Selects mass-related dataset
    Mb = h5_input['Input']['Mstars_bulge']
    Md = h5_input['Input']['Mstars_disk']
    Mg = h5_input['Input']['Mgas_disk']
    M = Md  #+ Mg + Mb
    if 'label' in args:
        custom_label = True
    else:
        custom_label = False
    if 'central' in h5_input['Input']:
        central = h5_input['Input']['central']

    # Gets index of the selected redshift
    zs = h5_input['Input']['z'][:]
    iz = np.argmin(abs(zs - redshift))

    # Selects the dataset for the chosen quantity
    if quantity in h5_input['Input']:
        data = h5_input['Input'][quantity]
    else:
        ngals, ngrid, nzs = h5_output['Output']['n'].shape
        if quantity in h5_output['Output']:
            data = h5_output['Output'][quantity]
        else:
            data = None # I.e. to be computed later!
        if pos_type == 'relative':
            rmax_rdisc = params.grid['P_RMAX_OVER_RDISK']
            i_target = int(ngrid/rmax_rdisc*position)
        elif pos_type =='absolute':
            raise NotImplementedError
        else:
            raise ValueError

    for mmin, mmax in zip(mass_bins[:-1], mass_bins[1:]):
        # Creates filter for the current bin
        ok = M[:,iz]>mmin
        ok *= M[:,iz]<mmax

        if select_only=='satellite' or select_only=='satellites':
            ok *= ~ central[:,iz].astype(bool)
        elif select_only=='central' or select_only=='centrals':
            ok *= central[:,iz].astype(bool)

        # Skips if there aren't enough galaxies in the bin
        if not len(ok[ok])>1:
            continue

        # Loads the values at the specified mass bin, radius and redshift
        if data is not None:
            if len(data.shape) == 3:
                values = data[ok,i_target,iz]
            else:
                values = data[ok,iz]
        else:
            data_dict = {}
            for x in (h5_output['Output'],h5_input['Input']):
              for d in x:
                data_dict[d] = x[d]

            values = compute_extra_quantity(quantity,data_dict,
                                            ok,i_target,iz)
            
        # Removes invalid data
        filter_invalid = h5_output['Output']['n'][ok,i_target,iz] > 0
        filter_invalid *= np.isfinite(values)
        values = values[filter_invalid]

        # Sets maximum and minimum values
        if vmax is None:
            values_max = values.max()
        else:
            values_max = vmax
        if vmin is None:
            values_min = values.min()
        else:
            values_min = vmin

        if pdf_type == 'log':
            values = np.log10(values[values>-100])
        elif pdf_type !='normal':
            raise ValueError

        # Ignores mass bin if not enough galaxies with valid values
        values = values[np.isfinite(values)]
        if not len(values)>10:
            continue

        # Uses gaussian kernel density estimator to evaluate the PDF
        kernel = stat.gaussian_kde(values)

        if pdf_type == 'normal':
            x = np.linspace(values_min,values_max, 200)
        else:
            x = np.linspace(np.log10(values_min),np.log10(values_max), 200)
        y = kernel.evaluate(x)

        if pdf_type == 'log':  x = 10**x

        if not custom_label:
            args['label'] = r'$10^{{ {0:.2f} }}<M/M_\odot<10^{{ {1:.2f} }}$, $N={2}$'.format(np.log10(mmin),np.log10(mmax), len(values))

        ax.plot(x,y, **args)
        if plot_histogram:
            ax.hist(values, normed=True)
        ax.set_title(r'$z={0:.2f}$'.format(abs(zs[iz]))
                  )

    if quantity in quantities_dict:
        name, units = quantities_dict[quantity]
    else:
        name, units = quantity, None

    if pdf_type == 'normal':
        ax.set_ylabel(r'$P({0})$'.format(name, units))
    else:
        ax.set_xscale('log')
        if units:
            ax.set_ylabel(r'$P(\log({0}/{1}))$'.format(name, units))
        else:
            ax.set_ylabel(r'$P(\log({0}))$'.format(name))
    if units:
        ax.set_xlabel(r'${0}\,[{1}]$'.format(name, units))
    else:
        ax.set_xlabel(r'${0}$'.format(name))

    ax.legend(frameon=False)

    return ax

if __name__ == "__main__"  :

    h = h5py.File('/data/nlfsr/magnetizer-runs/GON9_5000_oldxi.hdf5','r')
    #h = h5py.File('/data/nlfsr/magnetizer_runs/GON9_5000.hdf5','r')
    #h = h5py.File('/data/nlfsr/magnetizer_runs/GON9_5000_noB_altxi.hdf5','r')

    hi = h5py.File('GON9_5000.hdf5','r')
    prepare_PDF('h/r', h5_output=h, h5_input=hi, pdf_type='normal',
                position=0.5, vmax=2)
    plt.show()