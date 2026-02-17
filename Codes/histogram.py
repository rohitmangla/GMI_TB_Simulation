import datetime
import numpy as np
from netCDF4 import Dataset
import os
os.environ['PROJ_LIB'] = r'C:\Users\rohit\Anaconda3\pkgs\proj4-5.2.0-ha925a31_1\Library\share'
import matplotlib
import pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker, cm

matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['legend.fontsize'] = 14.
matplotlib.rcParams['font.size'] = 26.
matplotlib.rcParams['axes.labelsize'] = 24.
matplotlib.rcParams['xtick.labelsize'] = 24.
matplotlib.rcParams['ytick.labelsize'] = 24.
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.hspace'] = 0.2
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = "D:/Python_processing/cyclone/results/histogram/"
mie_tables = 26
channels_low =9
channels_high =4
nbins = 40
ncfile = Dataset(PATH+'count.nc', 'r')
print(ncfile)
obs_high = (ncfile.variables['obs_high'])[:]
sim_high= (ncfile.variables['sim_high'])[:][:]
obs_high_tp = np.transpose(obs_high)
sim_high_tp = np.transpose(sim_high)

cmap = matplotlib.cm.RdYlBu
z_min = 100
z_max = 300
res_colorbar = 5.
clevs = np.arange(0, (z_max - z_min) / res_colorbar + 1) * res_colorbar + z_min
xedges = np.linspace(100, 300, num=41)
yedges = np.linspace(0, 4, num=5)

for i in range (mie_tables):
    fig = plt.figure(figsize=(30, 10))
    ax1 = plt.subplot(1, 2, 1)
    im_obs_norm = ax1.imshow((obs_high_tp /np.nanmax(obs_high_tp)), interpolation='nearest', origin='low', \
                             extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto')
    ax1.set_title('(a) Observed Histogram')
    ax1.set_xlabel('Brightness Temp (K)')
    ax1.set_ylabel('Channel name')
    ax1.set_ylim(0, 4)
    ax1.set_yticks(np.arange(0, 4, step=1))
    ax1.set_xlim(100, 300)
    ax1.set_yticklabels(labels= ['166V','166H','$183\pm\ 3V$', '$183\pm\ 7V$'])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(im_obs_norm, cax=cax, extend='both')
    vmin = cbar.vmin
    vmax = cbar.vmax
    im_obs_norm.set_clim(vmin, vmax)
    cbar.set_label('counts in bin')
    plt.tight_layout()

    ax2 = plt.subplot(1, 2, 2)
    im_sim_norm = ax2.imshow((sim_high_tp[i, :, :] /np.nanmax(sim_high_tp[i, :, :])), interpolation='nearest', origin='low', \
                             extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto')
    ax2.set_title('(b) Simulated Histogram')
    ax2.set_xlabel('Brightness Temp (K)')
    ax2.set_ylabel('Channel name')
    ax2.set_ylim(0, 4)
    ax2.set_xlim(100, 300)
    ax2.set_yticks(np.arange(0.5, 4, step=1))
    ax2.set_yticklabels(labels= ['166V','166H','$183\pm\ 3V$', '$183\pm\ 7V$'])
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar = plt.colorbar(im_sim_norm, cax=cax, extend='both')
    vmin = cbar.vmin
    vmax = cbar.vmax
    im_sim_norm.set_clim(vmin, vmax)
    cbar.set_label('counts in bin')
    plt.tight_layout()
    plt.show()
    plt.savefig('D:/Python_processing/cyclone/results/histogram/hist_'+str(i)+'.png', dpi=300,bbox_inches='tight')