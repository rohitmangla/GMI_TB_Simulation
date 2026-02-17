import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.rcParams['axes.titlesize'] = 18
matplotlib.rcParams['legend.fontsize'] =12
matplotlib.rcParams['font.size']= 18
matplotlib.rcParams['axes.labelsize']= 17
matplotlib.rcParams['xtick.labelsize'] = 17
matplotlib.rcParams['ytick.labelsize'] = 17
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_low_final = np.zeros((1, channels_low), dtype=np.float64)
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_low_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)

for ifile in FILES:
    print(ifile)
    ncfile = Dataset(ifile,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_low  = ncfile.variables['Simulated_Tb_low'][:]   # low frequency Simulation
    Simulated_Tb_high = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    GMI_Tb_low        = ncfile.variables['GMI_gridded_low'][:]   # low frequency GMI Tb
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
    lsm               = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
    # Land masking
    ind_land = np.where(lsm ==1)
    Simulated_Tb_low[ind_land, :, :]= np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]      = np.nan
    GMI_Tb_high[ind_land, :]   = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)

GMI_low_final[0, :] = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
Sim_high_final[0, :, :] = np.nan

FG_low = GMI_low_final-Sim_low_final[:, 0:9, 4]
FG_high = GMI_high_final-Sim_high_final[:, 9:13, 4]

# occurrence percentage
count_19V  = np.zeros((100, 100), dtype = np.float64)
initial_sim = 100.
ini_err  = -50.
for i in range (100): # simulated radiances
    ind = np.where((Sim_low_final[:, 2, 4] >= initial_sim) * (Sim_low_final[:, 2, 4] < initial_sim+2))
    err = np.squeeze(FG_low[ind, 2])
    for j in range (100):  # error
        ind_err = np.where((err>=ini_err)*(err<ini_err+1))
        count_19V[i,j] = np.size(ind_err)
        ini_err = ini_err+1
    initial_sim = initial_sim+2
    ini_err = -50.
error



# 2D histogram plots for single channel and single itable
cmap = plt.cm.RdBu
# Scatter plot 19 V GHz
fig = plt.figure(figsize=(25, 12))
ax0= plt.subplot(2,3,1)
ax0.set_title('19V')
cs0 =ax0.hist2d(Sim_low_final[:, 2, 4], FG_low[:, 2], bins=50, range=[[150, 300],[150,300]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax0.set_xlim(150, 300)
ax0.set_ylim(150, 300)
divider = make_axes_locatable(ax0)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs0[3], cax=cax, orientation='vertical', extend='both')
ax0.set_xlabel('Simulated Tb [K]')
ax0.set_ylabel('Observed Tb [K]')

error
# scatter plot 23V
ax1= plt.subplot(2,3,2)
ax1.set_title('23V')
cs1 =ax1.hist2d(Sim_low_final[:, 4, 4], GMI_low_final[:, 4], bins=50, range=[[150, 300],[150,300]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax1.set_xlim(150, 300)
ax1.set_ylim(150, 300)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs1[3], cax=cax, orientation='vertical', extend='both')
ax1.set_xlabel('Simulated Tb [K]')
ax1.set_ylabel('Observed Tb [K]')

# scatter plot 37V
ax2= plt.subplot(2,3,3)
ax2.set_title('37V')
cs2 =ax2.hist2d(Sim_low_final[:, 5, 4], GMI_low_final[:, 5], bins=50, range=[[150, 300],[150,300]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax2.set_xlim(150, 300)
ax2.set_ylim(150, 300)
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs2[3], cax=cax, orientation='vertical', extend='both')
ax2.set_xlabel('Simulated Tb [K]')
ax2.set_ylabel('Observed Tb [K]')


# Scatter plot 89 V
ax3= plt.subplot(2,3,4)
ax3.set_title('89V')
cs3 =ax3.hist2d(Sim_low_final[:, 7, 4], GMI_low_final[:, 7], bins=50, range=[[150, 300],[150,300]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax3.set_xlim(150, 300)
ax3.set_ylim(150, 300)
divider = make_axes_locatable(ax3)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs3[3], cax=cax, orientation='vertical', extend='both')
ax3.set_xlabel('Simulated Tb [K]')
ax3.set_ylabel('Observed Tb [K]')

# Scatter plot 166 V
ax4= plt.subplot(2,3,5)
ax4.set_title('166 V')
cs4=ax4.hist2d(Sim_high_final[:, 9, 4], GMI_high_final[:, 0], bins=50, range=[[150, 300],[150,300]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax4.set_xlim(150, 300)
ax4.set_ylim(150, 300)
divider = make_axes_locatable(ax4)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs4[3], cax=cax, orientation='vertical', extend='both')
ax4.set_xlabel('Simulated Tb [K]')
ax4.set_ylabel('Observed Tb [K]')

# Scatter plot 183+7 V
ax5= plt.subplot(2,3,6)
ax5.set_title('$183\pm\ 7V$')
cs5=ax5.hist2d(Sim_high_final[:, 12, 4], GMI_high_final[:, 3], bins=50, range=[[150, 300],[150,300]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax5.set_xlim(150, 300)
ax5.set_ylim(150, 300)
divider = make_axes_locatable(ax5)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs5[3], cax=cax, orientation='vertical', extend='both')
ax5.set_xlabel('Simulated Tb [K]')
ax5.set_ylabel('Observed Tb [K]')
plt.show()
plt.savefig('D:/Python_processing/cyclone/results/scatter_plot.png', dpi=300, bbox_inches= 'tight')
