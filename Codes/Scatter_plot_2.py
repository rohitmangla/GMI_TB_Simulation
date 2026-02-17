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
spread = np.zeros((Sim_high_final.shape[0], channels_high), dtype=np.float64)
for i in range (Sim_high_final.shape[0]):
    for j in range (channels_high):
        Tb_val = [Sim_high_final[i, j, 4], Sim_high_final[i, j, 9], Sim_high_final[i, j, 5], Sim_high_final[i, j, 12],
                  Sim_high_final[i, j, 16], Sim_high_final[i, j, 24], Sim_high_final[i, j, 25]]
        spread[i,j] = np.nanmax(Tb_val)-np.nanmin(Tb_val)

# 2D histogram plots for single channel and single itable
cmap = plt.cm.jet
# Scatter plot 19 V GHz
fig = plt.figure(figsize=(12,12))
# Scatter plot 166 V
ax1= plt.subplot(2,2,1)
ax1.set_title('166 V')
cs1=ax1.hist2d(Sim_high_final[:, 9, 4], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs1 =ax1.hist2d(Sim_high_final[:, 9, 9], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs1 =ax1.hist2d(Sim_high_final[:, 9, 5], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs1 =ax1.hist2d(Sim_high_final[:, 9, 12], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs1 =ax1.hist2d(Sim_high_final[:, 9, 16], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs1 =ax1.hist2d(Sim_high_final[:, 9, 24], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs1 =ax1.hist2d(Sim_high_final[:, 9, 25], spread[:, 0], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
ax1.set_xlim(100, 300)
ax1.set_ylim(0, 15)
ax1.set_xlabel('Simulated Tb [K]')
ax1.set_ylabel('Spread')
cbar = fig.colorbar(cs1[3], ax=ax1, extend='both')

# Scatter plot 166 H
ax2= plt.subplot(2,2,2)
ax2.set_title('166 H')
cs2 =ax2.hist2d(Sim_high_final[:, 10, 4], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs2 =ax2.hist2d(Sim_high_final[:, 10, 9], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs2 =ax2.hist2d(Sim_high_final[:, 10, 5], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs2 =ax2.hist2d(Sim_high_final[:, 10, 12], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs2 =ax2.hist2d(Sim_high_final[:, 10, 16], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs2 =ax2.hist2d(Sim_high_final[:, 10, 24], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs2 =ax2.hist2d(Sim_high_final[:, 10, 25], spread[:, 1], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
ax2.set_xlim(100, 300)
ax2.set_ylim(0, 15)
ax2.set_xlabel('Simulated Tb [K]')
ax2.set_ylabel('Spread')
cbar = fig.colorbar(cs2[3], ax=ax2, extend='both')

# Scatter plot 183+3 V
ax3= plt.subplot(2,2,3)
ax3.set_title('$183\pm\ 3V$')
cs3 =ax3.hist2d(Sim_high_final[:, 11, 4], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs3 =ax3.hist2d(Sim_high_final[:, 11, 9], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs3 =ax3.hist2d(Sim_high_final[:, 11, 5], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs3 =ax3.hist2d(Sim_high_final[:, 11, 12], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs3 =ax3.hist2d(Sim_high_final[:, 11, 16], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs3 =ax3.hist2d(Sim_high_final[:, 11, 24], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs3 =ax3.hist2d(Sim_high_final[:, 11, 25], spread[:, 2], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
ax3.set_xlim(100, 300)
ax3.set_ylim(0, 15)
ax3.set_xlabel('Simulated Tb [K]')
ax3.set_ylabel('Spread')
cbar = fig.colorbar(cs3[3], ax=ax3, extend='both')

# Scatter plot 183+7 V
ax4= plt.subplot(2,2,4)
ax4.set_title('$183\pm\ 7V$')
cs4 =ax4.hist2d(Sim_high_final[:, 12, 4], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs4 =ax4.hist2d(Sim_high_final[:, 12, 9], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs4 =ax4.hist2d(Sim_high_final[:, 12, 5], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs4 =ax4.hist2d(Sim_high_final[:, 12, 12], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs4 =ax4.hist2d(Sim_high_final[:, 12, 16], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs4 =ax4.hist2d(Sim_high_final[:, 12, 24], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
cs4 =ax4.hist2d(Sim_high_final[:, 12, 25], spread[:, 3], bins=[80,15], range=[[100, 300],[0,15]], cmin=5, cmap= cmap, norm=mpl.colors.LogNorm())
ax4.set_xlim(100, 300)
ax4.set_ylim(0, 15)
ax4.set_xlabel('Simulated Tb [K]')
ax4.set_ylabel('Spread')
cbar = fig.colorbar(cs4[3], ax=ax4, extend='both')
plt.savefig('D:/Python_processing/cyclone/results/scatter_plot.png', dpi=300, bbox_inches= 'tight')
plt.show()

error
# Scatter plot 183+7 V
ax5= plt.subplot(2,3,6)
ax5.set_title('$183\pm\ 7V$')
cs5=ax5.hist2d(Sim_high_final[:, 12, 4], FG_high[:, 3], bins=50, range=[[150, 300],[-50,50]], cmin=10, cmap= cmap, norm=mpl.colors.LogNorm())
ax5.set_xlim(150, 300)
ax5.set_ylim(-100, 100)
divider = make_axes_locatable(ax5)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cs5[3], cax=cax, orientation='vertical', extend='both')
ax5.set_xlabel('Simulated Tb [K]')
ax5.set_ylabel('Observed Tb [K]')
plt.show()
plt.savefig('D:/Python_processing/cyclone/results/scatter_plot.png', dpi=300, bbox_inches= 'tight')
