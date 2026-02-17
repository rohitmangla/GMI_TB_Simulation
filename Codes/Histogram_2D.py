import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib
from matplotlib.colors import LogNorm
import ez_color

matplotlib.rcParams['axes.titlesize'] = 17
matplotlib.rcParams['legend.fontsize'] =12
matplotlib.rcParams['font.size']= 30
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.25

PATH = 'D:/Python_processing/cyclone/forecast_06hrs/'
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
    GMI_Tb_high[ind_land, :]     = np.nan
    ind_low = np.where(GMI_Tb_low[:, 0]>1)
    lat_low = lat[ind_low]
    lon_low = lon[ind_low]
    ind_high = np.where(GMI_Tb_high[:, 0] > 1)
    lat_high = lat[ind_high]
    lon_high = lon[ind_high]

    GMI_low   = np.squeeze(GMI_Tb_low[ind_low, :])              # (profiles, channels_low)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind_high, :])             # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind_low, :, :])     # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind_high, :, :])    # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)

GMI_low_final[0, :] = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
Sim_high_final[0, :, :] = np.nan
n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan

nbins = 80
bin_start = 100
dbin = 2.5
Tb_low_bin = np.zeros((nbins, channels_low, mie_tables), dtype=np.float64)
Tb_high_bin = np.zeros((nbins, channels_high, mie_tables), dtype=np.float64)
GMI_low_bin = np.zeros((nbins, channels_low), dtype=np.float64)
GMI_high_bin = np.zeros((nbins, channels_high), dtype=np.float64)
# Observed
for ichanel in range(channels_low):
    for ibin in range(nbins):
        ind_GMI_low_bin = np.where((GMI_low_final[:,ichanel]>=bin_start)*(GMI_low_final[:, ichanel]<bin_start+dbin))
        val = np.squeeze(GMI_low_final[ind_GMI_low_bin, ichanel])
        GMI_low_bin[ibin, ichanel] = np.size(val)
        bin_start = bin_start+dbin
    bin_start = 100
GMI_low_bin[GMI_low_bin==0]= np.nan

bin_start = 100
for ichanel in range(channels_high):
    for ibin in range(nbins):
        ind_GMI_high_bin = np.where((GMI_high_final[:,ichanel]>=bin_start)*(GMI_high_final[:, ichanel]<bin_start+dbin))
        val = np.squeeze(GMI_high_final[ind_GMI_high_bin, ichanel])
        GMI_high_bin[ibin, ichanel] = np.size(val)
        bin_start = bin_start+dbin
    bin_start = 100
GMI_high_bin[GMI_high_bin==0]= np.nan

bin_start = 100
# Simulations
for itable in range(mie_tables):
    for ichanel in range(channels_low):
        for ibin in range(nbins):
            ind_low_bin = np.where((Sim_low_final[:,ichanel,itable]>=bin_start)*(Sim_low_final[:, ichanel, itable]<bin_start+dbin))
            val = np.squeeze(Sim_low_final[ind_low_bin, ichanel, itable])
            Tb_low_bin[ibin, ichanel, itable] = np.size(val)
            bin_start = bin_start+dbin
        bin_start = 100
Tb_low_bin[Tb_low_bin==0]= np.nan
bin_start = 100
for itable in range(mie_tables):
    for ichanel in range(channels_high):
        for ibin in range(nbins):
            ind_high_bin = np.where((Sim_high_final[:,ichanel+9,itable]>=bin_start)*(Sim_high_final[:, ichanel+9, itable]<bin_start+dbin))
            val = np.squeeze(Sim_high_final[ind_high_bin, ichanel, itable])
            Tb_high_bin[ibin, ichanel, itable] = np.size(val)
            bin_start = bin_start+dbin
        bin_start = 100
Tb_high_bin[Tb_high_bin==0]= np.nan

# Tb_filter_Obs
Tb_filter_Obs = np.zeros((nbins, 7),dtype= np.float64)
Tb_filter_Obs[:, 0] = GMI_low_bin[:, 2]   # 19V
Tb_filter_Obs[:, 1] = GMI_low_bin[:, 4]   # 23V
Tb_filter_Obs[:, 2] = GMI_low_bin[:, 5]   # 37V
Tb_filter_Obs[:, 3] = GMI_low_bin[:, 7]   # 89V
Tb_filter_Obs[:, 4] = GMI_high_bin[:, 0]   # 166V
Tb_filter_Obs[:, 5] = GMI_high_bin[:, 2]   # 183+3V
Tb_filter_Obs[:, 6] = GMI_high_bin[:, 3]   # 183+7V

# Tb_filter_Sim
Tb_filter_Sim = np.zeros((nbins, 7, mie_tables),dtype= np.float64)
Tb_filter_Sim[:, 0, :] = Tb_low_bin[:, 2, :]   # 19V
Tb_filter_Sim[:, 1, :] = Tb_low_bin[:, 4, :]   # 23V
Tb_filter_Sim[:, 2, :] = Tb_low_bin[:, 5, :]   # 37V
Tb_filter_Sim[:, 3, :] = Tb_low_bin[:, 7, :]   # 89V
Tb_filter_Sim[:, 4, :] = Tb_high_bin[:, 0, :]   # 166V
Tb_filter_Sim[:, 5, :] = Tb_high_bin[:, 2, :]   # 183+3V
Tb_filter_Sim[:, 6, :] = Tb_high_bin[:, 3, :]   # 183+7V

cmap=ez_color.cmap(r=0)
# Create Spatial Grids
xedges = np.linspace(100, 300, num=81)
yedges = np.linspace(0, 7, num= 8)

#fig = plt.figure(figsize=(27,24))
fig = plt.figure(figsize=(27,12))
ax0  = plt.subplot(2,3,1)
cs0 = ax0.imshow(np.transpose(Tb_filter_Obs), interpolation='none', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
#ax0.set_title(' (a) Observed Histogram')
ax0.set_xlabel('Tb')
ax0.set_ylabel('channel')
ax0.set_ylim(0, 7)
ax0.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V','166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax0.set_yticks([0.5,1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax0.set_yticklabels(channels)
cbar =fig.colorbar(cs0, ax=ax0)
cs0.set_clim(1e0, 1e5)

ax1 = plt.subplot(2, 3, 2)
cs1 = ax1.imshow(np.transpose(Tb_filter_Sim[:, :, 0]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
#ax1.set_title(' (b) Simulation with Long Column')
ax1.set_xlabel('Tb')
#ax1.set_ylabel('channel')
ax1.set_ylim(0, 7)
ax1.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax1.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax1.set_yticklabels(channels)
cbar = fig.colorbar(cs1, ax=ax1)
cs1.set_clim(1e0, 1e5)

ax2 = plt.subplot(2, 3, 3)
cs2 = ax2.imshow(np.transpose(Tb_filter_Sim[:, :, 3]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
#ax2.set_title(' (c) Simulation with Thickplate')
ax2.set_xlabel('Tb')
#ax2.set_ylabel('channel')
ax2.set_ylim(0, 7)
ax2.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax2.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax2.set_yticklabels(channels)
cbar = fig.colorbar(cs2, ax=ax2)
cs2.set_clim(1e0, 1e5)

ax3 = plt.subplot(2, 3, 4)
cs3 = ax3.imshow(np.transpose(Tb_filter_Sim[:, :, 4]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
#ax3.set_title(' (d) Simulation with thin plate')
ax3.set_xlabel('Tb')
ax3.set_ylabel('channel')
ax3.set_ylim(0, 7)
ax3.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax3.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax3.set_yticklabels(channels)
cbar = fig.colorbar(cs3, ax=ax3)
cs3.set_clim(1e0, 1e5)

ax4 = plt.subplot(2, 3, 5)
cs4 = ax4.imshow(np.transpose(Tb_filter_Sim[:, :, 5]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
#ax4.set_title(' (e) Simulation with 3-bullet rosette')
ax4.set_xlabel('Tb')
#ax4.set_ylabel('channel')
ax4.set_ylim(0, 7)
ax4.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax4.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax4.set_yticklabels(channels)
cbar = fig.colorbar(cs4, ax=ax4)
cs4.set_clim(1e0, 1e5)

ax5 = plt.subplot(2, 3, 6)
cs5 = ax5.imshow(np.transpose(Tb_filter_Sim[:, :, 9]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
#ax5.set_title(' (f) Simulation with Sector Snowflake')
ax5.set_xlabel('Tb')
#ax5.set_ylabel('channel')
ax5.set_ylim(0, 7)
ax5.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax5.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax5.set_yticklabels(channels)
cbar = fig.colorbar(cs5, ax=ax5)
cs5.set_clim(1e0, 1e5)


plt.savefig('D:/Python_processing/cyclone/results/histogram/histogram_2.5bin_test.png', dpi=300, bbox_inches='tight')
plt.show()

ax6 = plt.subplot(4, 3, 7)
cs6 = ax6.imshow(np.transpose(Tb_filter_Sim[:, :, 13]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
ax6.set_title(' (g) Simulation with 6-bullet rosette')
ax6.set_xlabel('Tb')
ax6.set_ylabel('channel')
ax6.set_ylim(0, 7)
ax6.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax6.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax6.set_yticklabels(channels)
cbar = fig.colorbar(cs6, ax=ax6)
cs6.set_clim(1e0, 1e5)

ax7 = plt.subplot(4, 3, 8)
cs7 = ax7.imshow(np.transpose(Tb_filter_Sim[:, :, 16]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
ax7.set_title(' (h) Simulation with Icon CloudIce')
ax7.set_xlabel('Tb')
#ax5.set_ylabel('channel')
ax7.set_ylim(0, 7)
ax7.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax7.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax7.set_yticklabels(channels)
cbar = fig.colorbar(cs7, ax=ax7)
cs7.set_clim(1e0, 1e5)

ax8 = plt.subplot(4, 3, 9)
cs8 = ax8.imshow(np.transpose(Tb_filter_Sim[:, :, 18]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
ax8.set_title(' (i) Simulation with EvansSnow Aggregate')
ax8.set_xlabel('Tb')
#ax5.set_ylabel('channel')
ax8.set_ylim(0, 7)
ax8.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax8.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax8.set_yticklabels(channels)
cbar = fig.colorbar(cs8, ax=ax8)
cs8.set_clim(1e0, 1e5)

ax9 = plt.subplot(4, 3, 10)
cs9 = ax9.imshow(np.transpose(Tb_filter_Sim[:, :, 23]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
ax9.set_title(' (j) Simulation with Icon Snow')
ax9.set_xlabel('Tb')
ax9.set_ylabel('channel')
ax9.set_ylim(0, 7)
ax9.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax9.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax9.set_yticklabels(channels)
cbar = fig.colorbar(cs9, ax=ax9)
cs9.set_clim(1e0, 1e5)

ax10 = plt.subplot(4, 3, 11)
cs10 = ax10.imshow(np.transpose(Tb_filter_Sim[:, :, 24]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
ax10.set_title(' (k) Simulation with Icon Hail')
ax10.set_xlabel('Tb')
#ax5.set_ylabel('channel')
ax10.set_ylim(0, 7)
ax10.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax10.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax10.set_yticklabels(channels)
cbar = fig.colorbar(cs10, ax=ax10)
cs10.set_clim(1e0, 1e5)

ax11 = plt.subplot(4, 3, 12)
cs11 = ax11.imshow(np.transpose(Tb_filter_Sim[:, :, 25]), interpolation='nearest', origin='low', \
                 extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto', norm=LogNorm())
ax11.set_title(' (l) Simulation with Gem Graupel')
ax11.set_xlabel('Tb')
ax11.set_ylim(0, 7)
ax11.set_xlim(100, 300)
channels = ('19V', '23V', '37V', '89V', '166V', '$183\pm\ 3V$', '$183\pm\ 7V$')
ax11.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
ax11.set_yticklabels(channels)
cbar = fig.colorbar(cs11, ax=ax11)
cs11.set_clim(1e0, 1e5)

plt.savefig('D:/Python_processing/cyclone/results/histogram/histogram_2.5bin.png', dpi=300, bbox_inches='tight')
plt.show()
#plt.close()


