import numpy as np
from netCDF4 import Dataset
import glob
import pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors
import matplotlib as mpl

matplotlib.rcParams['axes.titlesize'] = 18
matplotlib.rcParams['legend.fontsize'] = 14.
matplotlib.rcParams['font.size'] = 22.
matplotlib.rcParams['axes.labelsize'] = 19.
matplotlib.rcParams['xtick.labelsize'] = 18.
matplotlib.rcParams['ytick.labelsize'] = 18.
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")
channels_high= 4
channels_low= 9
channels = 13
mie_tables= 26
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
GMI_low_final  = np.zeros((1, channels_low), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_low_final  = np.zeros((1, channels, mie_tables), dtype=np.float64)

for ifile in FILES:
    print(ifile)
    ncfile = Dataset(ifile,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_high = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    Simulated_Tb_low  = ncfile.variables['Simulated_Tb_low'][:]  # low frequency Simulation
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
    GMI_Tb_low        = ncfile.variables['GMI_gridded_low'][:]  # high frequency GMI Tb
    lsm               = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
    # Land masking
    ind_land = np.where(lsm == 1)
    Simulated_Tb_high[ind_land, :, :] = np.nan
    Simulated_Tb_low[ind_land, :, :] = np.nan
    GMI_Tb_high[ind_land, :] = np.nan
    GMI_Tb_low[ind_land, :] = np.nan
    ind_high = np.where(GMI_Tb_high[:, 0] > 1)
    lat_high = lat[ind_high]
    lon_high = lon[ind_high]
    ind_low = np.where(GMI_Tb_low[:, 0] > 1)
    lat_low = lat[ind_low]
    lon_low = lon[ind_low]

    GMI_high = np.squeeze(GMI_Tb_high[ind_high, :])  # (profiles, channels_high)
    GMI_low = np.squeeze(GMI_Tb_low[ind_low, :])  # (profiles, channels_low)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    Sim_high = np.squeeze(Simulated_Tb_high[ind_high, :, :])  # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)
    Sim_low  = np.squeeze(Simulated_Tb_low[ind_low, :, :])  # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)

GMI_high_final[0, :] = np.nan
Sim_high_final[0, :, :] = np.nan
GMI_low_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan

nbins = 80
delta_bin= 2.5
obs_high = np.zeros((nbins,channels_high), dtype= np.float64)
sim_high = np.zeros((nbins,channels_high, mie_tables), dtype= np.float64)
obs_low = np.zeros((nbins,channels_low), dtype= np.float64)
sim_low = np.zeros((nbins,channels_low, mie_tables), dtype= np.float64)

# for high range frequencies:
bin_start = 100
for ibin in range (nbins):
    for ichannel in range (channels_high):
        for itable in range (mie_tables):
            ind_obs= np.where((GMI_high_final[:, ichannel]>= bin_start)*(GMI_high_final[:, ichannel]<bin_start+delta_bin))
            ind_sim= np.where((Sim_high_final[:, ichannel+9, itable]>= bin_start)*(Sim_high_final[:, ichannel+9, itable]<bin_start+delta_bin))
            obs = GMI_high_final[:, ichannel][ind_obs]
            obs_count = np.count_nonzero(obs)
            obs_high[ibin, ichannel] = obs_count
            sim = Sim_high_final[:, ichannel+9, itable][ind_sim]
            sim_count = np.count_nonzero(sim)
            sim_high[ibin,ichannel,itable] = sim_count
    bin_start= bin_start+delta_bin

# for low range frequencies:
bin_start = 100
for ibin in range (nbins):
    for ichannel in range (channels_low):
        for itable in range (mie_tables):
            ind_obs= np.where((GMI_low_final[:, ichannel]>= bin_start)*(GMI_low_final[:, ichannel]<bin_start+delta_bin))
            ind_sim= np.where((Sim_low_final[:, ichannel, itable]>= bin_start)*(Sim_low_final[:, ichannel, itable]<bin_start+delta_bin))
            obs = GMI_low_final[:, ichannel][ind_obs]
            obs_count = np.count_nonzero(obs)
            obs_low[ibin, ichannel] = obs_count
            sim = Sim_low_final[:, ichannel, itable][ind_sim]
            sim_count = np.count_nonzero(sim)
            sim_low[ibin,ichannel,itable] = sim_count
    bin_start= bin_start+delta_bin

# Percentile Image high
Q1_high  = np.zeros ((nbins,channels_high), dtype= np.float64)
Q2_high  = np.zeros ((nbins,channels_high), dtype= np.float64)
Q3_high  = np.zeros ((nbins,channels_high), dtype= np.float64)
Q4_high  = np.zeros ((nbins,channels_high), dtype= np.float64)

# Percentile Image low
Q1_low  = np.zeros ((nbins,channels_low), dtype= np.float64)
Q2_low  = np.zeros ((nbins,channels_low), dtype= np.float64)
Q3_low  = np.zeros ((nbins,channels_low), dtype= np.float64)
Q4_low  = np.zeros ((nbins,channels_low), dtype= np.float64)

for ibin in range (nbins):
    for ichannel in range (channels_high):
        Q1_high[ibin, ichannel] = np.percentile(sim_high[ibin, ichannel, :], 5)
        Q2_high[ibin,ichannel]  = np.percentile(sim_high[ibin, ichannel, :], 25)
        Q3_high[ibin,ichannel]  = np.percentile(sim_high[ibin, ichannel, :], 75)
        Q4_high[ibin,ichannel]  = np.percentile(sim_high[ibin, ichannel, :], 95)
count_high   = np.zeros((nbins, channels_high), dtype= np.float64)-999.0

for ibin in range (nbins):
    for ichannel in range (channels_low):
        Q1_low[ibin, ichannel] = np.percentile(sim_low[ibin, ichannel, :], 5)
        Q2_low[ibin,ichannel]  = np.percentile(sim_low[ibin, ichannel, :], 25)
        Q3_low[ibin,ichannel]  = np.percentile(sim_low[ibin, ichannel, :], 75)
        Q4_low[ibin,ichannel]  = np.percentile(sim_low[ibin, ichannel, :], 95)
count_low= np.zeros((nbins, channels_low), dtype= np.float64)-999.0

for icount in range (nbins):
    for jcount in range (channels_high):
        if obs_high[icount,jcount]==0:
            count_high[icount, jcount] = np.nan
        elif obs_high[icount,jcount]>0 and obs_high[icount,jcount]<=Q1_high[icount,jcount]:
            count_high[icount,jcount]= 1
        elif obs_high[icount,jcount]>Q1_high[icount,jcount] and obs_high[icount,jcount]<=Q2_high[icount,jcount]:
            count_high[icount,jcount]= 2
        elif obs_high[icount,jcount]>Q2_high[icount,jcount] and obs_high[icount,jcount]<=Q3_high[icount,jcount]:
            count_high[icount,jcount]= 3
        elif obs_high[icount,jcount]>Q3_high[icount,jcount] and obs_high[icount,jcount]<=Q4_high[icount,jcount]:
            count_high[icount,jcount]= 4
        else:
            count_high[icount, jcount]= 5

count_high[count_high ==-999.0]= np.nan
count_high = np.transpose(count_high)

for icount in range (nbins):
    for jcount in range (channels_low):
        if obs_low[icount,jcount]==0:
            count_low[icount, jcount] = np.nan
        elif obs_low[icount,jcount]>0 and obs_low[icount,jcount]<=Q1_low[icount,jcount]:
            count_low[icount,jcount]= 1
        elif obs_low[icount,jcount]>Q1_low[icount,jcount] and obs_low[icount,jcount]<=Q2_low[icount,jcount]:
            count_low[icount,jcount]= 2
        elif obs_low[icount,jcount]>Q2_low[icount,jcount] and obs_low[icount,jcount]<=Q3_low[icount,jcount]:
            count_low[icount,jcount]= 3
        elif obs_low[icount,jcount]>Q3_low[icount,jcount] and obs_low[icount,jcount]<=Q4_low[icount,jcount]:
            count_low[icount,jcount]= 4
        else:
            count_low[icount, jcount]= 5

count_low[count_low ==-999.0]= np.nan
count_low = np.transpose(count_low)

#count_filter = np.zeros((4,nbins), dtype=np.float64)
count_filter = np.zeros((7,nbins), dtype=np.float64)
count_filter[0,:] = count_low[2,:]
count_filter[1,:] = count_low[4,:]
count_filter[2,:] = count_low[5,:]
count_filter[3,:] = count_low[7,:]
count_filter[4,:] = count_high[0,:]
count_filter[5,:] = count_high[2,:]
count_filter[6,:] = count_high[3,:]

fig  = plt.figure(figsize=(14,5))
cmap = mpl.colors.ListedColormap(['red','orange','yellow','#7B68EE', 'blue'])
bounds = [1,2,3,4,5,6]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
ax1  = plt.subplot(1,1,1)
#im_count = ax1.imshow(count_filter,interpolation= 'none',cmap = cmap, norm=norm, origin= 'low',aspect= 'auto')
im_count = ax1.imshow(count_filter,interpolation= 'none',cmap = cmap, norm=norm, origin= 'lower',aspect= 'auto')
ax1.set_xlabel('Brightness Temperature (K)')
ax1.set_ylabel('Channel name')
ax1.set_ylim(-0.5,6.5)
ax1.set_yticks(np.arange(0, 7, step=1))
ax1.set_yticklabels(labels= ['19V','23V','37V','89V','166V','$183\pm\ 3V$', '$183\pm\ 7V$'])
ax1.set_xlim(0, 80)
ax1.set_xticks(np.arange(0,81, step=10))
ax1.set_xticklabels(labels= ['100','125', '150','175', '200', '225', '250','275', '300'])
divider = make_axes_locatable(ax1)
cax = divider.append_axes ("right", size = "3%", pad = 0.05)
cbar= plt.colorbar(im_count, ticks =bounds, cax=  cax, extend = 'both')
cbar.ax.set_yticklabels(['0%','5%', '25%', '75%', '95%', '100%'])
plt.tight_layout()
plt.savefig('D:/Python_processing/cyclone/results/percentile_plot_2_12hrs.png', dpi=300, bbox_inches= 'tight')
plt.show()



