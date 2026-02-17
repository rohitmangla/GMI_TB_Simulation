import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 38
matplotlib.rcParams['axes.labelsize']= 22
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.25
matplotlib.rcParams['figure.subplot.wspace'] = 0.22
matplotlib.rcParams['lines.linewidth'] = 4

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)

for ifile in FILES:
    print(ifile)
    ncfile = Dataset(ifile,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_high = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
    lsm               = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
    # Land masking
    ind_land = np.where(lsm == 1)
    Simulated_Tb_high[ind_land, :, :] = np.nan
    GMI_Tb_high[ind_land, :] = np.nan
    ind_high = np.where(GMI_Tb_high[:, 0] > 1)
    lat_high = lat[ind_high]
    lon_high = lon[ind_high]

    GMI_high = np.squeeze(GMI_Tb_high[ind_high, :])  # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_high = np.squeeze(Simulated_Tb_high[ind_high, :, :])  # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)

GMI_high_final[0, :] = np.nan
Sim_high_final[0, :, :] = np.nan

bin_start = 100
bin_end   = 300
nbins     = 40
dbin = 5
log_ratio_high = np.zeros((nbins,mie_tables, channels_high), dtype=np.float64)
histogram_ratio_high = np.zeros((nbins,mie_tables,channels_high), dtype=np.float64)
count_sim_final= np.zeros((nbins,mie_tables, channels_high), dtype=np.float64)
count_obs_final= np.zeros((nbins,mie_tables, channels_high), dtype=np.float64)

for ibin in range (nbins):
    for itable in range (mie_tables):
        for ichannel in range(channels_high):
            ind_obs= np.where((GMI_high_final[:, ichannel]>= bin_start)*(GMI_high_final[:, ichannel]<bin_start+dbin))
            ind_sim= np.where((Sim_high_final[:, ichannel+9, itable]>= bin_start)*(Sim_high_final[:, ichannel+9, itable]<bin_start+dbin))
            obs = GMI_high_final[:, ichannel][ind_obs]
            count_obs = np.count_nonzero(obs)
            count_obs_final[ibin, itable, ichannel] = count_obs
            sim = Sim_high_final[:, ichannel+9, itable][ind_sim]
            count_sim = np.count_nonzero(sim)
            count_sim_final[ibin, itable, ichannel] = count_sim
            if (count_obs ==0) | (count_sim==0):
                count_obs= 0.1
                count_sim= 0.1
            log_ratio_high[ibin, itable, ichannel] = abs(np.log10(count_sim/count_obs))
            histogram_ratio_high[ibin, itable, ichannel] = np.log10(count_sim/count_obs)
    bin_start = bin_start+dbin

fig= plt.figure(figsize=(26,17))
ax1= plt.subplot(2,2,1)
x1= np.linspace(100, 300, 81)
ax1.set_title('(a) 166V')
ax1.plot(x1[1:], histogram_ratio_high[:,4, 0], 'r', label='Thinplate')
ax1.plot(x1[1:], histogram_ratio_high[:,5, 0], 'r--', label='3-bullet rosette')
ax1.plot(x1[1:], histogram_ratio_high[:,9, 0], 'r:', label='Sector Snowflake')
ax1.plot(x1[1:], histogram_ratio_high[:,16, 0],color='b', linestyle='-', label='IconCLoudIce')
ax1.plot(x1[1:], histogram_ratio_high[:,19, 0],color='b', linestyle='--', label='8-column Aggregate')
ax1.plot(x1[1:], histogram_ratio_high[:,24, 0], color ='b', linestyle = ':', label='Icon Hail')

ax1.set_ylim([-2.5,2])
ax1.set_yticks([-2, -1, 0, 1, 2])
ax1.set_xticks(np.arange(100, 325, step=25))
ax1.set_xticklabels(['100', '125', '150', '175', '200', '225', '250', '275', '300'], fontsize=16)
ax1.set_ylabel("$log_{10}$ (histogram ratio)")
ax1.set_xlabel("Brightness Temperature (K)")
ax1.legend()
ax3= plt.subplot(2,2,2)
ax3.set_title('(b) $183\pm\ 3V$')
ax3.plot(x1[1:], histogram_ratio_high[:,4, 2], 'r', label='Thinplate')
ax3.plot(x1[1:], histogram_ratio_high[:,5, 2], 'r--', label='3-bullet rosette')
ax3.plot(x1[1:], histogram_ratio_high[:,9, 2], 'r:', label='Sector Snowflake')
ax3.plot(x1[1:], histogram_ratio_high[:,16, 2],color='b', linestyle='-', label='IconCLoudIce')
ax3.plot(x1[1:], histogram_ratio_high[:,19, 2],color='b', linestyle='--', label='8-column Aggregate')
ax3.plot(x1[1:], histogram_ratio_high[:,24, 2], color ='b', linestyle = ':', label='Icon Hail')
ax3.set_ylim([-2.5,2])
ax3.set_yticks([-2, -1, 0, 1, 2])
ax3.set_xticks(np.arange(100, 325, step=25))
ax3.set_xticklabels(['100', '125', '150', '175', '200', '225', '250', '275', '300'], fontsize=16)
ax3.set_ylabel("$log_{10}$ (histogram ratio)")
ax3.set_xlabel("Brightness Temperature (K)")


ax4= plt.subplot(2,2,3)
ax4.set_title('(c) $183\pm\ 7V$')
ax4.plot(x1[1:], histogram_ratio_high[:,4, 3], 'r', label='Thinplate')
ax4.plot(x1[1:], histogram_ratio_high[:,5, 3], 'r--', label='3-bullet rosette')
ax4.plot(x1[1:], histogram_ratio_high[:,9, 3], 'r:', label='Sector Snowflake')
ax4.plot(x1[1:], histogram_ratio_high[:,16, 3],color='b', linestyle='-', label='IconCLoudIce')
ax4.plot(x1[1:], histogram_ratio_high[:,19, 3],color='b', linestyle='--', label='8-column Aggregate')
ax4.plot(x1[1:], histogram_ratio_high[:,24, 3], color ='b', linestyle = ':', label='Icon Hail')
ax4.set_ylim([-2.5,2])
ax4.set_yticks([-2, -1, 0, 1, 2])
ax4.set_xticks(np.arange(100, 325, step=25))
ax4.set_xticklabels(['100', '125', '150', '175', '200', '225', '250', '275', '300'], fontsize=16)
ax4.set_ylabel("$log_{10}$ (histogram ratio)")
ax4.set_xlabel("Brightness Temperature (K)")

ax5= plt.subplot(2,2,4)
ax5.set_title('(d) counts')
ax5.plot(x1[1:], count_sim_final[:,4, 3], 'r', label='Thinplate')
ax5.plot(x1[1:], count_sim_final[:,5, 3], 'r--', label='3-bullet rosette')
ax5.plot(x1[1:], count_sim_final[:,9, 3], 'r:', label='Sector Snowflake')
ax5.plot(x1[1:], count_sim_final[:,16, 3],color='b', linestyle='-', label='IconCLoudIce')
ax5.plot(x1[1:], count_sim_final[:,19, 3],color='b', linestyle='--', label='8-column Aggregate')
ax5.plot(x1[1:], count_sim_final[:,24, 3], color ='b', linestyle = ':', label='Icon Hail')
ax5.plot(x1[1:], count_obs_final[:,25, 3], 'k', label='Obs')
plt.yscale("log")
ax5.set_xticks(np.arange(100, 325, step=25))
ax5.set_xticklabels(['100', '125', '150', '175', '200', '225', '250', '275', '300'], fontsize=16)
ax5.set_ylabel("Number in bins")
ax5.set_xlabel("Brightness Temperature (K)")

plt.savefig('D:/Python_processing/cyclone/results/histogram_ratio.png', dpi=300, bbox_inches= 'tight')



