import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 14
matplotlib.rcParams['legend.fontsize'] =12
matplotlib.rcParams['font.size']= 22
matplotlib.rcParams['axes.labelsize']= 12
matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
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
lat_low_final = np.zeros((1), dtype=np.float64)
lon_low_final = np.zeros((1), dtype=np.float64)

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
    lat_low_final = np.concatenate((lat_low, lat_low_final), 0)
    lon_low_final = np.concatenate((lon_low, lon_low_final), 0)

GMI_low_final[0, :] = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
Sim_high_final[0, :, :] = np.nan
lat_low_final[0] = np.nan
lon_low_final[0] = np.nan
n= np.where((GMI_low_final[:,5]>160)*(GMI_low_final[:,5]<200))
GMI_low_final[n, 5] = np.nan


# density plot with shade
fig = plt.figure(figsize=(12, 17))
ax0 = plt.subplot(4,2,1)
ax0.set_title("(a) 19 V")
kwargs = dict(histtype='stepfilled', alpha=0.3, normed=True, bins=40)
ax0.hist(GMI_low_final[:,2], **kwargs)
ax0.hist(Sim_low_final[:,2,4], **kwargs)

#ax0 = sns.kdeplot(GMI_low_final[:, 2], shade=True, color="r", label = 'Observed')
#ax0 = sns.kdeplot(Sim_low_final[:, 2, 9], shade=True,color="b", label = 'Simulated')
ax0.set_xlim([100, 300])
ax0.set(yscale="log")
ax0.set_xlabel('Brightness Temperature (K)')
ax0.set_ylabel('PDF')
ax0.legend(loc= 'upper left')

ax1 = plt.subplot(4,2,2)
ax1.set_title("(b) 23 V")
ax1 = sns.kdeplot(GMI_low_final[:, 4], shade=True, color="r", label = 'Observed')
ax1 = sns.kdeplot(Sim_low_final[:, 4, 4], shade=True,color="b", label = 'Simulated')
ax1.set_xlim([100, 300])
ax1.set(yscale="log")
ax1.set_xlabel('Brightness Temperature (K)')
ax1.set_ylabel('PDF')
ax1.legend(loc= 'upper left')

ax2 = plt.subplot(4,2,3)
ax2.set_title("(c) 37 V")
ax2 = sns.kdeplot(GMI_low_final[:, 5], shade=True, color="r", label = 'Observed')
ax2 = sns.kdeplot(Sim_low_final[:, 5, 4], shade=True,color="b", label = 'Simulated')
ax2.set_xlim([100, 300])
ax2.set(yscale="log")
ax2.set_xlabel('Brightness Temperature (K)')
ax2.set_ylabel('PDF')
ax2.legend(loc= 'upper left')

ax3 = plt.subplot(4,2,4)
ax3.set_title("(d) 89 V")
ax3 = sns.kdeplot(GMI_low_final[:, 7], shade=True, color="r", label = 'Observed')
ax3 = sns.kdeplot(Sim_low_final[:, 7, 4], shade=True,color="b", label = 'Simulated')
ax3.set_xlim([100, 300])
ax3.set(yscale="log")
ax3.set_xlabel('Brightness Temperature (K)')
ax3.set_ylabel('PDF')
ax3.legend(loc= 'upper left')

ax4 = plt.subplot(4,2,5)
ax4.set_title("(e) 166 V")
ax4 = sns.kdeplot(GMI_high_final[:, 0], shade=True, color="r", label = 'Observed')
ax4 = sns.kdeplot(Sim_high_final[:, 9, 4], shade=True,color="b", label = 'Simulated')
ax4.set_xlim([100, 300])
ax4.set(yscale="log")
ax4.set_xlabel('Brightness Temperature (K)')
ax4.set_ylabel('PDF')
ax4.legend(loc= 'upper left')

ax5 = plt.subplot(4,2,6)
ax5.set_title("(f) $183\pm\ 3V$")
ax5 = sns.kdeplot(GMI_high_final[:, 2], shade=True, color="r", label = 'Observed')
ax5 = sns.kdeplot(Sim_high_final[:, 11, 4], shade=True, color="b", label = 'Simulated')
ax5.set_xlim([100, 300])
ax5.set(yscale="log")
ax5.set_xlabel('Brightness Temperature (K)')
ax5.set_ylabel('PDF')
ax5.legend(loc= 'upper left')

ax6 = plt.subplot(4,2,7)
ax6.set_title("(g) $183\pm\ 7V$")
#ax6 = sns.kdeplot(GMI_high_final[:, 3], shade=True, color="r", label = 'Observed')
#ax6 = sns.kdeplot(Sim_high_final[:, 12, 9], shade=True, color="b", label = 'Simulated')
kwargs = dict(histtype='stepfilled', alpha=0.3, normed=True, bins=100)

ax6.hist(GMI_high_final[:,3], **kwargs)
ax6.hist(Sim_high_final[:,12,9], **kwargs)

ax6.set_xlim([100, 300])
ax6.set(yscale="log")
ax6.set_xlabel('Brightness Temperature (K)')
ax6.set_ylabel('PDF')
ax6.legend(loc= 'upper left')

plt.savefig('D:/Python_processing/cyclone/results/Observed_sim_hist_sector_snowflake.png', dpi=300, bbox_inches= 'tight')
plt.show()

