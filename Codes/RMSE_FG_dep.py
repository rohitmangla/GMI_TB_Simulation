import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 15
matplotlib.rcParams['legend.fontsize'] =14
matplotlib.rcParams['font.size']= 20
matplotlib.rcParams['axes.labelsize']= 14
matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 17
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.4
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
GMI_low_final = np.zeros((1, channels_low), dtype=np.float64)
Sim_low_final = np.zeros((1, channels, mie_tables), dtype=np.float64)

for ifile in FILES:
    print(ifile)
    ncfile = Dataset(ifile,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_high = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    Simulated_Tb_low  = ncfile.variables['Simulated_Tb_low'][:]  # low frequency Simulation
    GMI_Tb_low        = ncfile.variables['GMI_gridded_low'][:]  # low frequency GMI Tb
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
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
    Sim_high     = np.squeeze(Simulated_Tb_high[ind_high, :, :])  # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)
    Sim_low = np.squeeze(Simulated_Tb_low[ind_low, :, :])  # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)

GMI_high_final= GMI_high_final[1:]
Sim_high_final= Sim_high_final[1:]
GMI_low_final= GMI_low_final[1:]
Sim_low_final= Sim_low_final[1:]

n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan

# calculate RMSE
rmse_high         = np.zeros((channels_high, mie_tables), dtype=np.float64)
rmse_low         = np.zeros((channels_low, mie_tables), dtype=np.float64)

for itable in range (mie_tables):
    for ichannel in range (channels_high):
        diff_high = GMI_high_final[:, ichannel] - Sim_high_final[:, ichannel+9, itable]
        diff_high = diff_high[~np.isnan(diff_high)]
        diff_sqr_high = diff_high**2
        Sum_val_high = np.sum(diff_sqr_high)
        n = diff_high.shape[0]
        rmse_high[ichannel, itable] = np.sqrt(Sum_val_high/n)

for itable in range (mie_tables):
    for ichannel in range (channels_low):
        diff_low = GMI_low_final[:, ichannel] - Sim_low_final[:, ichannel, itable]
        diff_low = diff_low[~np.isnan(diff_low)]
        diff_sqr_low = diff_low**2
        Sum_val_low = np.sum(diff_sqr_low)
        n = diff_low.shape[0]
        rmse_low[ichannel, itable] = np.sqrt(Sum_val_low/n)

rmse_filter  = np.zeros((7, 26), dtype = np.float64)
rmse_filter[0, :] = rmse_low[2, :]
rmse_filter[1, :] = rmse_low[4, :]
rmse_filter[2, :] = rmse_low[5, :]
rmse_filter[3, :] = rmse_low[7, :]
rmse_filter[4, :] = rmse_high[0, :]
rmse_filter[5, :] = rmse_high[2, :]
rmse_filter[6, :] = rmse_high[3, :]

rmse_mean = np.nanmean (rmse_filter, axis=0)
x_axis_labels = ['long\ncolumn', 'short\ncolumn', 'block\ncolumn','thick\nplate', 'thin\nplate', '3-bullet\nrosette', \
                 '4-bullet\nrosette', '5-bullet\nrosette', '6-bullet\nrosette', 'sector\nsnowflake', 'dendrite\nsnowflake',\
                 'PlateType1', 'ColumnType1', '6-Bullet\nRosette', 'Perpendicular\n4-Bullet Rosette', 'Flat3-\nBulletrosette',\
                 'IconCloudIce', 'Sector\nSnowflake', 'EvansSnow\nAggregate', '8-column\nAggregate','Large-Plate\nAggregate', \
                 'LargeColumn\nAggregate', 'LargeBlock\nAggregate', 'IconSnow', 'IconHail','GemGraupel'] # labels for x-axis

x_val = np.linspace(1,26, 26)
fig= plt.figure(figsize=(17,6))
plt.plot(x_val, rmse_low[2, :], dashes=[3, 10, 1, 10], color= 'black', label= '19 V')
plt.plot(x_val, rmse_low[4, :], dashes=[1, 10, 1,10], color= 'black', label= '23 V')
plt.plot(x_val, rmse_low[5, :], dashes=[3, 5, 1, 5], color= 'black', label= '37 V')
plt.plot(x_val, rmse_low[7, :], dashes=[10, 5, 10, 5], color= 'black', label= '89 V')
plt.plot(x_val, rmse_high[0, :], 'k--', label= '166 V')
plt.plot(x_val, rmse_high[2, :], 'k-.', label= '$183\pm\ 3V$')
plt.plot(x_val, rmse_high[3, :], 'k:', label= '$183\pm\ 7V$')
plt.plot(x_val, rmse_mean, 'k', label= 'Mean')

plt.ylim([2, 35])
plt.xticks(x_val, x_axis_labels, rotation=45, fontname="Times New Roman",fontweight="bold")
plt.ylabel("FG departure RMS (K)")
plt.ylabel("FG departure RMS (K)")
plt.legend(loc= 'upper right', ncol=2)
plt.savefig('D:/Python_processing/cyclone/results/rmse_FGdep.png', dpi=300, bbox_inches= 'tight')
