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

bin_start = 100
nbins     = 80
dbin= 2.5
log_ratio_high = np.zeros((nbins,mie_tables, channels_high), dtype=np.float64)
histogram_ratio_high = np.zeros((nbins,mie_tables,channels_high), dtype=np.float64)

for ibin in range (nbins):
    for itable in range (mie_tables):
        for ichannel in range(channels_high):
            ind_obs= np.where((GMI_high_final[:, ichannel]>= bin_start) & (GMI_high_final[:, ichannel]<bin_start+dbin))
            ind_sim= np.where((Sim_high_final[:, ichannel+9, itable]>= bin_start) & (Sim_high_final[:, ichannel+9, itable]<bin_start+dbin))
            obs = GMI_high_final[:, ichannel][ind_obs]
            count_obs = np.count_nonzero(obs)
            sim = Sim_high_final[:, ichannel+9, itable][ind_sim]
            count_sim = np.count_nonzero(sim)
            if (count_obs ==0) | (count_sim==0):
                count_obs= 0.1
                count_sim= 0.1
            log_ratio_high[ibin, itable, ichannel] = abs(np.log10(count_sim/count_obs))
            histogram_ratio_high[ibin, itable, ichannel] = np.log10(count_sim/count_obs)
    bin_start = bin_start+dbin

log_ratio_low = np.zeros((nbins,mie_tables, channels_low), dtype=np.float64)
histogram_ratio_low = np.zeros((nbins,mie_tables,channels_low), dtype=np.float64)

bin_start=100
for ibin in range (nbins):
    for itable in range (mie_tables):
        for ichannel in range(channels_low):
            ind_obs= np.where((GMI_low_final[:, ichannel]>= bin_start) & (GMI_low_final[:, ichannel]<bin_start+dbin))
            ind_sim= np.where((Sim_low_final[:, ichannel, itable]>= bin_start) & (Sim_low_final[:, ichannel, itable]<bin_start+dbin))
            obs = GMI_low_final[:, ichannel][ind_obs]
            count_obs = np.count_nonzero(obs)
            sim = Sim_low_final[:, ichannel, itable][ind_sim]
            count_sim = np.count_nonzero(sim)
            if (count_obs ==0) | (count_sim==0):
                count_obs= 0.1
                count_sim= 0.1
            log_ratio_low[ibin, itable, ichannel] = abs(np.log10(count_sim/count_obs))
            histogram_ratio_low[ibin, itable, ichannel] = np.log10(count_sim/count_obs)
    bin_start = bin_start+dbin

h_high = np.zeros((mie_tables,channels_high),  dtype= np.float64)
for itable in range(mie_tables):
    for ichannel in range (channels_high):
        summation = sum(log_ratio_high[:, itable, ichannel])
        total_bins = np.count_nonzero(log_ratio_high[:, itable, ichannel])
        h_high[itable, ichannel] = summation/total_bins

h_low = np.zeros((mie_tables,channels_low),  dtype= np.float64)
for itable in range(mie_tables):
    for ichannel in range (channels_low):
        summation = sum(log_ratio_low[:, itable, ichannel])
        total_bins = np.count_nonzero(log_ratio_low[:, itable, ichannel])
        h_low[itable, ichannel] = summation/total_bins

x_axis_labels = ['long\ncolumn', 'short\ncolumn', 'block\ncolumn','thick\nplate', 'thin\nplate', '3-bullet\nrosette', \
                 '4-bullet\nrosette', '5-bullet\nrosette', '6-bullet\nrosette', 'sector\nsnowflake', 'dendrite\nsnowflake',\
                 'PlateType1', 'ColumnType1', '6-Bullet\nRosette', 'Perpendicular\n4-Bullet Rosette', 'Flat3-\nBulletrosette',\
                 'IconCloudIce', 'Sector\nSnowflake', 'EvansSnow\nAggregate', '8-column\nAggregate','Large-Plate\nAggregate', \
                 'LargeColumn\nAggregate', 'LargeBlock\nAggregate', 'IconSnow', 'IconHail','GemGraupel'] # labels for x-axis

h_filter  = np.zeros((26, 7), dtype = np.float64)
h_filter[:, 0] = h_low[:, 2]
h_filter[:, 1] = h_low[:, 4]
h_filter[:, 2] = h_low[:, 5]
h_filter[:, 3] = h_low[:, 7]
h_filter[:, 4] = h_high[:, 0]
h_filter[:, 5] = h_high[:, 2]
h_filter[:, 6] = h_high[:, 3]

h_mean = np.nanmean (h_filter, axis=1)
x_val = np.linspace(1,26, 26)
fig= plt.figure(figsize=(17,6))
plt.plot(x_val, h_low[:, 2], dashes=[3, 10, 1, 10], color= 'black', label= '19 V')
plt.plot(x_val, h_low[:, 4], dashes=[1, 10, 1,10], color= 'black', label= '23 V')
plt.plot(x_val, h_low[:, 5], dashes=[3, 5, 1, 5], color= 'black', label= '37 V')
plt.plot(x_val, h_low[:, 7], dashes=[10, 5, 10, 5], color= 'black', label= '89 V')
plt.plot(x_val, h_high[:, 0], 'k--', label= '166 V')
plt.plot(x_val, h_high[:, 2], 'k-.', label= '$183\pm\ 3V$')
plt.plot(x_val, h_high[:, 3], 'k:', label= '$183\pm\ 7V$')
plt.plot(x_val, h_mean, 'k', label= 'Mean')

plt.ylim([0, 0.9])
#plt.xticks(x_val, x_axis_labels, rotation=50)
plt.xticks(x_val, x_axis_labels, rotation=45, fontname="Times New Roman",fontweight="bold")
#plt.title("(a) Tb Histogram fit")
plt.ylabel("h-values")
plt.legend(loc= 'upper right', ncol=2)
plt.savefig('D:/Python_processing/cyclone/results/h_values_2.png', dpi=300, bbox_inches= 'tight')




