import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 24
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.35
matplotlib.rcParams['lines.linewidth'] = 2.5

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
PATH_Clr = 'D:/Python_processing/cyclone/forecast_clear_sky/'
FILES= glob.glob(PATH+"*.nc")
FILES_Clr= glob.glob(PATH_Clr+"*.nc")
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_low_final = np.zeros((1, channels_low), dtype=np.float64)
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_low_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)

Sim_low_clr_final = np.zeros((1, channels), dtype=np.float64)
Sim_high_clr_final = np.zeros((1, channels), dtype=np.float64)

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

    ind_high = np.where(GMI_Tb_low[:, 0] > 1)
    lat_high = lat[ind_high]
    lon_high = lon[ind_high]

    GMI_low   = np.squeeze(GMI_Tb_low[ind_high, :])              # (profiles, channels_low)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind_high, :])             # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind_high, :, :])     # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind_high, :, :])    # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)

GMI_low_final[0, :] = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
Sim_high_final[0, :, :] = np.nan

# import Clear sky radiances
for ifile_clr in FILES_Clr:
    print(ifile_clr)
    ncfile = Dataset(ifile_clr,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_low_clr  = ncfile.variables['Simulated_Tb_low'][:]   # low frequency Simulation
    Simulated_Tb_high_clr = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    GMI_Tb_low        = ncfile.variables['GMI_gridded_low'][:]   # low frequency GMI Tb
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
    lsm               = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
    # Land masking
    ind_land = np.where(lsm ==1)
    Simulated_Tb_low_clr[ind_land, :] = np.nan
    Simulated_Tb_high_clr[ind_land, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind_high = np.where(GMI_Tb_low[:, 0]>1)
    Sim_low_clr   = np.squeeze(Simulated_Tb_low_clr[ind_high, :])     # (profiles, channels)
    Sim_low_clr_final = np.concatenate((Sim_low_clr_final, Sim_low_clr), 0)
    Sim_high_clr  = np.squeeze(Simulated_Tb_high_clr[ind_high, :])    # (profiles, channels)
    Sim_high_clr_final = np.concatenate((Sim_high_clr_final, Sim_high_clr), 0)

Sim_low_clr_final[0, :]  = np.nan
Sim_high_clr_final[0, :] = np.nan

# polarization differences at 183 Ghz
PD_sim = (Sim_low_final[:, 5, 4]- Sim_low_final[:, 6, 4])/(Sim_low_clr_final[:, 5]-Sim_low_clr_final[:, 6])
CA_sim = 1-PD_sim
PD_obs = (GMI_low_final[:, 5]- GMI_low_final[:, 6])/(Sim_low_clr_final[:, 5]-Sim_low_clr_final[:, 6])
CA_obs = 1-PD_obs
CA_avg = (CA_obs+CA_sim)/2

cloud_amount_std = np.zeros((22, 3, 13), dtype=np.float64)
for i in range (9):
    # for low frequency
    FG = GMI_low_final[:, i]-Sim_low_final[:,i,4]
    cmin = -0.075
    for ibin in range (22):
        ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+0.05))
        FG_sub = FG[ind]
        count = np.count_nonzero(FG_sub)
        cloud_amount_std[ibin, 0, i] = (cmin+cmin+0.05)/2
        cloud_amount_std[ibin, 1, i] = np.nanstd(FG_sub)
        cloud_amount_std[ibin, 2, i] = count
        cmin = cmin+0.05
for i in range (4):
    # for high frequency
    FG = GMI_high_final[:, i]-Sim_high_final[:,i+9,4]
    cmin = -0.075
    for ibin in range (22):
        ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+0.05))
        FG_sub = FG[ind]
        count = np.count_nonzero(FG_sub)
        cloud_amount_std[ibin, 0, i+9] = (cmin+cmin+0.05)/2
        cloud_amount_std[ibin, 1, i+9] = np.nanstd(FG_sub)
        cloud_amount_std[ibin, 2, i+9] = count
        cmin = cmin+0.05
# calculate percentage
total_samples = np.sum(cloud_amount_std[:, 2,2])
cloud_per= (cloud_amount_std[:,2,2]*100)/total_samples
fig = plt.figure(figsize=(8, 5))
ax1= plt.subplot(1,1,1)
ax1.plot(cloud_amount_std[0:18,0, 2], cloud_amount_std[0:18,1, 2], 'k', label = '19 V')
ax1.plot(cloud_amount_std[0:18,0, 4], cloud_amount_std[0:18,1, 4], 'k--', label = '23 V')
ax1.plot(cloud_amount_std[0:18,0, 5], cloud_amount_std[0:18,1, 5], 'k:', label = '37 V')
ax1.set_xlabel('$ C_{37avg} $')
ax1.set_ylabel('Std (K)')
ax1.set_ylim([0, 40])
ax1.legend(loc='upper left', ncol=2)

ax2 = ax1.twinx()
ax2.plot(cloud_amount_std[0:18,0, 2], cloud_per[0:18], color = 'tan', linestyle='--', label= 'number of samples\n(% total)')
ax2.set_ylabel("Number of Samples(% total)")
ax2.legend(loc='upper right')
ax2.set_ylim([0, 60])
plt.savefig('D:/Python_processing/cyclone/results/error_model/std_dev_low.png', dpi=300, bbox_inches='tight')

