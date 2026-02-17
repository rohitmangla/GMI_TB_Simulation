import numpy as np
from netCDF4 import Dataset
import glob
from scipy.stats import skew
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 16
matplotlib.rcParams['legend.fontsize'] =8
matplotlib.rcParams['font.size']= 16
matplotlib.rcParams['axes.labelsize']= 14
matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.25
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



# first guess departures
FG_low = GMI_low_final-Sim_low_final[:, 0:9, 4]
FG_high = GMI_high_final-Sim_high_final[:, 9:13, 4]
FG_final = np.concatenate((FG_low, FG_high), 1)
RMSE = np.zeros(13, dtype=np.float64)
# RMSE error
for i in range (13):
    FG_sq = FG_final[:, i]*FG_final[:, i]
    FG_sq = FG_sq[~np.isnan(FG_sq)]
    FG_sq_sum = np.sum(FG_sq)
    num = len(FG_sq)
    MSE = FG_sq_sum/num
    RMSE[i]= np.sqrt(MSE)

# Skewness
Skew = np.zeros(13, dtype=np.float64)
for i in range(13):
    FG = FG_final[:, i]
    FG = FG[~np.isnan(FG)]
    Skew[i] = skew(FG)
Skew_filter = np.zeros((6,1), dtype= np.float64)
RMSE_filter = np.zeros((6,1), dtype= np.float64)

Skew_filter[0] = Skew[2]
Skew_filter[1] = Skew[4]
Skew_filter[2] = Skew[5]
Skew_filter[3] = Skew[7]
Skew_filter[4] = Skew[9]
Skew_filter[5] = Skew[12]

RMSE_filter[0] = RMSE[2]
RMSE_filter[1] = RMSE[4]
RMSE_filter[2] = RMSE[5]
RMSE_filter[3] = RMSE[7]
RMSE_filter[4] = RMSE[9]
RMSE_filter[5] = RMSE[12]

fig= plt.figure(figsize=(10,6))
x_val = np.linspace(1,13, 13)
ax0= plt.subplot(1,2,1)
ax0.set_title("(a) Skewness")
ax0.plot(x_val, Skew, 'k', )
ax0.set_ylim([-2,3])
ax0.set_xticks(np.arange(1, 14, step=1))
#ax0.set_xticklabels(['19V','23V', '37V','89V','166V','$183\pm\ 7V$'])
ax0.set_xticklabels(['10V', '10H', '19V','19H', '23V', '37V', '37H', '89V','89H', '166V','166H', '$183\pm\ 3V$','$183\pm\ 7V$'])
ax0.set_xlabel("channel name")
ax0.set_ylabel("skewness")

ax1= plt.subplot(1,2,2)
ax1.set_title("(b) RMSE")
ax1.plot(x_val, RMSE, 'k', )
ax1.set_ylim([0, 30])
ax1.set_xticks(np.arange(1, 14, step=1))
ax1.set_xticklabels(['10V', '10H', '19V','19H', '23V', '37V', '37H', '89V','89H', '166V','166H', '$183\pm\ 3V$','$183\pm\ 7V$'])
#ax1.set_xticklabels(['19V','23V', '37V','89V','166V','$183\pm\ 7V$'])
ax1.set_xlabel("channel name")
ax1.set_ylabel("RMSE")
plt.savefig('D:/Python_processing/cyclone/results/RMSE_Skew_full.png', dpi=300, bbox_inches= 'tight')

plt.show()