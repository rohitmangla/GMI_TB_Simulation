import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import matplotlib
from scipy.stats import norm

matplotlib.rcParams['axes.titlesize'] = 11
matplotlib.rcParams['legend.fontsize'] =8
matplotlib.rcParams['font.size']= 18
matplotlib.rcParams['axes.labelsize']= 8
matplotlib.rcParams['xtick.labelsize'] = 7
matplotlib.rcParams['ytick.labelsize'] = 7
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.4
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/Forecast_12hrs/All/'
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
FG_low_final  = np.zeros((GMI_low_final.shape[0], 9, mie_tables), dtype=np.float64)
FG_high_final = np.zeros((GMI_high_final.shape[0], 4, mie_tables), dtype=np.float64)
rmse_low          = np.zeros((9, 26), dtype=np.float64)
rmse_high         = np.zeros((4, 26), dtype=np.float64)

for i in range (26):
    for j in range (9):
        diff_low= GMI_low_final[:, j]-Sim_low_final[:,j,i]
        FG_low_final[:, j, i]= diff_low
        diff_low = diff_low[~np.isnan(diff_low)]
        diff_sqr_low = diff_low*diff_low
        Sum_val_low = np.sum(diff_sqr_low)
        n = diff_low.shape[0]
        rmse_low[j,i] = np.sqrt(Sum_val_low/n)

for i in range (26):
    for j in range (4):
        diff_high = GMI_high_final[:, j] - Sim_high_final[:, j, i]
        FG_high_final[:, j, i] = diff_high
        diff_high = diff_high[~np.isnan(diff_high)]
        diff_sqr_high = diff_high*diff_high
        Sum_val_high = np.sum(diff_sqr_high)
        n = diff_high.shape[0]
        rmse_high[j, i] = np.sqrt(Sum_val_high/n)
rmse = np.concatenate((rmse_low, rmse_high), 0)

shapes = np.linspace(1,26,26)
channels = np.linspace(1,13,1)

fig = plt.figure(figsize=(10,6))
ax1 = plt.subplot(1,1,1)
ax1.plot(shapes, rmse[12,:], 'b', label = '183+7 Ghz')
ax1.set_title("RMSE w.r.t. shapes")
ax1.set_xlabel("DDA shapes")
ax1.set_ylabel("RMSE [K]")
#ax1.set_yticks(np.arange(1, 14, step=1))
ax1.set_xlim(1,26)
#ax1.set_yticklabels(['10V', '10H', '19V', '19H', '24V', '37V', '37H', '89V', '89H', '166V', '166H', '183+3', '183+7'])
ax1.legend
plt.savefig('D:/Python_processing/cyclone/results/rmse.png', dpi=300, bbox_inches= 'tight')
plt.show()






