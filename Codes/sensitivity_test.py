import numpy as np
from netCDF4 import Dataset
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 22
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.2
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/forecast_24hrs/'
FILES= glob.glob(PATH+"*.nc")
channels_high =4
channels_low =9

ifile = FILES[1]
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

ind_high = np.where(GMI_Tb_high[:,0]>1)
ind_low = np.where(GMI_Tb_low[:,0]>1)

GMI_high = np.squeeze(GMI_Tb_high[ind_high, :])
GMI_low  = np.squeeze(GMI_Tb_low[ind_low, :])

Sim_low = np.squeeze(Simulated_Tb_low[ind_low, :, :])    # (profiles, channels, mietables)
Sim_high = np.squeeze(Simulated_Tb_high[ind_high, :, :])    # (profiles, channels, mietables)

FG_high   = GMI_high-Sim_high[:,9:13 ,4]
FG_low    = GMI_low-Sim_low[:,0:9 ,4]

Mean_high = np.zeros(channels_high, dtype = np.float64)
SD_high   = np.zeros(channels_high, dtype = np.float64)
Mean_low  = np.zeros(channels_low, dtype = np.float64)
SD_low    = np.zeros(channels_low, dtype = np.float64)
RMSE_low =  np.zeros(channels_low, dtype = np.float64)
RMSE_high =  np.zeros(channels_high, dtype = np.float64)

for ichannel in range(channels_high):
    FG = FG_high[:, ichannel]
    FG = FG[~np.isnan(FG)]
    FG_sq = FG**2
    n = np.size(FG_sq)
    FG_sq_sum = (np.sum(FG_sq))/n
    RMSE_high[ichannel] = np.sqrt(FG_sq_sum)
    Mean_high[ichannel] = np.nanmean(FG)
    SD_high[ichannel]   = np.nanstd(FG)
for ichannel in range(channels_low):
    FG = FG_low[:, ichannel]
    FG = FG[~np.isnan(FG)]
    FG_sq = FG**2
    n = np.size(FG_sq)
    FG_sq_sum = (np.sum(FG_sq) )/n
    RMSE_low[ichannel] = np.sqrt(FG_sq_sum)
    Mean_low[ichannel] = np.nanmean(FG)
    SD_low[ichannel]   = np.nanstd(FG)



#plt.savefig('D:/Python_processing/cyclone/results/Spatial_distribution/final_thinplate'+ str(i)+'.png', dpi=300, bbox_inches='tight')


