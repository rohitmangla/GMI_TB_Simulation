import numpy as np
from netCDF4 import Dataset
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import mpl_scatter_density
from matplotlib.colors import LogNorm
import matplotlib as mpl
import glob
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib
from scipy.stats import norm
import matplotlib.ticker as ticker
import matplotlib.ticker as mtick

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

PATH_06hrs = 'D:/Python_processing/cyclone/forecast_06hrs/'
PATH_12hrs = 'D:/Python_processing/cyclone/forecast_12hrs/'
PATH_18hrs = 'D:/Python_processing/cyclone/forecast_18hrs/'
PATH_24hrs = 'D:/Python_processing/cyclone/forecast_24hrs/'
PATH_30hrs = 'D:/Python_processing/cyclone/forecast_30hrs/'

FILES_06hrs= glob.glob(PATH_06hrs+"*.nc")
FILES_12hrs= glob.glob(PATH_12hrs+"*.nc")
FILES_18hrs= glob.glob(PATH_18hrs+"*.nc")
FILES_24hrs= glob.glob(PATH_24hrs+"*.nc")
FILES_30hrs= glob.glob(PATH_30hrs+"*.nc")

forecast_type = 5
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_low_final = np.zeros((1, channels_low), dtype=np.float64)
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_low_final_06hrs= np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final_06hrs = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_low_final_12hrs= np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final_12hrs = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_low_final_18hrs= np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final_18hrs = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_low_final_24hrs= np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final_24hrs = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_low_final_30hrs= np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final_30hrs = np.zeros((1, channels, mie_tables), dtype=np.float64)

for ifile in FILES_06hrs:
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
    Simulated_Tb_low[ind_land, :, :] = np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final_06hrs = np.concatenate((Sim_low_final_06hrs, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final_06hrs = np.concatenate((Sim_high_final_06hrs, Sim_high), 0)

Sim_low_final_06hrs[0, :, :]  = np.nan
Sim_high_final_06hrs[0, :, :] = np.nan

for ifile in FILES_12hrs:
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
    Simulated_Tb_low[ind_land, :, :] = np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final_12hrs = np.concatenate((Sim_low_final_12hrs, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final_12hrs = np.concatenate((Sim_high_final_12hrs, Sim_high), 0)

Sim_low_final_12hrs[0, :, :]  = np.nan
Sim_high_final_12hrs[0, :, :] = np.nan

for ifile in FILES_18hrs:
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
    Simulated_Tb_low[ind_land, :, :] = np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final_18hrs = np.concatenate((Sim_low_final_18hrs, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final_18hrs = np.concatenate((Sim_high_final_18hrs, Sim_high), 0)

Sim_low_final_18hrs[0, :, :]  = np.nan
Sim_high_final_18hrs[0, :, :] = np.nan

for ifile in FILES_24hrs:
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
    Simulated_Tb_low[ind_land, :, :] = np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final_24hrs = np.concatenate((Sim_low_final_24hrs, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final_24hrs = np.concatenate((Sim_high_final_24hrs, Sim_high), 0)

GMI_low_final[0, :]  = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final_24hrs[0, :, :]  = np.nan
Sim_high_final_24hrs[0, :, :] = np.nan

for ifile in FILES_30hrs:
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
    Simulated_Tb_low[ind_land, :, :] = np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final_30hrs = np.concatenate((Sim_low_final_30hrs, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final_30hrs = np.concatenate((Sim_high_final_30hrs, Sim_high), 0)

GMI_low_final[0, :]  = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final_30hrs[0, :, :]  = np.nan
Sim_high_final_30hrs[0, :, :] = np.nan

x_06= Sim_high_final_06hrs[:,12,4]
x_06 = x_06[~np.isnan(x_06)]

x_12= Sim_high_final_12hrs[:,12,4]
x_12 = x_12[~np.isnan(x_12)]

x_18= Sim_high_final_18hrs[:,12,4]
x_18 = x_18[~np.isnan(x_18)]

x_24= Sim_high_final_24hrs[:,12,4]
x_24 = x_24[~np.isnan(x_24)]

x_30= Sim_high_final_30hrs[:,12,4]
x_30 = x_30[~np.isnan(x_30)]

data_to_plot = [x_06, x_12, x_18, x_24, x_30]

# box plot

fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot(data_to_plot)

plt.savefig('D:/Python_processing/cyclone/results/box_plot.png', dpi=300, bbox_inches='tight')
