import numpy as np
from netCDF4 import Dataset
import glob

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

nprofiles = GMI_low_final.shape[0]
OmFDep_low    = np.zeros((nprofiles, channels_low, mie_tables), dtype=np.float64)
OmFmean_low   = np.zeros((channels_low, mie_tables), dtype=np.float64)
OmFstd_low    = np.zeros((channels_low, mie_tables), dtype=np.float64)

OmFDep_high    = np.zeros((nprofiles, channels_high, mie_tables), dtype=np.float64)
OmFmean_high   = np.zeros((channels_high, mie_tables), dtype=np.float64)
OmFstd_high    = np.zeros((channels_high, mie_tables), dtype=np.float64)

for ichannel in range (channels_low):
    for itable in range (mie_tables):
        OmFDep_low[:, ichannel, itable] = GMI_low_final[:, ichannel]-Sim_low_final[:, ichannel, itable]
        OmFmean_low[ichannel, itable] = np.nanmean(OmFDep_low[:, ichannel, itable], axis=0)
        OmFstd_low[ichannel, itable]  = np.nanstd(OmFDep_low[:, ichannel, itable], axis=0)

for ichannel in range (channels_high):
    for itable in range (mie_tables):
        OmFDep_high[:, ichannel, itable] = GMI_high_final[:, ichannel]-Sim_high_final[:, ichannel+9, itable]
        OmFmean_high[ichannel, itable] = np.nanmean(OmFDep_high[:, ichannel, itable], axis=0)
        OmFstd_high[ichannel, itable]  = np.nanstd(OmFDep_high[:, ichannel, itable], axis=0)

ncfile = Dataset('D:/Python_processing/cyclone/First_Guess/FG_dep_30hrs.nc', 'w')
ncfile.createDimension('channels_low', size=9)
ncfile.createDimension('channels_high', size=4)
ncfile.createDimension('mie_tables', size= 26)
ncfile.createDimension('nprofiles', size= nprofiles)

datanc = ncfile.createVariable('OmFmean_low', np.float32, ('channels_low', 'mie_tables'))
datanc[:, :]= np.float32(OmFmean_low)
datanc = ncfile.createVariable('OmFstd_low', np.float32, ('channels_low', 'mie_tables'))
datanc[:, :]= np.float32(OmFstd_low)
datanc = ncfile.createVariable('OmFDep_low', np.float32, ('nprofiles','channels_low', 'mie_tables'))
datanc[:,:,:]= np.float32(OmFDep_low)
datanc = ncfile.createVariable('OmFmean_high', np.float32, ('channels_high', 'mie_tables'))
datanc[:, :]= np.float32(OmFmean_high)
datanc = ncfile.createVariable('OmFstd_high', np.float32, ('channels_high', 'mie_tables'))
datanc[:, :]= np.float32(OmFstd_high)
datanc = ncfile.createVariable('OmFDep_high', np.float32, ('nprofiles','channels_high', 'mie_tables'))
datanc[:,:,:]= np.float32(OmFDep_high)




