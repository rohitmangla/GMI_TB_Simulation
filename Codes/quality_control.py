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

# import datasets
PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
PATH_Clr = 'D:/Python_processing/cyclone/forecast_clear_sky/'
FILES= glob.glob(PATH+"*.nc")
FILES_Clr = glob.glob(PATH_Clr + "*.nc")
channels_low = 9
channels_high = 4
channels = 13
mie_tables = 26
ifile = FILES[1]
print(ifile)
ncfile = Dataset(ifile,'r')
lat = ncfile.variables['lat'][:]    # latitude
lon = ncfile.variables['lon'][:]    # longitude

Simulated_Tb_low = ncfile.variables['Simulated_Tb_low'][:]  # low frequency Simulation
Simulated_Tb_high = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
GMI_Tb_low = ncfile.variables['GMI_gridded_low'][:]  # low frequency GMI Tb
GMI_Tb_high = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
lsm = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
# Land masking
ind_land = np.where(lsm == 1)
Simulated_Tb_low[ind_land, :, :] = np.nan
Simulated_Tb_high[ind_land, :, :] = np.nan
GMI_Tb_low[ind_land, :] = np.nan
GMI_Tb_high[ind_land, :] = np.nan
ind = np.where(GMI_Tb_low[:, 0] > 0)
lat = lat[ind]
lon = lon[ind]
GMI_low = np.squeeze(GMI_Tb_low[ind, :])  # (profiles, channels_low)
GMI_high = np.squeeze(GMI_Tb_high[ind, :])  # (profiles, channels_high)
Sim_low = np.squeeze(Simulated_Tb_low[ind, :, :])  # (profiles, channels, mietables)
Sim_high = np.squeeze(Simulated_Tb_high[ind, :, :])  # (profiles, channels, mietables)

FG_low = GMI_low - Sim_low[:, 0:9, 4]
FG_high = GMI_high - Sim_high[:, 9:13, 4]

# import Clear sky radiances
ifile_clr =FILES_Clr[1]
print(ifile_clr)
ncfile = Dataset(ifile_clr, 'r')
lat = ncfile.variables['lat'][:]  # latitude
lon = ncfile.variables['lon'][:]  # longitude
Simulated_Tb_low_clr = ncfile.variables['Simulated_Tb_low'][:]  # low frequency Simulation
Simulated_Tb_high_clr = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
GMI_Tb_low = ncfile.variables['GMI_gridded_low'][:]  # low frequency GMI Tb
GMI_Tb_high = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
lsm = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
# Land masking
ind_land = np.where(lsm == 1)
Simulated_Tb_low_clr[ind_land, :] = np.nan
Simulated_Tb_high_clr[ind_land, :] = np.nan
GMI_Tb_low[ind_land, :] = np.nan
GMI_Tb_high[ind_land, :] = np.nan
ind = np.where(GMI_Tb_low[:, 0] > 0)
lat = lat[ind]
lon = lon[ind]
Sim_low_clr = np.squeeze(Simulated_Tb_low_clr[ind, :])  # (profiles, channels)
Sim_high_clr = np.squeeze(Simulated_Tb_high_clr[ind, :])  # (profiles, channels)

# polarization differences at 183 Ghz

PD_sim = (Sim_low[:, 5, 4] - Sim_low[:, 6, 4]) / (Sim_low_clr[:, 5] - Sim_low_clr[:, 6])
CA_sim = 1 - PD_sim

PD_obs = (GMI_low[:, 5] - GMI_low[:, 6]) / (Sim_low_clr[:, 5] - Sim_low_clr[:, 6])
CA_obs = 1 - PD_obs

CA_avg = (CA_obs + CA_sim) / 2

cloud_amount_std = np.zeros((22, 2, 13), dtype=np.float64)

for i in range(9):
    FG = GMI_low[:, i] - Sim_low[:, i, 4]
    cmin = -0.05
    for ibin in range(22):
        ind = np.where((CA_avg > cmin) * (CA_avg <= cmin + 0.05))
        FG_sub = FG[ind]
        cloud_amount_std[ibin, 0, i] = cmin + 0.05
        cloud_amount_std[ibin, 1, i] = np.nanstd(FG_sub)
        cmin = cmin + 0.05

for i in range(4):
    FG = GMI_high[:, i] - Sim_high[:, i + 9, 4]
    cmin = -0.05
    for ibin in range(22):
        ind = np.where((CA_avg >= cmin) * (CA_avg < cmin + 0.05))
        FG_sub = FG[ind]
        cloud_amount_std[ibin, 0, i + 9] = cmin + 0.05
        cloud_amount_std[ibin, 1, i + 9] = np.nanstd(FG_sub)
        cmin = cmin + 0.05

# Normalized PDF
# for 166 V
clr = 0.0
cld_166V = 0.5
gclr_166V = 10.29
gcld_166V = 43.09

ind_166V_1 = np.where(CA_avg < clr)
FG_166V_1 = np.squeeze(FG_high[ind_166V_1, 0])
lat_1     = lat[ind_166V_1]
lon_1     = lon[ind_166V_1]

#FG_166V_1 = FG_166V_1[~np.isnan(FG_166V_1)]
FG_166V_norm_1 = FG_166V_1/gclr_166V

ind_166V_2 = np.where(CA_avg > cld_166V)
FG_166V_2 = np.squeeze(FG_high[ind_166V_2, 0])
lat_2     = lat[ind_166V_2]
lon_2     = lon[ind_166V_2]

#FG_166V_2 = FG_166V_2[~np.isnan(FG_166V_2)]
FG_166V_norm_2 = FG_166V_2/gcld_166V

ind_166V_3 = np.where((CA_avg >= clr) * (CA_avg <= cld_166V))
FG_166V_3 = np.squeeze(FG_high[ind_166V_3, 0])
lat_3     = lat[ind_166V_3]
lon_3     = lon[ind_166V_3]

CA_166V_3 = CA_avg[ind_166V_3]
SD_166V = gclr_166V + ((gcld_166V - gclr_166V) / (cld_166V - clr)) * (CA_166V_3 - clr)
FG_166V_norm_3 = FG_166V_3 / SD_166V
#FG_166V_norm_3 = FG_166V_norm_3[~np.isnan(FG_166V_norm_3)]

FG_166V_norm = np.concatenate((FG_166V_norm_1, FG_166V_norm_2, FG_166V_norm_3), 0)
lat_norm = np.concatenate((lat_1, lat_2, lat_3), 0)
lon_norm = np.concatenate((lon_1, lon_2, lon_3), 0)
ind_rem = np.where((FG_166V_norm>-2.5)*(FG_166V_norm<2.5))
FG_166V_norm_2 = FG_166V_norm[ind_rem]
lat_norm_2 = lat_norm[ind_rem]
lon_norm_2 = lon_norm[ind_rem]

mu_4 = np.mean(FG_166V_norm)
sigma_4 = np.std(FG_166V_norm)
median_4 = np.median(FG_166V_norm)

m = Basemap(resolution='c', projection='cyl', llcrnrlon=75.00, llcrnrlat=5.00, urcrnrlon=100.00, urcrnrlat=20.00)
#cmap = plt.cm.coolwarm
cmap = plt.cm.coolwarm

fig = plt.figure(figsize=(14, 20))
# 183+7 V
ax1 = plt.subplot(3, 1, 1)
ax1.set_title("(a) FG departure")
m.scatter(lon, lat, c=FG_high[:,0], s=5.0, cmap=cmap, edgecolors='face',linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
m.colorbar(location='right')
plt.clim([-50, 50])

# Normalized FG departures
m = Basemap(resolution='c', projection='cyl', llcrnrlon=75.00, llcrnrlat=5.00, urcrnrlon=100.00, urcrnrlat=20.00)
ax2 = plt.subplot(3, 1, 2)
ax2.set_title("(b) Normalized FG departure-all")
m.scatter(lon_norm, lat_norm, c=FG_166V_norm, s=5.0, cmap=cmap, edgecolors='face',linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
m.colorbar(location='right')
plt.clim([-5, 5])

m = Basemap(resolution='c', projection='cyl', llcrnrlon=75.00, llcrnrlat=5.00, urcrnrlon=100.00, urcrnrlat=20.00)
ax3 = plt.subplot(3, 1, 3)
ax3.set_title("(c) Normalized FG departure-screened")
m.scatter(lon_norm_2, lat_norm_2, c=FG_166V_norm_2, s=5.0, cmap=cmap, edgecolors='face',linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
m.colorbar(location='right')
plt.clim([-5, 5])
plt.savefig('D:/Python_processing/cyclone/results/error_model/quality_control.png', dpi=300, bbox_inches='tight')




