import numpy as np
from netCDF4 import Dataset
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from make_plt import make_plt
import glob
PATH = 'D:/Python_processing/cyclone/clear_sky/'
FILES= glob.glob(PATH+"*.nc")
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
    Simulated_Tb_low[ind_land,:]= np.nan
    Simulated_Tb_high[ind_land,:]= np.nan
    GMI_Tb_low[ind_land, :]      = np.nan
    GMI_Tb_high[ind_land, :]   = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :])
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :])
    FG_low    = GMI_low-Sim_low[:, 0:9]
    FG_high   = GMI_high-Sim_high[:,9:13]

m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
cmap= plt.cm.coolwarm
fig = plt.figure()
ax = fig.add_subplot(131)
ax.set_title("GMI Radiances")
m.scatter(lon, lat, c=GMI_high[:, 3], s= 5.0, cmap= cmap,\
          edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
m.colorbar(location='right')
plt.clim([150, 300])

ax = fig.add_subplot(132)
ax.set_title("Simulated Radiances")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
m.scatter(lon, lat, c=Sim_high[:, 12], s= 5.0, cmap= cmap,\
          edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
m.colorbar(location='right')
#plt.clim([225, 300])

ax = fig.add_subplot(133)
ax.set_title("First Guess Departures")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
m.scatter(lon, lat, c=FG_high[:, 3], s= 5.0, cmap= cmap,\
          edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
m.colorbar(location='right')
plt.clim([-50, 50])
plt.show()
