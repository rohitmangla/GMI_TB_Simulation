import numpy as np
from netCDF4 import Dataset
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
import matplotlib
import ez_color

matplotlib.rcParams['axes.titlesize'] = 27
matplotlib.rcParams['legend.fontsize'] =20
matplotlib.rcParams['font.size']= 26
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 24
matplotlib.rcParams['ytick.labelsize'] = 24
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.05
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")
i=1

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
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])
    FG_low    = GMI_low-Sim_low[:, 0:9, 24]
    FG_high   = GMI_high-Sim_high[:,9:13 ,24]

    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 70.00, llcrnrlat= 0.00, urcrnrlon= 99.00, urcrnrlat= 32.00)
    cmap= plt.cm.jet
    fig = plt.figure(figsize=(14, 18))
    ax1 = plt.subplot(3, 1, 1)
    ax1.set_title("37V-Obs")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 70.00, llcrnrlat= 0.00, urcrnrlon= 99.00, urcrnrlat= 32.00)
    m.scatter(lon, lat, c=GMI_low[:, 5], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([170, 280])
    fig.tight_layout()

    ax2= plt.subplot(3, 1, 2)
    ax2.set_title("37V-Sim")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 70.00, llcrnrlat= 0.00, urcrnrlon= 99.00, urcrnrlat= 32.00)
    m.scatter(lon, lat, c=Sim_low[:, 5, 24], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([170, 280])
    fig.tight_layout()

    ax3 = plt.subplot(3,1,3)
    ax3.set_title("37V-FG Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 70.00, llcrnrlat= 0.00, urcrnrlon= 99.00, urcrnrlat= 32.00)
    m.scatter(lon, lat, c=FG_low[:, 5], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()
    plt.savefig('D:/Python_processing/cyclone/results/Spatial_distribution/forward_sim/final_37'+ str(i)+'.png', dpi=300, bbox_inches='tight')
    i=i+1
    plt.close()

