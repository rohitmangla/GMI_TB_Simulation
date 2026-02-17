import numpy as np
from netCDF4 import Dataset
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
import matplotlib

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

PATH = 'D:/Python_processing/cyclone/forecast_06hrs/'
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
    FG_low    = GMI_low-Sim_low[:, 0:9, 4]
    FG_high   = GMI_high-Sim_high[:,9:13 ,4]
    n= np.where((GMI_low[:,5]>160)*(GMI_low[:,5]<200))
    GMI_low[n,5]=np.nan
    Sim_low[n,5,:]=np.nan
    FG_low[n,5]=np.nan

    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    cmap= plt.cm.jet
    fig = plt.figure(figsize=(30, 45))
    # 19 V
    ax1 = plt.subplot(7, 3, 1)
    ax1.set_title("19V-Obs")
    m.scatter(lon, lat, c=GMI_low[:, 2], s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([200, 280])
    fig.tight_layout()

    ax2 = plt.subplot(7, 3, 2)
    ax2.set_title("19V-Sim")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_low[:,2, 4], s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([200, 280])
    fig.tight_layout()

    ax3 = plt.subplot(7, 3, 3)
    ax3.set_title("19V-FG Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_low[:, 2], s= 5.0, cmap= cmap,edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()

    # 23 V
    ax5 = plt.subplot(7, 3, 4)
    ax5.set_title("23V-Obs")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=GMI_low[:, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([240, 280])
    fig.tight_layout()

    ax6 = plt.subplot(7, 3, 5)
    ax6.set_title("23V-Sim")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_low[:, 4, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([240, 280])
    fig.tight_layout()

    ax7 = plt.subplot(7, 3, 6)
    ax7.set_title("23V-FG Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_low[:, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()

    # 37V
    ax9 = plt.subplot(7, 3, 7)
    ax9.set_title("37V-Obs")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=GMI_low[:, 5], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([210, 280])
    fig.tight_layout()

    ax10= plt.subplot(7, 3, 8)
    ax10.set_title("37V-Sim")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_low[:, 5, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([210, 280])
    fig.tight_layout()

    ax11 = plt.subplot(7, 3, 9)
    ax11.set_title("37V-FG Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_low[:, 5], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()

    # 89 GHz
    ax13 = plt.subplot(7,3,10)
    ax13.set_title("89V-Obs")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=GMI_low[:, 7], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([120, 290])
    fig.tight_layout()

    ax14 = plt.subplot(7, 3, 11)
    ax14.set_title("89V-Sim")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_low[:, 7, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([120, 290])
    fig.tight_layout()

    ax15 = plt.subplot(7,3,12)
    ax15.set_title("89V-FG Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_low[:, 7], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()

    # 166V
    ax17 = plt.subplot(7, 3, 13)
    ax17.set_title("166V-Obs")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=GMI_high[:, 0], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([90, 285])
    fig.tight_layout()

    ax18 = plt.subplot(7, 3, 14)
    ax18.set_title("166V-Sim")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_high[:, 9, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([90, 285])
    fig.tight_layout()

    ax19 = plt.subplot(7, 3, 15)
    ax19.set_title("166V-FG Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_high[:, 0], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()

    # 183+3V
    ax21 = plt.subplot(7, 3, 16)
    ax21.set_title("$183\pm\ 3V-Obs$")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=GMI_high[:, 2], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([90, 265])
    fig.tight_layout()

    ax22 = plt.subplot(7, 3, 17)
    ax22.set_title("$183\pm\ 3V-Sim$")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_high[:, 11, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([90, 265])
    fig.tight_layout()

    ax23= plt.subplot(7, 3, 18)
    ax23.set_title("$183\pm\ 3V-FG Departures$")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_high[:, 2], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()


    # 183+7V
    ax25 = plt.subplot(7, 3, 19)
    ax25.set_title("$183\pm\ 7V-Obs$")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=GMI_high[:, 3], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([90, 275])
    fig.tight_layout()

    ax26= plt.subplot(7, 3, 20)
    ax26.set_title("$183\pm\ 7V-Sim$")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_high[:, 12, 4], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([90, 275])
    fig.tight_layout()

    ax27 = plt.subplot(7, 3, 21)
    ax27.set_title("$183\pm\ 7V-FG Departures$")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_high[:, 3], s=5.0, cmap=cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])
    fig.tight_layout()

    plt.savefig('D:/Python_processing/cyclone/results/Spatial_distribution/forward_sim/thinplate/6hrs/final'+ str(i)+'.png', dpi=300, bbox_inches='tight')
    i=i+1
    plt.close()

