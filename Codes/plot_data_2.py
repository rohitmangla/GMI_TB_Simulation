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
    FG_low    = GMI_low-Sim_low[:, 0:9, 4]
    FG_high   = GMI_high-Sim_high[:,9:13 ,4]

    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 77.00, llcrnrlat= 5.00, urcrnrlon= 97.00, urcrnrlat= 20.00)
    cmap= plt.cm.jet
    fig = plt.figure(figsize=(40, 45))
    # 19 V
    ax1 = plt.subplot(7, 4, 1)
    ax1.set_title("19V-Obs")
    m.scatter(lon, lat, c=GMI_low[:, 2], s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([200, 280])
    fig.tight_layout()

    ax2 = plt.subplot(7, 4, 2)
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

    ax3 = plt.subplot(7, 4, 3)
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

    ax4 = plt.subplot(7, 4, 4)
    ax4.set_title("Histogram 19V-FG")
    FG_hist_19V = FG_low[:, 2]
    mu = np.nanmean(FG_hist_19V)
    std = np.nanstd(FG_hist_19V)
    bin = np.linspace(-50, 50)
    ax4.hist(FG_hist_19V, bins=bin, facecolor='gold')
    ax4.text(-10, 5000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax4.set_xlabel('FG departures')
    ax4.set_xlim([-50, 50])
    ax4.set_ylabel('frequency')
    fig.tight_layout()

    # 23 V
    ax5 = plt.subplot(7, 4, 5)
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

    ax6 = plt.subplot(7, 4, 6)
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

    ax7 = plt.subplot(7, 4, 7)
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

    ax8 = plt.subplot(7, 4, 8)
    ax8.set_title("Histogram 23V-FG")
    FG_hist_23V = FG_low[:, 4]
    mu = np.nanmean(FG_hist_23V)
    std = np.nanstd(FG_hist_23V)
    bin = np.linspace(-50, 50)
    ax8.hist(FG_hist_23V, bins=bin, facecolor='gold')
    ax8.text(-10, 5000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax8.set_xlabel('FG departures')
    ax8.set_xlim([-50, 50])
    ax8.set_ylabel('frequency')
    fig.tight_layout()

    # 37V
    ax9 = plt.subplot(7, 4, 9)
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

    ax10= plt.subplot(7, 4, 10)
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

    ax11 = plt.subplot(7, 4, 11)
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

    ax12 = plt.subplot(7, 4, 12)
    ax12.set_title("Histogram 37V-FG")
    FG_hist_37V = FG_low[:, 5]
    mu = np.nanmean(FG_hist_37V)
    std = np.nanstd(FG_hist_37V)
    bin = np.linspace(-50, 50)
    ax12.hist(FG_hist_37V, bins=bin, facecolor='gold')
    ax12.text(-10, 5000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax12.set_xlabel('FG departures')
    ax12.set_xlim([-50, 50])
    ax12.set_ylabel('frequency')
    fig.tight_layout()

    # 89 GHz
    ax13 = plt.subplot(7,4,13)
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

    ax14 = plt.subplot(7, 4, 14)
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

    ax15 = plt.subplot(7,4,15)
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

    ax16 = plt.subplot(7, 4, 16)
    ax16.set_title("Histogram 89V-FG")
    FG_hist_89V = FG_low[:, 7]
    mu = np.nanmean(FG_hist_89V)
    std = np.nanstd(FG_hist_89V)
    bin = np.linspace(-50, 50)
    ax16.hist(FG_hist_89V, bins=bin, facecolor='gold')
    ax16.text(-10, 5000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax16.set_xlabel('FG departures')
    ax16.set_xlim([-50, 50])
    ax16.set_ylabel('frequency')
    fig.tight_layout()

    # 166V
    ax17 = plt.subplot(7, 4, 17)
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

    ax18 = plt.subplot(7, 4, 18)
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

    ax19 = plt.subplot(7, 4, 19)
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

    ax20 = plt.subplot(7, 4, 20)
    ax20.set_title("Histogram 166V-FG")
    FG_hist_166V = FG_high[:, 0]
    mu = np.nanmean(FG_hist_166V)
    std = np.nanstd(FG_hist_166V)
    bin = np.linspace(-50, 50)
    ax20.hist(FG_hist_166V, bins=bin, facecolor='gold')
    ax20.text(-10, 4000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax20.set_xlabel('FG departures')
    ax20.set_xlim([-50, 50])
    ax20.set_ylabel('frequency')
    fig.tight_layout()

    # 183+3V
    ax21 = plt.subplot(7, 4, 21)
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

    ax22 = plt.subplot(7, 4, 22)
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

    ax23= plt.subplot(7, 4, 23)
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

    ax24 = plt.subplot(7, 4, 24)
    ax24.set_title("$Histogram 183\pm\ 3V-FG$")
    FG_hist_183_3V = FG_high[:, 2]
    mu = np.nanmean(FG_hist_183_3V)
    std = np.nanstd(FG_hist_183_3V)
    bin = np.linspace(-50, 50)
    ax24.hist(FG_hist_183_3V, bins=bin, facecolor='gold')
    ax24.text(-10, 4000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax24.set_xlabel('FG departures')
    ax24.set_xlim([-50, 50])
    ax24.set_ylabel('frequency')
    fig.tight_layout()

    # 183+7V
    ax25 = plt.subplot(7, 4, 25)
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

    ax26= plt.subplot(7, 4, 26)
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

    ax27 = plt.subplot(7, 4, 27)
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

    ax28 = plt.subplot(7, 4, 28)
    ax28.set_title("$Histogram 183\pm\ 7V-FG$")
    FG_hist_183_7V = FG_high[:, 3]
    mu = np.nanmean(FG_hist_183_7V)
    std = np.nanstd(FG_hist_183_7V)
    bin = np.linspace(-50, 50)
    ax28.hist(FG_hist_183_7V, bins=bin, facecolor='gold')
    ax28.text(-10, 4000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax28.set_xlabel('FG departures')
    ax28.set_xlim([-50, 50])
    ax28.set_ylabel('frequency')
    fig.tight_layout()

    plt.savefig('D:/Python_processing/cyclone/results/Spatial_distribution/forward_sim/final'+ str(i)+'.png', dpi=300, bbox_inches='tight')
    i=i+1
    plt.show()
    error

