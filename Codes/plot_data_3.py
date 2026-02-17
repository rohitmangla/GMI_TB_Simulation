import numpy as np
from netCDF4 import Dataset
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from make_plt import make_plt
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
    ind = np.where(GMI_Tb_high[:, 0]>0)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])
    FG_low    = GMI_low-Sim_low[:, 0:9, 8]
    FG_high   = GMI_high-Sim_high[:,9:13 ,8]

    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    cmap= plt.cm.coolwarm
    fig = plt.figure(figsize=(30, 13))
    # 183+7 V
    ax1 = plt.subplot(2, 4, 1)
    ax1.set_title("(a) 37V-Obs")
    m.scatter(lon, lat, c=GMI_low[:, 4], s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([150, 300])

    ax2 = plt.subplot(2, 4, 2)
    ax2.set_title(" 37V-Sim ")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=Sim_low[:,4,9], s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([150, 300])

    ax3 = plt.subplot(2, 4, 3)
    ax3.set_title(" (c) 37V-FG dep")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    m.scatter(lon, lat, c=FG_low[:, 4], s= 5.0, cmap= cmap,edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])

    ax4 = plt.subplot(2,4, 4)
    ax4.set_title(" (d) Histogram of FG departures")
    FG_low_37V = FG_low[:, 4]
    mu = np.nanmean(FG_low_37V)
    std = np.nanstd(FG_low_37V)
    bin = np.linspace(-50, 50)
    # bin = [-30, -20, -15, -10, -5, 0, 5, 10, 15, 20, 30, 35, 40, 45]
    ax4.hist(FG_low_37V, bins=bin, facecolor='gold')
    ax4.text(-10, 5000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax4.set_xlabel('FG departures')
    ax4.set_xlim([-50, 50])
    ax4.set_ylabel('frequency')
    plt.show()
    error

    ax4 = plt.subplot(2, 2, 4)
    ax4.set_title(" (d) Histogram of FG departures")
    FG_hist_183_7V = FG_high[:, 3]
    mu =  np.nanmean(FG_hist_183_7V)
    std = np.nanstd(FG_hist_183_7V)
    bin = np.linspace(-50, 50)
    #bin = [-30, -20, -15, -10, -5, 0, 5, 10, 15, 20, 30, 35, 40, 45]
    ax4.hist(FG_hist_183_7V, bins=bin, facecolor='gold')
    ax4.text(-10, 5000, 'mu,sigma(%s, %s)' % (mu, std), fontsize=18)
    ax4.set_xlabel('FG departures')
    ax4.set_xlim([-50, 50])
    ax4.set_ylabel('frequency')

    plt.savefig('D:/Python_processing/cyclone/results/Spatial_distribution/new_183/sector/final_sector'+ str(i)+'.png', dpi=300, bbox_inches='tight')
    i=i+1

