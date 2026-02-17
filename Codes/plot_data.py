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

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")

for ifile in FILES:
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
    #lat_min = 1.1863
    #lon_min = 75.3898
    delta= 0.1352
    lat_new = np.zeros((174, 209), dtype = np.float64)
    lon_new = np.zeros((174, 209), dtype = np.float64)
    Tb_high = np.zeros((174, 209), dtype = np.float64)
    Tb_low  = np.zeros((174, 209), dtype = np.float64)
    Sim_high = np.zeros((174, 209), dtype = np.float64)
    Sim_low  = np.zeros((174, 209), dtype = np.float64)

    for i in range(174):  # lat
        for j in range(209): # lon
            ind_data = np.where((lat>=lat_min)*(lat<lat_min+delta)*(lon>=lon_min)*(lon<lon_min+delta))
            a= lat[ind_data]
            b= lon[ind_data]
            c = GMI_Tb_high[ind_data, 0]
            d = GMI_Tb_low[ind_data, 7]
            e = Simulated_Tb_high[ind_data, 12, 4]
            f = Simulated_Tb_low[ind_data, 7, 4]

            lat_new[i,j] = np.mean(a)
            lon_new[i,j] = np.mean(b)
            Tb_high[i,j] = np.nanmean(c)
            Tb_low[i,j]  = np.nanmean(d)
            Sim_high[i, j] = np.nanmean(e)
            Sim_low[i, j] = np.nanmean(f)

            lon_min = lon_min+delta
        lat_min = lat_min+delta
        lon_min = 75.3898


    ind_high = np.where(Tb_high>1)
    GMI_high = np.squeeze(GMI_Tb_high[ind_high, :])
    #lat_low = lat[ind_high]
    #lon_low = lon[ind_high]

    #Sim_high = np.squeeze(Simulated_Tb_high[ind_high, :, :])    # (profiles, channels, mietables)
    #FG_high   = GMI_high-Sim_high[:,9:13 ,4]

    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    cmap= plt.cm.coolwarm
    fig = plt.figure(figsize=(20, 13))
    # 183+7 V
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_title("(a) GMI Radiances")
    m.scatter(lon_new.flatten(), lat_new.flatten(), c=Tb_high.flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([150, 300])

    ax2 = plt.subplot(2, 2, 2)
    ax2.set_title(" (b) Simulated Radiances")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    m.scatter(lon_new.flatten(), lat_new.flatten(), c=Tb_low.flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([150, 300])

    ax3 = plt.subplot(2, 2, 3)
    ax3.set_title("(a) GMI Radiances")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    m.scatter(lon_new.flatten(), lat_new.flatten(), c=Sim_high.flatten(), s=5.0, cmap=cmap, edgecolors='face',
              linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([150, 300])

    ax4 = plt.subplot(2, 2, 4)
    ax4.set_title(" (b) Simulated Radiances")
    m = Basemap(resolution='c', projection='cyl', llcrnrlon=75.00, llcrnrlat=5.00, urcrnrlon=100.00, urcrnrlat=20.00)
    m.scatter(lon_new.flatten(), lat_new.flatten(), c=Sim_low.flatten(), s=5.0, cmap=cmap, edgecolors='face',
              linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([150, 300])

    plt.show()
    error

    ax3 = plt.subplot(2, 2, 3)
    ax3.set_title(" (c) First Guess Departures")
    m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 75.00, llcrnrlat= 5.00, urcrnrlon= 100.00, urcrnrlat= 20.00)
    m.scatter(lon_low, lat_low, c=FG_high[:, 3], s= 5.0, cmap= cmap,edgecolors='face', linewidth=0)
    m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    m.colorbar(location='right')
    plt.clim([-50, 50])

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
    plt.savefig('D:/Python_processing/cyclone/results/Spatial_distribution/final_thinplate'+ str(i)+'.png', dpi=300, bbox_inches='tight')
    i=i+1
    error


