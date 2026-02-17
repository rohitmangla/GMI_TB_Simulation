import numpy as np
import h5py
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
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

PATH = 'D:/Python_processing/cyclone/GMI/'
GMI_file = '1B.GPM.GMI.TB2016.20141009-S045030-E062300.003478.V05A.HDF5'
filename = PATH+GMI_file

with h5py.File(filename, 'r') as f:
    # read dataset
    name_low = '/S1/Tb'
    name_high = '/S2/Tb'
    Tb_gmi_low = f[name_low][:, :, :]
    Tb_gmi_high = f[name_high][:, :, :]
    units = f[name_low].attrs['units']
    FillValue = f[name_low].attrs['_FillValue']
    Tb_gmi_low = np.ma.masked_where(np.isnan(Tb_gmi_low), Tb_gmi_low)
    Tb_gmi_high = np.ma.masked_where(np.isnan(Tb_gmi_high), Tb_gmi_high)
    # geo located data
    lat_gmi_low = f['/S1/Latitude'][:]
    lon_gmi_low = f['/S1/Longitude'][:]
    azimuth_low = f['/S1/satAzimuthAngle'][:]

    lat_gmi_high = f['/S2/Latitude'][:]
    lon_gmi_high = f['/S2/Longitude'][:]
    azimuth_high = f['/S2/satAzimuthAngle'][:]

# Create Plot
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
cmap= plt.cm.jet
fig = plt.figure(figsize=(35, 18))
# 19 V
ax1 = plt.subplot(2, 4, 1)
ax1.set_title("(a) 19 V")
m.scatter(lon_gmi_low.flatten(), lat_gmi_low.flatten(), c=Tb_gmi_low[:,:, 2].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar = m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

# 23 V
ax2 = plt.subplot(2, 4, 2)
ax2.set_title("(b) 23 V")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
m.scatter(lon_gmi_low.flatten(), lat_gmi_low.flatten(), c=Tb_gmi_low[:,:, 4].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar =m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

# 37 V
ax3 = plt.subplot(2, 4, 3)
ax3.set_title("(c) 37 V")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
m.scatter(lon_gmi_low.flatten(), lat_gmi_low.flatten(), c=Tb_gmi_low[:,:, 5].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar =m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

# 89 V
ax4 = plt.subplot(2, 4, 4)
ax4.set_title("(d) 89 V")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
m.scatter(lon_gmi_low.flatten(), lat_gmi_low.flatten(), c=Tb_gmi_low[:,:, 7].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar = m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

# 166 V
ax5 = plt.subplot(2, 4, 5)
ax5.set_title("(e) 166 V")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
m.scatter(lon_gmi_high.flatten(), lat_gmi_high.flatten(), c=Tb_gmi_high[:,:, 0].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar = m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

# 183+3 V
ax6 = plt.subplot(2, 4, 6)
ax6.set_title('(f) $183\pm\ 3V$')
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
m.scatter(lon_gmi_high.flatten(), lat_gmi_high.flatten(), c=Tb_gmi_high[:,:, 2].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar =m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

# 183+7 V
ax7 = plt.subplot(2, 4, 7)
ax7.set_title('(f) $183\pm\ 7V$')
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
m.scatter(lon_gmi_high.flatten(), lat_gmi_high.flatten(), c=Tb_gmi_high[:,:, 3].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
cbar =m.colorbar(location='right')
cbar.set_label('Tb (K)')
plt.clim([150, 300])

plt.savefig('D:/Python_processing/cyclone/results/multi_freq_dispaly_2.png', dpi=300,bbox_inches='tight')
