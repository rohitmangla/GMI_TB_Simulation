import numpy as np
from  mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 25
matplotlib.rcParams['legend.fontsize'] =18
matplotlib.rcParams['font.size']= 18
matplotlib.rcParams['axes.labelsize']= 17
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.2
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

lat_hudhud= [11.7, 13.2,  13.9, 15.5,  16.1]
lon_hudhud= [94.8, 90.2,  88.8,  86.4,  85.1]
lat_vardah= [9.5, 11.5, 11.7, 12.7, 13.2]
lon_vardah= [90.5, 90.5, 90.5, 88.0, 81.2]
lat_kyant= [13.7, 15.2, 16.0, 17.0, 16.7,  15.6,  15.4]
lon_kyant= [89.3, 92.5, 93.2, 91.2, 89.8, 85.0, 83.0]

# Create Plot
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 78.00, llcrnrlat= 8.00, urcrnrlon= 93.00, urcrnrlat= 25.00)
cmap= plt.cm.jet
fig = plt.figure(figsize=(15, 6))
# 37 V
ax1 = plt.subplot(1,1,1)
#ax1.set_title("(a) 37 V")
m= Basemap(resolution= 'c', projection= 'cyl',llcrnrlon= 70.00, llcrnrlat= 0.00, urcrnrlon= 99.00, urcrnrlat= 32.00)
#m.scatter(lon_gmi_low.flatten(), lat_gmi_low.flatten(), c=Tb_gmi_low[:,:, 5].flatten(), s= 5.0, cmap= cmap, edgecolors='face', linewidth=0)
cs1 = m.scatter(lon_hudhud, lat_hudhud, c= 'k', s=10.0)
cs2 = m.scatter(lon_vardah, lat_vardah, c= 'r', s=10.0)
cs3 = m.scatter(lon_kyant, lat_kyant, c= 'b', s=10.0)

m.drawparallels(np.arange(-90., 90., 5.), labels= [1, 0, 0, 0], linewidth=0)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0,0,0,1], linewidth=0)
m.drawmapboundary(fill_color='w')
m.drawcoastlines()
plt.show()
plt.savefig('D:/Python_processing/cyclone/results/cyclone_track.png', dpi=300,bbox_inches='tight')
