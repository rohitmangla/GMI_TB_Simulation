# define map function
import numpy as np
import pylab as plt
import matplotlib.pyplot as plt
def make_plt(m, lat, lon, data, TITLESTRING):
    x, y = m(lon, lat)
    cs = m.scatter(x,y, c=data, s=2.0, cmap=plt.cm.coolwarm)
    m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 0, 0], fontsize=16, linewidth=0)
    m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 0, 0, 1], fontsize=16, linewidth=0)
    m.drawmapboundary(fill_color='w')
    m.drawcoastlines()
    #cbar = m.colorbar()
    plt.title(' %s' % (TITLESTRING), fontsize=13, weight='normal', \
              style='normal', stretch='normal', family='Times New Roman')



