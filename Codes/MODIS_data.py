import numpy as np
from pyhdf.SD import SD, SDC
import os
os.environ['PROJ_LIB'] = r'C:\Users\Rohit\.conda\pkgs\proj4-5.2.0-ha925a31_1\Library\share'
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
from matplotlib.colors import BoundaryNorm

PATH_MODIS = "D:/Datasets/MODIS_L3/daily/"
MODIS_FILE = glob.glob(PATH_MODIS+"*.hdf")

AOD_daily = np.zeros((180, 360,4), dtype = np.float64)
k=0
for ifile in MODIS_FILE:
    print (ifile)
    file= SD (ifile, SDC.READ)
    sds_long  = file.select('XDim')
    Lon       = sds_long.get()
    sds_lat   = file.select('YDim')
    Lat       = sds_lat.get()
    sds_AOD   = file.select ('AOD_550_Dark_Target_Deep_Blue_Combined_Mean')
    AOD_550   = sds_AOD.get()
    AOD_550   = AOD_550 * 0.001
    ind = np.where(AOD_550<0)
    AOD_550[ind] = np.nan
    cmap = plt.cm.jet
    AOD_min = 0
    AOD_max = 1
    res_colorbar = 0.1
    clevs = np.arange(0, (AOD_max-AOD_min)/res_colorbar+1)*res_colorbar+AOD_min
    m= Basemap(resolution = 'c', projection= 'cyl', llcrnrlon = 64., llcrnrlat = 5.0, urcrnrlon =96.0, urcrnrlat= 40.0)
    #m= Basemap(resolution = 'c', projection= 'cyl', llcrnrlon = 73.6753792663, llcrnrlat = 18.197700914, urcrnrlon =135.026311477, urcrnrlat= 53.4588044297)
    xx, yy = np.meshgrid(Lon, Lat)
    fig= plt.figure()
    norm = BoundaryNorm(clevs, ncolors=cmap.N, clip=True)
    #cs = m.contourf(xx, yy, AOD_550, 5, levels=clevs, cmap=cmap)
    cs = m.pcolor(xx, yy, AOD_550, cmap= cmap, norm= norm)
    m.drawcoastlines()
    m.drawmapboundary()
    parallels = np.arange(0., 90, 8.)
    m.drawparallels(parallels, labels= [1,0,0,0], fontsize= 16, linewidth= 0.0)
    meridians = np.arange(64., 104., 8.)
    m.drawmeridians(meridians, labels= [0,0,0,1], fontsize= 16, linewidth= 0.0)
    cbar = m.colorbar(cs, location= 'right')
    cbar.set_label('Aerosol optical depth')
    plt.savefig('D:/Datasets/MODIS_L3/daily/India_'+str(k)+'.png', dpi= 300)
    k=k+1






