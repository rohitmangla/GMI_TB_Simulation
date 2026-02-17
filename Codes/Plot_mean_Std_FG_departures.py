import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import glob
matplotlib.rcParams['axes.titlesize'] = 16
matplotlib.rcParams['legend.fontsize'] =17
matplotlib.rcParams['font.size']= 25
matplotlib.rcParams['axes.labelsize']= 19
matplotlib.rcParams['xtick.labelsize'] = 19
matplotlib.rcParams['ytick.labelsize'] = 19
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.2

PATH = 'D:/Python_processing/cyclone/First_Guess/'
FILES = glob.glob(PATH+"*.nc")

OmFmean_total_final = np.zeros((13, 26, 5), dtype= np.float64)
OmFstd_total_final  = np.zeros((13, 26, 5), dtype= np.float64)
OmFDep_total_final  = np.zeros((290210,13, 26, 5), dtype= np.float64)

i=0
for ifile in FILES:
    print (ifile)
    ncfile = Dataset(ifile, 'r')
    OmFmean_low = (ncfile.variables['OmFmean_low'])[:]     # (channels_low, mie_tables)
    OmFstd_low = (ncfile.variables['OmFstd_low'])[:]       # (channels_low, mie_tables)
    OmFDep_low =  (ncfile.variables['OmFDep_low'])[:]      # (nprofiles, channels_low, mie_tables)
    OmFmean_high = (ncfile.variables['OmFmean_high'])[:]   #  (channels_high, mie_tables)
    OmFstd_high = (ncfile.variables['OmFstd_high'])[:]     #  (channels_high, mie_tables)
    OmFDep_high =  (ncfile.variables['OmFDep_high'])[:]    #  (nprofiles, channels_high, mie_tables)
    OmFmean_total = np.concatenate((OmFmean_low, OmFmean_high), 0)
    OmFstd_total = np.concatenate((OmFstd_low, OmFstd_high), 0)
    OmFmean_total_final[:, :, i] = OmFmean_total
    OmFstd_total_final[:, :, i]  = OmFstd_total
    OmFDep_total_final[:, :, :, i] = np.concatenate((OmFDep_low, OmFDep_high), 1)
    i=i+1

OmFmean_total_filter = np.zeros((7, 26, 5), dtype = np.float64)
OmFstd_total_filter = np.zeros((7, 26, 5), dtype = np.float64)

# RMSE
#RMSE = np.zeros((13, 5), dtype = np.float64)
#for itime in range(5):
#    for ichannel in range(13):
#        OmFsq = (OmFDep_total_final[:,ichannel, 4, itime])**2
#        OmFsq = OmFsq[~np.isnan(OmFsq)]
#        n = np.size(OmFsq)
#        OmFsq_sum = (np.sum(OmFsq))/n
#        RMSE[ichannel, itime] = np.sqrt(OmFsq_sum)

#error


OmFmean_total_filter[0, :, :] = OmFmean_total_final[2, :, :]
OmFmean_total_filter[1, :, :] = OmFmean_total_final[4, :, :]
OmFmean_total_filter[2, :, :] = OmFmean_total_final[5, :, :]
OmFmean_total_filter[3, :, :] = OmFmean_total_final[7, :, :]
OmFmean_total_filter[4, :, :] = OmFmean_total_final[9, :, :]
OmFmean_total_filter[5, :, :] = OmFmean_total_final[11, :, :]
OmFmean_total_filter[6, :, :] = OmFmean_total_final[12, :, :]

OmFstd_total_filter[0, :, :] = OmFstd_total_final[2, :, :]
OmFstd_total_filter[1, :, :] = OmFstd_total_final[4, :, :]
OmFstd_total_filter[2, :, :] = OmFstd_total_final[5, :, :]
OmFstd_total_filter[3, :, :] = OmFstd_total_final[7, :, :]
OmFstd_total_filter[4, :, :] = OmFstd_total_final[9, :, :]
OmFstd_total_filter[5, :, :] = OmFstd_total_final[11, :, :]
OmFstd_total_filter[6, :, :] = OmFstd_total_final[12, :, :]

channels = np.linspace(1,7, num=7)
fig = plt.figure(figsize=(20,8))
ax0 = plt.subplot(1,2,1)
ax0.plot(OmFmean_total_filter[:, 4, 0], channels, 'k--', label = 'T+6')
ax0.plot(OmFmean_total_filter[:, 4, 1], channels, 'k', label = 'T+12')
ax0.plot(OmFmean_total_filter[:, 4, 2], channels, 'k-.', label = 'T+18')
ax0.plot(OmFmean_total_filter[:, 4, 3], channels, 'k:', label = 'T+24')
ax0.plot(OmFmean_total_filter[:, 4, 4], channels, 'k.', label = 'T+30')

ax0.plot(np.zeros(7), channels, color= 'silver', linestyle= '--')
ax0.set_xlabel("Mean FG departure [K]")
ax0.set_ylabel("Channel name")
ax0.set_xlim(-5,8)
ax0.set_yticks(np.arange(1, 8, step=1))
ax0.set_yticklabels(['19V','23V', '37V', '89V','166V','$183\pm\ 3V$','$183\pm\ 7V$'])
ax0.legend(loc= 'upper right')

ax1 = plt.subplot(1,2,2)
ax1.plot(OmFstd_total_filter[:, 4, 0], channels, 'k--',   label = 'T+6')
ax1.plot(OmFstd_total_filter[:, 4, 1], channels, 'k', label = 'T+12')
ax1.plot(OmFstd_total_filter[:, 4, 2], channels, 'k-.', label = 'T+18')
ax1.plot(OmFstd_total_filter[:, 4, 3], channels, 'k:', label = 'T+24')
ax1.plot(OmFstd_total_filter[:, 4, 4], channels, 'k.', label = 'T+30')

ax1.set_xlabel(" Std. dev. of FG departure [K]")
ax1.set_yticks(np.arange(1, 8, step=1))
ax1.set_xlim(5,35)
ax1.set_yticklabels(['19V','23V', '37V','89V','166V','$183\pm\ 3V$','$183\pm\ 7V$'])
#ax1.legend(loc= 'upper right')
plt.savefig('D:/Python_processing/cyclone/results/mean_std_fG_Departures.png', dpi=300, bbox_inches= 'tight')
plt.show()

# Distribution of FG departures
