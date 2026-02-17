import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import glob
matplotlib.rcParams['axes.titlesize'] = 16
matplotlib.rcParams['legend.fontsize'] =11
matplotlib.rcParams['font.size']= 25
matplotlib.rcParams['axes.labelsize']= 19
matplotlib.rcParams['xtick.labelsize'] = 19
matplotlib.rcParams['ytick.labelsize'] = 19
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.2
matplotlib.rcParams['figure.subplot.wspace'] = 0.2

PATH = 'D:/Python_processing/cyclone/First_Guess_2/'
FILES = glob.glob(PATH+"*.nc")


channels = np.linspace(1,7, num=7)

SD_19V = [12.832, 10.301, 10.508, 11.630, 12.771]
SD_23V = [5.970, 4.979, 5.208, 5.584, 6.044]
SD_37V = [11.579,9.956, 10.152, 10.626, 11.588]
SD_89V = [18.345, 11.802, 14.263, 17.641, 15.053,]
SD_166V = [29.990, 18.077, 19.771, 25.139, 26.028]
SD_183_3V = [14.319, 9.400, 11.580, 13.925, 13.422]
SD_183_7V = [23.483, 13.872, 16.186, 20.570, 21.111]
Mean_19V = [4.385, 2.876, 2.577, 1.981, 4.948]
Mean_23V = [1.209, 0.273, 0.211, 0.040, 1.557]
Mean_37V = [5.553, 3.598, 3.724, 3.193, 6.217]
Mean_89V = [1.135, -0.064, 1.148, 2.451, -1.227]
Mean_166V = [0.099, -1.079, 0.422, 2.845, -4.907]
Mean_183_3V = [0.285, -0.261, 0.535, 0.627, -1.544]
Mean_183_7V = [-0.020, -1.230, 0.214, 1.586, -3.520]

fig = plt.figure(figsize=(14,5))
ax0 = plt.subplot(1,2,1)
# for 166V
x_val = np.linspace(1, 5, num=5)
ax0.plot(x_val, SD_19V, color='r', linestyle= '--',  label ='19V')
ax0.plot(x_val, SD_23V, color='r', linestyle= '-.',  label ='23V')
ax0.plot(x_val, SD_37V, color='r', linestyle= ':',  label ='37V')
ax0.plot(x_val, SD_89V, color='r', linestyle= '-',  label ='89V')
ax0.plot(x_val, SD_166V, 'b--', label ='166V')
ax0.plot(x_val, SD_183_3V, 'b:', label ='$183\pm\ 3V$')
ax0.plot(x_val, SD_183_7V, 'b', label ='$183\pm\ 7V$')

ax0.set_title('(a)')

ax0.set_xlabel(" Forecast time")
ax0.set_ylabel("Std. dev. (Obs-FG) [K]")
ax0.set_ylim(8, 32)
ax0.set_xticks(np.arange(1, 6, step=1))
ax0.set_xticklabels(['T+6','T+12', 'T+18','T+24','T+30'])
ax0.legend()

ax1 = plt.subplot(1,2,2)
ax1.set_title('(b)')
ax1.plot(x_val, Mean_19V, color='r', linestyle= '--',  label ='19V')
ax1.plot(x_val, Mean_23V, color='r', linestyle= '-.',  label ='23V')
ax1.plot(x_val, Mean_37V, color='r', linestyle= ':',  label ='37V')
ax1.plot(x_val, Mean_89V, color='r', linestyle= '-',  label ='89V')
ax1.plot(x_val, Mean_166V, 'b--', label ='166V')
ax1.plot(x_val, Mean_183_3V, 'b:', label ='$183\pm\ 3V$')
ax1.plot(x_val, Mean_183_7V, 'b', label ='$183\pm\ 7V$')
ax1.plot(x_val, np.linspace(0, 0, num=5), 'k--')

ax1.set_xlabel(" Forecast time")
ax1.set_ylabel("Mean Bias(Obs-FG) [K]")
ax1.set_ylim(-5,5)
ax1.set_xticks(np.arange(1, 6, step=1))
ax1.set_xticklabels(['T+6','T+12', 'T+18','T+24','T+30'])

plt.savefig('D:/Python_processing/cyclone/results/mean_SD_T.png', dpi=300, bbox_inches= 'tight')
plt.show()

# Distribution of FG departures
