import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import matplotlib
from scipy.stats import norm

matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['legend.fontsize'] =8
matplotlib.rcParams['font.size']= 14
matplotlib.rcParams['axes.labelsize']= 12
matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.25
matplotlib.rcParams['figure.subplot.wspace'] = 0.3


PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
FILES= glob.glob(PATH+"*.nc")
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_low_final = np.zeros((1, channels_low), dtype=np.float64)
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_low_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)

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
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)

GMI_low_final[0, :] = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
Sim_high_final[0, :, :] = np.nan
FG_high_final = GMI_high_final-Sim_high_final[:, 9:13, 4]
FG_low_final  = GMI_low_final-Sim_low_final[:, 0:9, 4]

# Compute a histogram of the sample
bins = np.linspace(-100, 100, 50)
# channel 10
fig = plt.figure(figsize=(5,5))
plt.yscale("log")
histogram_19V, bins_19V = np.histogram(FG_low_final[:, 2], bins=bins, normed=True)
histogram_23V, bins_23V = np.histogram(FG_low_final[:, 4], bins=bins, normed=True)
histogram_37V, bins_37V = np.histogram(FG_low_final[:, 5], bins=bins, normed=True)
histogram_89V, bins_89V = np.histogram(FG_low_final[:, 7], bins=bins, normed=True)
histogram_166V, bins_166V = np.histogram(FG_high_final[:, 0], bins=bins, normed=True)
histogram_183_3V, bins_183_3V = np.histogram(FG_high_final[:, 2], bins=bins, normed=True)
histogram_183_7V, bins_183_7V = np.histogram(FG_high_final[:, 3], bins=bins, normed=True)

bin_centers_19V = 0.5*(bins_19V[1:] + bins_19V[:-1])
bin_centers_23V = 0.5*(bins_23V[1:] + bins_23V[:-1])
bin_centers_37V = 0.5*(bins_37V[1:] + bins_37V[:-1])
bin_centers_89V = 0.5*(bins_89V[1:] + bins_89V[:-1])
bin_centers_166V = 0.5*(bins_166V[1:] + bins_166V[:-1])
bin_centers_183_3V = 0.5*(bins_183_3V[1:] + bins_183_3V[:-1])
bin_centers_183_7V = 0.5*(bins_183_7V[1:] + bins_183_7V[:-1])

# Compute the PDF on the bin centers from scipy distribution object
pdf_19V = norm.pdf(bin_centers_19V)
pdf_23V = norm.pdf(bin_centers_23V)
pdf_37V = norm.pdf(bin_centers_37V)
pdf_89V = norm.pdf(bin_centers_89V)
pdf_166V = norm.pdf(bin_centers_166V)
pdf_183_3V = norm.pdf(bin_centers_183_3V)
pdf_183_7V = norm.pdf(bin_centers_183_7V)

plt.plot(bin_centers_19V, histogram_19V, 'k--', label= '19V')
plt.plot(bin_centers_23V, histogram_23V, 'k-.', label = '23V')
plt.plot(bin_centers_37V, histogram_37V, 'k:', label= '37V')
plt.plot(bin_centers_89V, histogram_89V, 'k*', label= '89V')
plt.plot(bin_centers_166V, histogram_166V, 'kx', label = '166V')
plt.plot(bin_centers_183_3V, histogram_183_3V, 'k+', label= '$183\pm\ 3V$')
plt.plot(bin_centers_183_7V, histogram_183_7V, 'k', label= '$183\pm\ 7V$')

plt.xlim([-100, 100])
plt.ylim([1e-4, 1e0])
plt.xlabel("FG Departures")
plt.ylabel("PDF")
plt.legend(loc='upper right')
plt.savefig('D:/Python_processing/cyclone/results/FG_Departures.png', dpi=300, bbox_inches= 'tight')
plt.show()


