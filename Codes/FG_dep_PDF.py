import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import matplotlib
from scipy.stats import norm

matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 22
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
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
# channle 10
fig = plt.figure(figsize=(15, 20))
ax1= plt.subplot(3,2,1)
ax1.set_title(" 19 V")
ax1.set(yscale="log")
histogram, bins = np.histogram(FG_low_final[:, 2], bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax1.plot(bin_centers, histogram, 'k')
ax1.set_xlim([-100, 100])
ax1.set_ylim([1e-4, 1e0])
ax1.set_xlabel("FG Departures")
ax1.set_ylabel("PDF")

ax2= plt.subplot(3,2,2)
ax2.set_title("23 V")
ax2.set(yscale="log")
bin_centers = 0.5*(bins[1:] + bins[:-1])
histogram, bins = np.histogram(FG_low_final[:, 4], bins=bins, normed=True)
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax2.plot(bin_centers, histogram, 'k')
ax2.set_xlim([-100, 100])
ax2.set_ylim([1e-4, 1e0])
ax2.set_xlabel("FG Departures")
ax2.set_ylabel("PDF")

ax3= plt.subplot(3,2,3)
ax3.set_title("37 V")
ax3.set(yscale="log")
bin_centers = 0.5*(bins[1:] + bins[:-1])
histogram, bins = np.histogram(FG_low_final[:, 5], bins=bins, normed=True)
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax3.plot(bin_centers, histogram, 'k')
ax3.set_xlim([-100, 100])
ax3.set_ylim([1e-4, 1e0])
ax3.set_xlabel("FG Departures")
ax3.set_ylabel("PDF")

ax4= plt.subplot(3,2,4)
ax4.set_title("89 V")
ax4.set(yscale="log")
bin_centers = 0.5*(bins[1:] + bins[:-1])
histogram, bins = np.histogram(FG_low_final[:, 7], bins=bins, normed=True)
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax4.plot(bin_centers, histogram, 'k')
ax4.set_xlim([-100, 100])
ax4.set_ylim([1e-4, 1e0])
ax4.set_xlabel("FG Departures")
ax4.set_ylabel("PDF")

ax5= plt.subplot(3,2,5)
ax5.set_title("166 V")
ax5.set(yscale="log")
bin_centers = 0.5*(bins[1:] + bins[:-1])
histogram, bins = np.histogram(FG_high_final[:, 0], bins=bins, normed=True)
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax5.plot(bin_centers, histogram, 'k')
ax5.set_xlim([-100, 100])
ax5.set_ylim([1e-4, 1e0])
ax5.set_xlabel("FG Departures")
ax5.set_ylabel("PDF")

ax6= plt.subplot(3,2,6)
ax6.set_title("183+7 V")
ax6.set(yscale="log")
bin_centers = 0.5*(bins[1:] + bins[:-1])
histogram, bins = np.histogram(FG_high_final[:, 3], bins=bins, normed=True)
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax6.plot(bin_centers, histogram, 'k')
ax6.set_xlim([-100, 100])
ax6.set_ylim([1e-4, 1e0])
ax6.set_xlabel("FG Departures")
ax6.set_ylabel("PDF")

plt.savefig('D:/Python_processing/cyclone/results/FG_Departures.png', dpi=300, bbox_inches= 'tight')
plt.show()


