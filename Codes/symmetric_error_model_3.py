import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib
from scipy.stats import norm

matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 26
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

PATH = 'D:/Python_processing/cyclone/forecast_12hrs/'
PATH_Clr = 'D:/Python_processing/cyclone/forecast_clear_sky/'
FILES= glob.glob(PATH+"*.nc")
FILES_Clr= glob.glob(PATH_Clr+"*.nc")
channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_low_final = np.zeros((1, channels_low), dtype=np.float64)
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
Sim_low_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)

Sim_low_clr_final = np.zeros((1, channels), dtype=np.float64)
Sim_high_clr_final = np.zeros((1, channels), dtype=np.float64)

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
    Simulated_Tb_low[ind_land, :, :] = np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>1)
    lat = lat[ind]
    lon = lon[ind]
    GMI_low   = np.squeeze(GMI_Tb_low[ind, :])              # (profiles, channels_low)
    GMI_low_final = np.append(GMI_low_final, GMI_low, axis= 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind, :])             # (profiles, channels_high)
    GMI_high_final = np.append(GMI_high_final, GMI_high, axis=0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind, :, :])     # (profiles, channels, mietables)
    Sim_low_final = np.append(Sim_low_final, Sim_low, axis=0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind, :, :])    # (profiles, channels, mietables)
    Sim_high_final = np.append(Sim_high_final, Sim_high, axis=0)

GMI_low_final[0, :]  = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :]  = np.nan
Sim_high_final[0, :, :] = np.nan
n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan

FG_89V = GMI_low_final[:, 7]-Sim_low_final[:, 7, 4]
mu_89V = np.nanmean(FG_89V)
std_89V = np.nanstd(FG_89V)
#FG_89V[FG_89V>2*std_89V]= np.nan
#FG_89V[FG_89V<-2*std_89V]= np.nan

FG_166V = GMI_high_final[:, 0]-Sim_high_final[:, 9, 4]
mu_166V = np.nanmean(FG_166V)
std_166V = np.nanstd(FG_166V)
#FG_166V[FG_166V>2*std_166V]= np.nan
#FG_166V[FG_166V<-2*std_166V]= np.nan

FG_183_3V = GMI_high_final[:, 2]-Sim_high_final[:, 11, 4]
mu_183_3V = np.nanmean(FG_183_3V)
std_183_3V = np.nanstd(FG_183_3V)
#FG_183_3V[FG_183_3V>2*std_183_3V]= np.nan
#FG_183_3V[FG_183_3V<-2*std_183_3V]= np.nan

FG_183_7V = GMI_high_final[:, 3]-Sim_high_final[:, 12, 4]
mu_183_7V = np.nanmean(FG_183_7V)
std_183_7V = np.nanstd(FG_183_7V)
#FG_183_7V[FG_183_7V>2*std_183_7V]= np.nan
#FG_183_7V[FG_183_7V<-2*std_183_7V]= np.nan

# import Clear sky radiances
for ifile_clr in FILES_Clr:
    print(ifile_clr)
    ncfile = Dataset(ifile_clr,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_low_clr  = ncfile.variables['Simulated_Tb_low'][:]   # low frequency Simulation
    Simulated_Tb_high_clr = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    GMI_Tb_low        = ncfile.variables['GMI_gridded_low'][:]   # low frequency GMI Tb
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
    lsm               = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
    # Land masking
    ind_land = np.where(lsm ==1)
    Simulated_Tb_low_clr[ind_land, :] = np.nan
    Simulated_Tb_high_clr[ind_land, :]= np.nan
    GMI_Tb_low[ind_land, :]          = np.nan
    GMI_Tb_high[ind_land, :]         = np.nan
    ind = np.where(GMI_Tb_low[:, 0]>1)
    lat = lat[ind]
    lon = lon[ind]
    Sim_low_clr   = np.squeeze(Simulated_Tb_low_clr[ind, :])     # (profiles, channels)
    Sim_low_clr_final = np.concatenate((Sim_low_clr_final, Sim_low_clr), 0)
    Sim_high_clr  = np.squeeze(Simulated_Tb_high_clr[ind, :])    # (profiles, channels)
    Sim_high_clr_final = np.concatenate((Sim_high_clr_final, Sim_high_clr), 0)

Sim_low_clr_final[0, :]  = np.nan
Sim_high_clr_final[0, :] = np.nan

# polarization differences at 183 Ghz

PD_sim = (Sim_low_final[:, 5, 4]- Sim_low_final[:, 6, 4])/(Sim_low_clr_final[:, 5]-Sim_low_clr_final[:, 6])
CA_sim = 1-PD_sim

PD_obs = (GMI_low_final[:, 5]- GMI_low_final[:, 6])/(Sim_low_clr_final[:, 5]-Sim_low_clr_final[:, 6])
CA_obs = 1-PD_obs
CA_avg = (CA_obs+CA_sim)/2

nbins = 22
dbin = 0.05
cmin = -0.05
cloud_amount_std_89V    = np.zeros((nbins, 2), dtype=np.float64)
cloud_amount_std_166V   = np.zeros((nbins, 2), dtype=np.float64)
cloud_amount_std_183_3V = np.zeros((nbins, 2), dtype=np.float64)
cloud_amount_std_183_7V = np.zeros((nbins, 2), dtype=np.float64)

for ibin in range (nbins):
    ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+dbin))
    FG_sub = FG_89V[ind]
    cloud_amount_std_89V[ibin, 0] = cmin
    cloud_amount_std_89V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dbin


cmin = -0.05
for ibin in range (nbins):
    ind = np.where((CA_avg>=cmin)*(CA_avg<cmin+dbin))
    FG_sub = FG_166V[ind]
    cloud_amount_std_166V[ibin, 0] = cmin
    cloud_amount_std_166V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dbin


cmin = -0.05
for ibin in range (nbins):
    ind = np.where((CA_avg>=cmin)*(CA_avg<cmin+dbin))
    FG_sub = FG_183_3V[ind]
    cloud_amount_std_183_3V[ibin, 0] = cmin
    cloud_amount_std_183_3V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dbin

cmin = -0.05
for ibin in range (nbins):
    ind = np.where((CA_avg>=cmin)*(CA_avg<cmin+dbin))
    FG_sub = FG_183_7V[ind]
    cloud_amount_std_183_7V[ibin, 0] = cmin
    cloud_amount_std_183_7V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dbin



# Normalized PDF
# for 89 V
clr = 0.0
cld_89V = 0.55
gclr_89V= 3.44
gcld_89V = 52.56

ind_89V_1 = np.where(CA_avg<clr)
FG_89V_1 = FG_89V[ind_89V_1]
FG_89V_1 = FG_89V_1[~np.isnan(FG_89V_1)]
FG_89V_norm_1 = FG_89V_1/gclr_89V

ind_89V_2 = np.where(CA_avg>cld_89V)
FG_89V_2 = FG_89V[ind_89V_2]
FG_89V_2 = FG_89V_2[~np.isnan(FG_89V_2)]
FG_89V_norm_2 = FG_89V_2/gcld_89V

ind_89V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_89V))
FG_89V_3 = FG_89V[ind_89V_3]
CA_89V_3 = CA_avg[ind_89V_3]
SD_89V = gclr_89V+ ((gcld_89V-gclr_89V)/(cld_89V-clr))*(CA_89V_3-clr)
FG_89V_norm_3 = FG_89V_3/SD_89V
FG_89V_norm_3 = FG_89V_norm_3[~np.isnan(FG_89V_norm_3)]

FG_89V_norm = np.concatenate((FG_89V_norm_1, FG_89V_norm_2, FG_89V_norm_3), 0)
mu_89_norm= np.nanmean(FG_89V_norm)
sigma_89_norm = np.nanstd(FG_89V_norm)
#------------------------------------------------------------------#
# 166V
cld_166V = 0.55
gclr_166V= 4.0
gcld_166V = 68.51

ind_166V_1 = np.where(CA_avg<clr)
FG_166V_1 = FG_166V[ind_166V_1]
FG_166V_1 = FG_166V_1[~np.isnan(FG_166V_1)]
FG_166V_norm_1 = FG_166V_1/gclr_166V

ind_166V_2 = np.where(CA_avg>cld_166V)
FG_166V_2 = FG_166V[ind_166V_2]
FG_166V_2 = FG_166V_2[~np.isnan(FG_166V_2)]
FG_166V_norm_2 = FG_166V_2/gcld_166V

ind_166V_3 = np.where((CA_avg>0.2)*(CA_avg<=cld_166V))
FG_166V_3 = FG_166V[ind_166V_3]
CA_166V_3 = CA_avg[ind_166V_3]
SD_166V = gclr_166V+((gcld_166V-gclr_166V)/(cld_166V-clr))*(CA_166V_3-clr)
FG_166V_norm_3 = FG_166V_3/SD_166V
FG_166V_norm_3 = FG_166V_norm_3[~np.isnan(FG_166V_norm_3)]

FG_166V_norm = np.concatenate((FG_166V_norm_1, FG_166V_norm_2, FG_166V_norm_3), 0)

mu_166V_norm= np.mean(FG_166V_norm)
sigma_166V_norm = np.std(FG_166V_norm)
# ---------------------------------------------------------------------------------#

# for 183+3 V
cld_183_3V = 0.55
gclr_183_3V= 3.802
gcld_183_3V = 43.477

ind_183_3V_1 = np.where(CA_avg<clr)
FG_183_3V_1 = FG_183_3V[ind_183_3V_1]
FG_183_3V_1 = FG_183_3V_1[~np.isnan(FG_183_3V_1)]
FG_183_3V_norm_1 = FG_183_3V_1/gclr_183_3V

ind_183_3V_2 = np.where(CA_avg> cld_183_3V)
FG_183_3V_2 = FG_183_3V[ind_183_3V_2]
FG_183_3V_2 = FG_183_3V_2[~np.isnan(FG_183_3V_2)]
FG_183_3V_norm_2 = FG_183_3V_2/gcld_183_3V

ind_183_3V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_183_3V))
FG_183_3V_3= FG_183_3V[ind_183_3V_3]
CA_183_3V_3 = CA_avg[ind_183_3V_3]
SD_183_3V = gclr_183_3V+ ((gcld_183_3V-gclr_183_3V)/(cld_183_3V-clr))*(CA_183_3V_3-clr)
FG_183_3V_norm_3 = FG_183_3V_3/SD_183_3V
FG_183_3V_norm_3 = FG_183_3V_norm_3[~np.isnan(FG_183_3V_norm_3)]

FG_183_3V_norm = np.concatenate((FG_183_3V_norm_1, FG_183_3V_norm_2, FG_183_3V_norm_3), 0)
mu_183_3_norm= np.mean(FG_183_3V_norm)
sigma_183_3_norm = np.std(FG_183_3V_norm)
#-----------------------------------------------------------------------------#
# for 183+7 V
cld_183_7V = 0.55
gclr_183_7V= 6.43406039
gcld_183_7V = 21.528

ind_183_7V_1 = np.where(CA_avg<clr)
FG_183_7V_1 = FG_183_7V[ind_183_7V_1]
FG_183_7V_1 = FG_183_7V_1[~np.isnan(FG_183_7V_1)]
FG_183_7V_norm_1 = FG_183_7V_1/gclr_183_7V

ind_183_7V_2 = np.where(CA_avg> cld_183_7V)
FG_183_7V_2 = FG_183_7V[ind_183_7V_2]
FG_183_7V_2 = FG_183_7V_2[~np.isnan(FG_183_7V_2)]
FG_183_7V_norm_2 = FG_183_7V_2/gcld_183_7V

ind_183_7V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_183_7V))
FG_183_7V_3= FG_183_7V[ind_183_7V_3]
CA_183_7V_3 = CA_avg[ind_183_7V_3]
SD_183_7V = gclr_183_7V+ ((gcld_183_7V-gclr_183_7V)/(cld_183_7V-clr))*(CA_183_7V_3-clr)
FG_183_7V_norm_3 = FG_183_7V_3/SD_183_7V
FG_183_7V_norm_3 = FG_183_7V_norm_3[~np.isnan(FG_183_7V_norm_3)]

FG_183_7V_norm = np.concatenate((FG_183_7V_norm_1, FG_183_7V_norm_2, FG_183_7V_norm_3), 0)
mu_183_7_norm= np.mean(FG_183_7V_norm)
sigma_183_7_norm = np.std(FG_183_7V_norm)

# Absolute FG departures

fig = plt.figure(figsize=(28, 15))
ax1= plt.subplot(2,4,1)
ax1.set_title("(a) 89V, absolute")
bins_FG = np.linspace(-100,100,1000)
ax1.set(yscale="log")
histogram, bins = np.histogram(FG_89V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax1.plot(bin_centers, histogram, 'k')
data_89V = FG_89V
data_89V= data_89V[~np.isnan(data_89V)]
mu_1= np.nanmean(data_89V)
sigma_1 = np.nanstd(data_89V)
ax1.plot(bins, 1/(sigma_1 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_1)**2 / (2 * sigma_1**2) ), 'k--')
ax1.set_xlim([-100,100])
ax1.set_ylim([1e-3, 1e0])
ax1.set_xlabel("FG Departure [K]")
ax1.set_ylabel("PDF")

ax2= plt.subplot(2,4,5)
ax2.set(yscale="log")
ax2.set_title("(e) 89V, normalized")
bins = np.linspace(-6,6,240)
histogram, bins = np.histogram(FG_89V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax2.plot(bin_centers, histogram, 'k')
ax2.plot(bins, 1/(sigma_89_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_89_norm)**2 / (2 * sigma_89_norm**2) ), 'k--')
ax2.set_xlim([-10, 10])
ax2.set_ylim([1e-3, 1e0])
ax2.set_xlabel("Normalized FG Departure [K]")
ax2.set_ylabel("PDF")

ax3= plt.subplot(2,4,2)
ax3.set(yscale="log")
ax3.set_title('(b) 166V,  absolute')
histogram, bins = np.histogram(FG_166V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax3.plot(bin_centers, histogram, 'k')
data_166V = FG_166V
data_166V= data_166V[~np.isnan(data_166V)]
mu_2= np.mean(data_166V)
sigma_2 = np.std(data_166V)
median_2 = np.median(data_166V)
ax3.plot(bins, 1/(sigma_2 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_2)**2 / (2 * sigma_2**2) ), 'k--')
ax3.set_xlim([-100,100])
ax3.set_ylim([1e-3, 1e0])
ax3.set_xlabel("FG Departure [K]")
ax3.set_ylabel("PDF")

ax4= plt.subplot(2,4,6)
ax4.set(yscale="log")
ax4.set_title('(f) 166V, normalized')
bins = np.linspace(-6,6,240)
histogram, bins = np.histogram(FG_166V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax4.plot(bin_centers, histogram, 'k')
ax4.plot(bins, 1/(sigma_166V_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_166V_norm)**2 / (2 * sigma_166V_norm**2) ), 'k--')
ax4.set_xlim([-10, 10])
ax4.set_ylim([1e-3, 1e0])
ax4.set_xlabel("Normalized FG Departure [K]")
ax4.set_ylabel("PDF")

ax5= plt.subplot(2,4,3)
ax5.set(yscale="log")
ax5.set_title('(c) $183\pm\ 3V$, absolute')
histogram, bins = np.histogram(FG_183_3V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax5.plot(bin_centers, histogram, 'k')
data_183_3V = FG_183_3V
data_183_3V= data_183_3V[~np.isnan(data_183_3V)]
mu_3= np.mean(data_183_3V)
sigma_3 = np.std(data_183_3V)
ax5.plot(bins, 1/(sigma_3 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_3)**2 / (2 * sigma_3**2) ), 'k--')
ax5.set_xlim([-100,100])
ax5.set_ylim([1e-3, 1e0])
ax5.set_xlabel("FG Departure [K]")
ax5.set_ylabel("PDF")

ax6= plt.subplot(2,4,7)
ax6.set(yscale="log")
ax6.set_title('(g) $183\pm\ 3V$, normalized')
bins = np.linspace(-6,6,240)
histogram, bins = np.histogram(FG_183_3V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax6.plot(bin_centers, histogram, 'k')
ax6.plot(bins, 1/(sigma_183_3_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_183_3_norm)**2 / (2 * sigma_183_3_norm**2) ), 'k--')
ax6.set_xlim([-10, 10])
ax6.set_ylim([1e-3, 1e0])
ax6.set_xlabel("Normalized FG Departure [K]")
ax6.set_ylabel("PDF")

ax7= plt.subplot(2,4,4)
ax7.set(yscale="log")
ax7.set_title('(d) $183\pm\ 7V$, absolute')
histogram, bins = np.histogram(FG_183_7V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax7.plot(bin_centers, histogram, 'k')
data_183_7V = FG_183_7V
data_183_7V= data_183_7V[~np.isnan(data_183_7V)]
mu_4= np.mean(data_183_7V)
sigma_4 = np.std(data_183_7V)
ax7.plot(bins, 1/(sigma_4 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_4)**2 / (2 * sigma_4**2) ), 'k--')
ax7.set_xlim([-100,100])
ax7.set_ylim([1e-3, 1e0])
ax7.set_xlabel("FG Departure [K]")
ax7.set_ylabel("PDF")

ax8= plt.subplot(2,4,8)
ax8.set(yscale="log")
ax8.set_title('(g) $183\pm\ 3V$, normalized')
bins = np.linspace(-6,6,240)
histogram, bins = np.histogram(FG_183_7V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax8.plot(bin_centers, histogram, 'k')
ax8.plot(bins, 1/(sigma_183_7_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_183_7_norm)**2 / (2 * sigma_183_7_norm**2) ), 'k--')
ax8.set_xlim([-10, 10])
ax8.set_ylim([1e-3, 1e0])
ax8.set_xlabel("Normalized FG Departure [K]")
ax8.set_ylabel("PDF")

plt.savefig('D:/Python_processing/cyclone/results/error_model/PDF_high_freq.png', dpi=300, bbox_inches='tight')
plt.show()

















