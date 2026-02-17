import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib
from scipy.stats import norm
from scipy.stats import kurtosis, skew

matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 26
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.25
matplotlib.rcParams['figure.subplot.wspace'] = 0.3
matplotlib.rcParams['lines.linewidth'] = 2.5

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

GMI_low_final[0, :]  = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :]  = np.nan
Sim_high_final[0, :, :] = np.nan
n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan

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
    ind = np.where(GMI_Tb_low[:, 0]>0)
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

cloud_amount_std = np.zeros((22, 2, 13), dtype=np.float64)

FG_low_final= GMI_low_final-Sim_low_final[:, 0:9, 4]
FG_high_final= GMI_high_final-Sim_high_final[:, 9:13, 4]

for i in range (9):
    # for low frequency
    FG = GMI_low_final[:, i]-Sim_low_final[:,i,4]
    cmin = -0.05
    for ibin in range (22):
        ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+0.05))
        FG_sub = FG[ind]
        cloud_amount_std[ibin, 0, i] = cmin+0.05
        cloud_amount_std[ibin, 1, i] = np.nanstd(FG_sub)
        cmin = cmin+0.05

for i in range (4):
    # for low frequency
    FG = GMI_high_final[:, i]-Sim_high_final[:,i+9,4]
    cmin = -0.05
    for ibin in range (22):
        ind = np.where((CA_avg>=cmin)*(CA_avg<cmin+0.05))
        FG_sub = FG[ind]
        cloud_amount_std[ibin, 0, i+9] = cmin+0.05
        cloud_amount_std[ibin, 1, i+9] = np.nanstd(FG_sub)
        cmin = cmin+0.05

# Normalized PDF
# for 19 V
clr = 0.0
cld_19V = 0.55
gclr_19V= 1.824
gcld_19V = 28.204
ind_19V_1 = np.where(CA_avg<clr)
FG_19V_1 = FG_low_final[ind_19V_1, 2]
FG_19V_1 = FG_19V_1[~np.isnan(FG_19V_1)]
FG_19V_norm_1 = FG_19V_1/gclr_19V

ind_19V_2 = np.where(CA_avg>cld_19V)
FG_19V_2 = FG_low_final[ind_19V_2, 2]
FG_19V_2 = FG_19V_2[~np.isnan(FG_19V_2)]
FG_19V_norm_2 = FG_19V_2/gcld_19V

ind_19V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_19V))
FG_19V_3 = FG_low_final[ind_19V_3, 2]
CA_19V_3 = CA_avg[ind_19V_3]
SD_19V = gclr_19V+ ((gcld_19V-gclr_19V)/(cld_19V-clr))*(CA_19V_3-clr)
FG_19V_norm_3 = FG_19V_3/SD_19V
FG_19V_norm_3 = FG_19V_norm_3[~np.isnan(FG_19V_norm_3)]

FG_19V_norm = np.concatenate((FG_19V_norm_1, FG_19V_norm_2, FG_19V_norm_3), 0)

# for 23 V
cld_23V = 0.4
gclr_23V= 2.54
gcld_23V = 9.873

ind_23V_1 = np.where(CA_avg<clr)
FG_23V_1 = FG_low_final[ind_23V_1, 4]
FG_23V_1 = FG_23V_1[~np.isnan(FG_23V_1)]
FG_23V_norm_1 = FG_23V_1/gclr_23V

ind_23V_2 = np.where(CA_avg>cld_23V)
FG_23V_2 = FG_low_final[ind_23V_2, 4]
FG_23V_2 = FG_23V_2[~np.isnan(FG_23V_2)]
FG_23V_norm_2 = FG_23V_2/gcld_23V

ind_23V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_23V))
FG_23V_3 = FG_low_final[ind_23V_3, 4]
CA_23V_3 = CA_avg[ind_23V_3]
SD_23V = gclr_23V+ ((gcld_23V-gclr_23V)/(cld_23V-clr))*(CA_23V_3-clr)
FG_23V_norm_3 = FG_23V_3/SD_23V
FG_23V_norm_3 = FG_23V_norm_3[~np.isnan(FG_23V_norm_3)]

FG_23V_norm = np.concatenate((FG_23V_norm_1, FG_23V_norm_2, FG_23V_norm_3), 0)

# for 37 V
cld_37V  = 0.40
gclr_37V = 1.979
gcld_37V = 15.98739955

ind_37V_1 = np.where(CA_avg<clr)
FG_37V_1 = FG_low_final[ind_37V_1, 5]
FG_37V_1 = FG_37V_1[~np.isnan(FG_37V_1)]
FG_37V_norm_1 = FG_37V_1/gclr_37V

ind_37V_2 = np.where(CA_avg>cld_37V)
FG_37V_2 = FG_low_final[ind_37V_2, 5]
FG_37V_2 = FG_37V_2[~np.isnan(FG_37V_2)]
FG_37V_norm_2 = FG_37V_2/gcld_37V

ind_37V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_37V))
FG_37V_3 = FG_low_final[ind_37V_3, 5]
CA_37V_3 = CA_avg[ind_37V_3]
SD_37V = gclr_37V+ ((gcld_37V-gclr_37V)/(cld_37V-clr))*(CA_37V_3-clr)
FG_37V_norm_3 = FG_37V_3/SD_37V
FG_37V_norm_3 = FG_37V_norm_3[~np.isnan(FG_37V_norm_3)]

FG_37V_norm = np.concatenate((FG_37V_norm_1, FG_37V_norm_2, FG_37V_norm_3), 0)

#-----------Absolute Values-------------------------#
FG_19V = FG_low_final[:,2]
FG_19V = FG_19V[~np.isnan(FG_19V)]
FG_23V = FG_low_final[:,4]
FG_23V = FG_23V[~np.isnan(FG_23V)]
FG_37V = FG_low_final[:,5]
FG_37V= FG_37V[~np.isnan(FG_37V)]

#----------Normalized Values-------------------------#
FG_19V_norm[FG_19V_norm>2.5]= np.nan
FG_19V_norm[FG_19V_norm<-2.5]= np.nan
FG_19V_norm = FG_19V_norm[~np.isnan(FG_19V_norm)]

FG_23V_norm[FG_23V_norm>2.5]= np.nan
FG_23V_norm[FG_23V_norm<-2.5]= np.nan
FG_23V_norm = FG_23V_norm[~np.isnan(FG_23V_norm)]

FG_37V_norm[FG_37V_norm>2.5]= np.nan
FG_37V_norm[FG_37V_norm<-2.5]= np.nan
FG_37V_norm = FG_37V_norm[~np.isnan(FG_37V_norm)]

# 19 absolute and normalized
mu_19= np.nanmean(FG_19V)
sigma_19 = np.nanstd(FG_19V)
median_19 = np.nanmedian(FG_19V)

mu_19_norm= np.nanmean(FG_19V_norm)
sigma_19_norm = np.nanstd(FG_19V_norm)
median_19_norm = np.nanmedian(FG_19V_norm)

# 23 absolute and normalized
mu_23= np.nanmean(FG_23V)
sigma_23 = np.nanstd(FG_23V)
median_23 = np.nanmedian(FG_23V)

mu_23_norm= np.nanmean(FG_23V_norm)
sigma_23_norm = np.nanstd(FG_23V_norm)
median_23_norm = np.nanmedian(FG_23V_norm)

# 37 absolute and normalized
mu_37= np.nanmean(FG_37V)
sigma_37 = np.nanstd(FG_37V)
median_37 = np.nanmedian(FG_37V)

mu_37_norm= np.nanmean(FG_37V_norm)
sigma_37_norm = np.nanstd(FG_37V_norm)
median_37_norm = np.nanmedian(FG_37V_norm)

#---------Plot---------------------#
fig = plt.figure(figsize=(21,13))
ax1= plt.subplot(2,3,1)
ax1.set(yscale="log")
ax1.set_title("(a) 19V, absolute")
bins_FG = np.linspace(-50, 50, 100)
histogram, bins = np.histogram(FG_low_final[:,2], bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax1.plot(bin_centers, histogram, 'k')
ax1.plot(bins, 1/(sigma_19 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_19)**2 / (2 * sigma_19**2) ), 'k--')
ax1.set_xlim([-50,50])
ax1.set_ylim([1e-3, 1e0])
ax1.set_xlabel("FG Departure [K]")
ax1.set_ylabel("PDF")
skew_19V = skew(FG_19V)
kurtosis_19V = kurtosis(FG_19V)

ax2= plt.subplot(2,3,4)
ax2.set(yscale="log")
ax2.set_title("(d) 19V, normalized")
bins = np.linspace(-6,6,100)
histogram, bins = np.histogram(FG_19V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax2.plot(bin_centers, histogram, 'k')
ax2.plot(bins, 1/(sigma_19_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_19_norm)**2 / (2 * sigma_19_norm**2) ), 'k--')
ax2.set_xlim([-6, 6])
ax2.set_ylim([1e-3, 1e0])
ax2.set_xlabel("Normalized FG Departure [K]")
ax2.set_ylabel("PDF")
skew_19V_norm = skew(FG_19V_norm)
kurtosis_19V_norm = kurtosis(FG_19V_norm)

ax3= plt.subplot(2,3,2)
ax3.set(yscale="log")
ax3.set_title("(b) 23V, absolute")
histogram, bins = np.histogram(FG_low_final[:,4], bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax3.plot(bin_centers, histogram, 'k')
ax3.plot(bins, 1/(sigma_23 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_23)**2 / (2 * sigma_23**2) ), 'k--')
ax3.set_xlim([-50,50])
ax3.set_ylim([1e-3, 1e0])
ax3.set_xlabel("FG Departure [K]")
ax3.set_ylabel("PDF")
skew_23V = skew(FG_23V)
kurtosis_23V = kurtosis(FG_23V)

ax4= plt.subplot(2,3,5)
ax4.set(yscale="log")
ax4.set_title("(e) 23V, normalized")
bins = np.linspace(-6,6,100)
histogram, bins = np.histogram(FG_23V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax4.plot(bin_centers, histogram, 'k')
ax4.plot(bins, 1/(sigma_23_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_23_norm)**2 / (2 * sigma_23_norm**2) ), 'k--')
ax4.set_xlim([-6, 6])
ax4.set_ylim([1e-3, 1e0])
ax4.set_xlabel("Normalized FG Departure [K]")
ax4.set_ylabel("PDF")
skew_23V_norm = skew(FG_23V_norm)
kurtosis_23V_norm = kurtosis(FG_23V_norm)

ax5= plt.subplot(2,3,3)
ax5.set(yscale="log")
ax5.set_title("(c) 37V, absolute")
histogram, bins = np.histogram(FG_low_final[:,5], bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax5.plot(bin_centers, histogram, 'k')
ax5.plot(bins, 1/(sigma_37 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_37)**2 / (2 * sigma_37**2) ), 'k--')
ax5.set_xlim([-50,50])
ax5.set_ylim([1e-3, 1e0])
ax5.set_xlabel("FG Departure [K]")
ax5.set_ylabel("PDF")
skew_37V = skew(FG_37V)
kurtosis_37V = kurtosis(FG_37V)

ax6= plt.subplot(2,3,6)
ax6.set(yscale="log")
ax6.set_title("(f) 37V, normalized")
bins = np.linspace(-6,6,100)
histogram, bins = np.histogram(FG_37V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax6.plot(bin_centers, histogram, 'k')
ax6.plot(bins, 1/(sigma_37_norm * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_37_norm)**2 / (2 * sigma_37_norm**2) ), 'k--')
ax6.set_xlim([-6, 6])
ax6.set_ylim([1e-3, 1e0])
ax6.set_xlabel("Normalized FG Departure [K]")
ax6.set_ylabel("PDF")
skew_37V_norm = skew(FG_37V_norm)
kurtosis_37V_norm = kurtosis(FG_37V_norm)

plt.savefig('D:/Python_processing/cyclone/results/error_model/PDF_low_freq.png', dpi=300, bbox_inches='tight')
plt.show()

















