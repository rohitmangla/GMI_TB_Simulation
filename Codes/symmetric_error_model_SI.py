import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib
from scipy.stats import norm, kurtosis, skew

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
    ind = np.where(GMI_Tb_high[:, 0]>1)
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
FG_low_final= GMI_low_final-Sim_low_final[:, 0:9, 4]
FG_high_final= GMI_high_final-Sim_high_final[:, 9:13, 4]

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
    ind = np.where(GMI_Tb_high[:, 0]>1)
    lat = lat[ind]
    lon = lon[ind]
    Sim_low_clr   = np.squeeze(Simulated_Tb_low_clr[ind, :])     # (profiles, channels)
    Sim_low_clr_final = np.concatenate((Sim_low_clr_final, Sim_low_clr), 0)
    Sim_high_clr  = np.squeeze(Simulated_Tb_high_clr[ind, :])    # (profiles, channels)
    Sim_high_clr_final = np.concatenate((Sim_high_clr_final, Sim_high_clr), 0)

Sim_low_clr_final[0, :]  = np.nan
Sim_high_clr_final[0, :] = np.nan

# Symmetric Index
SIOB = (GMI_low_final[:, 8]-GMI_high_final[:, 0])-(Sim_low_clr_final[:, 8]-Sim_high_clr_final[:, 9])
SIFG = (Sim_low_final[:, 8, 4]-Sim_high_final[:, 9, 4])-(Sim_low_clr_final[:, 8]-Sim_high_clr_final[:, 9])
Csym = (SIOB+SIFG)/2

cloud_amount_std_89V    = np.zeros((36,2), dtype=np.float64)
cloud_amount_std_166V   = np.zeros((36,2), dtype=np.float64)
cloud_amount_std_183_3V = np.zeros((36,2), dtype=np.float64)
cloud_amount_std_183_7V = np.zeros((36,2), dtype=np.float64)
dkel = 2.5
# for high frequency
FG_89V = GMI_low_final[:, 7]-Sim_low_final[:,7,4]
FG_166V = GMI_high_final[:, 0]-Sim_high_final[:,9,4]
FG_183_3V = GMI_high_final[:, 2]-Sim_high_final[:,11,4]
FG_183_7V = GMI_high_final[:, 3]-Sim_high_final[:,12,4]


cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_89V[ind]
    cloud_amount_std_89V[ibin, 0] = cmin
    cloud_amount_std_89V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel
cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_166V[ind]
    cloud_amount_std_166V[ibin, 0] = cmin
    cloud_amount_std_166V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel
cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_183_3V[ind]
    cloud_amount_std_183_3V[ibin, 0] = cmin
    cloud_amount_std_183_3V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel
cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_183_7V[ind]
    cloud_amount_std_183_7V[ibin, 0] = cmin
    cloud_amount_std_183_7V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel

# Standard deviation diagram
fig = plt.figure(figsize=(8, 5))
ax1= plt.subplot(1,1,1)
ax1.plot(cloud_amount_std_89V[0:32,0], cloud_amount_std_89V[0:32,1], dashes=[3, 10, 1, 10], color= 'black', label='89V')
ax1.plot(cloud_amount_std_166V[0:32,0], cloud_amount_std_166V[0:32,1], 'k--', label= '166V')
ax1.plot(cloud_amount_std_183_3V[0:32,0], cloud_amount_std_183_3V[0:32,1], 'k-.', label= '$183\pm\ 3V$')
ax1.plot(cloud_amount_std_183_7V[0:32,0], cloud_amount_std_183_7V[0:32,1], 'k:',label= '$183\pm\ 7V$')
ax1.set_xlabel('$ SI_{avg} $')
ax1.set_ylabel('Std (K)')
ax1.legend(loc= 'upper right', ncol=2)
ax1.set_ylim([0, 80])
plt.savefig('D:/Python_processing/cyclone/results/error_model/Sdv_curve.png', dpi=300, bbox_inches= 'tight')

#-----------------------------------#
clr =  3
cld_89V =  25
#gclr_89V = 6.27
gclr_89V= 3.35
gcld_89V = 45.15

cld_166V = 25.0
gclr_166V= 3.51
gcld_166V = 70.77

cld_183_3V = 30.0
gclr_183_3V= 2.81
gcld_183_3V = 42.43

cld_183_7V = 27.5
gclr_183_7V= 3.14
gcld_183_7V = 61.0
#------------------------------------------------------#
ind_89V_1 = np.where(Csym<=clr)
FG_89V_1 = FG_89V[ind_89V_1]
FG_89V_1 = FG_89V_1[~np.isnan(FG_89V_1)]
FG_89V_norm_1 = FG_89V_1/gclr_89V

ind_166V_1 = np.where(Csym<=clr)
FG_166V_1 = FG_166V[ind_166V_1]
FG_166V_1 = FG_166V_1[~np.isnan(FG_166V_1)]
FG_166V_norm_1 = FG_166V_1/gclr_166V

ind_183_3V_1 = np.where(Csym<=clr)
FG_183_3V_1 = FG_183_3V[ind_183_3V_1]
FG_183_3V_1 = FG_183_3V_1[~np.isnan(FG_183_3V_1)]
FG_183_3V_norm_1 = FG_183_3V_1/gclr_183_3V

ind_183_7V_1 = np.where(Csym<=clr)
FG_183_7V_1 = FG_183_7V[ind_183_7V_1]
FG_183_7V_1 = FG_183_7V_1[~np.isnan(FG_183_7V_1)]
FG_183_7V_norm_1 = FG_183_7V_1/gclr_183_7V
#-------------------------#
ind_89V_2 = np.where(Csym>=cld_166V)
FG_89V_2 = FG_89V[ind_89V_2]
FG_89V_2 = FG_89V_2[~np.isnan(FG_89V_2)]
FG_89V_norm_2 = FG_89V_2/gcld_89V

ind_166V_2 = np.where(Csym>=cld_166V)
FG_166V_2 = FG_166V[ind_166V_2]
FG_166V_2 = FG_166V_2[~np.isnan(FG_166V_2)]
FG_166V_norm_2 = FG_166V_2/gcld_166V

ind_183_3V_2 = np.where(Csym>=cld_183_3V)
FG_183_3V_2 = FG_183_3V[ind_183_3V_2]
FG_183_3V_2 = FG_183_3V_2[~np.isnan(FG_183_3V_2)]
FG_183_3V_norm_2 = FG_183_3V_2/gcld_183_3V

ind_183_7V_2 = np.where(Csym>=cld_183_7V)
FG_183_7V_2 = FG_183_7V[ind_183_7V_2]
FG_183_7V_2 = FG_183_7V_2[~np.isnan(FG_183_7V_2)]
FG_183_7V_norm_2 = FG_183_7V_2/gcld_183_7V
#-------------------------#

ind_89V_3 = np.where((Csym>clr)*(Csym<cld_89V))
FG_89V_3= FG_89V[ind_89V_3]
CA_89V_3 = Csym[ind_89V_3]
SD_89V = gclr_89V+ (gcld_89V-gclr_89V)*(((CA_89V_3-clr)/(cld_89V-clr))**2)
FG_89V_norm_3 = FG_89V_3/SD_89V
FG_89V_norm_3 = FG_89V_norm_3[~np.isnan(FG_89V_norm_3)]

ind_166V_3 = np.where((Csym>clr)*(Csym<cld_166V))
FG_166V_3= FG_166V[ind_166V_3]
CA_166V_3 = Csym[ind_166V_3]
SD_166V = gclr_166V+ (gcld_166V-gclr_166V)*(((CA_166V_3-clr)/(cld_166V-clr))**2)
FG_166V_norm_3 = FG_166V_3/SD_166V
FG_166V_norm_3 = FG_166V_norm_3[~np.isnan(FG_166V_norm_3)]

ind_183_3V_3 = np.where((Csym>clr)*(Csym<cld_183_3V))
FG_183_3V_3= FG_183_3V[ind_183_3V_3]
CA_183_3V_3 = Csym[ind_183_3V_3]
SD_183_3V = gclr_183_3V+ (gcld_183_3V-gclr_183_3V)*(((CA_183_3V_3-clr)/(cld_183_3V-clr))**2)
FG_183_3V_norm_3 = FG_183_3V_3/SD_183_3V
FG_183_3V_norm_3 = FG_183_3V_norm_3[~np.isnan(FG_183_3V_norm_3)]

ind_183_7V_3 = np.where((Csym>clr)*(Csym<cld_183_7V))
FG_183_7V_3= FG_183_7V[ind_183_7V_3]
CA_183_7V_3 = Csym[ind_183_7V_3]
SD_183_7V = gclr_183_7V+ (gcld_183_7V-gclr_183_7V)*(((CA_183_7V_3-clr)/(cld_183_7V-clr))**2)
FG_183_7V_norm_3 = FG_183_7V_3/SD_183_7V
FG_183_7V_norm_3 = FG_183_7V_norm_3[~np.isnan(FG_183_7V_norm_3)]
#---------------------#

FG_89V_norm = np.concatenate((FG_89V_norm_1, FG_89V_norm_2, FG_89V_norm_3), 0)
FG_166V_norm = np.concatenate((FG_166V_norm_1, FG_166V_norm_2, FG_166V_norm_3), 0)
FG_183_3V_norm = np.concatenate((FG_183_3V_norm_1, FG_183_3V_norm_2, FG_183_3V_norm_3), 0)
FG_183_7V_norm = np.concatenate((FG_183_7V_norm_1, FG_183_7V_norm_2, FG_183_7V_norm_3), 0)

error

#-----------Absolute Values-------------------------#
FG_89V = FG_89V[~np.isnan(FG_89V)]
FG_166V = FG_166V[~np.isnan(FG_166V)]
FG_183_3V = FG_183_3V[~np.isnan(FG_183_3V)]
FG_183_7V = FG_183_7V[~np.isnan(FG_183_7V)]

#-----------Normalized values----------------------#
FG_89V_norm[FG_89V_norm>2.5]= np.nan
FG_89V_norm[FG_89V_norm<-2.5]= np.nan
FG_89V_norm = FG_89V_norm[~np.isnan(FG_89V_norm)]

FG_166V_norm[FG_166V_norm>2.5]= np.nan
FG_166V_norm[FG_166V_norm<-2.5]= np.nan
FG_166V_norm = FG_166V_norm[~np.isnan(FG_166V_norm)]

FG_183_3V_norm[FG_183_3V_norm>2.5]= np.nan
FG_183_3V_norm[FG_183_3V_norm<-2.5]= np.nan
FG_183_3V_norm = FG_183_3V_norm[~np.isnan(FG_183_3V_norm)]

FG_183_7V_norm[FG_183_7V_norm>2.5]= np.nan
FG_183_7V_norm[FG_183_7V_norm<-2.5]= np.nan
FG_183_7V_norm = FG_183_7V_norm[~np.isnan(FG_183_7V_norm)]

error
# 89 absolute and normalized
mu_0= np.nanmean(FG_89V)
sigma_0 = np.nanstd(FG_89V)
median_0 = np.nanmedian(FG_89V)

mu_1= np.nanmean(FG_89V_norm)
sigma_1 = np.nanstd(FG_89V_norm)
median_1 = np.nanmedian(FG_89V_norm)

# 166V absolute and normalized
mu_2= np.nanmean(FG_166V)
sigma_2 = np.nanstd(FG_166V)
median_2 = np.nanmedian(FG_166V)

mu_3= np.nanmean(FG_166V_norm)
sigma_3 = np.nanstd(FG_166V_norm)
median_3 = np.nanmedian(FG_166V_norm)

# 183_3 V absolute and normalized
mu_4     = np.nanmean(FG_183_3V)
sigma_4  = np.nanstd(FG_183_3V)
median_4 = np.nanmedian(FG_183_3V)

mu_5= np.nanmean(FG_183_3V_norm)
sigma_5 = np.nanstd(FG_183_3V_norm)
median_5 = np.nanmedian(FG_183_3V_norm)

# 183_7 V absolute and normalized
mu_6     = np.nanmean(FG_183_7V)
sigma_6  = np.nanstd(FG_183_7V)
median_6 = np.nanmedian(FG_183_7V)

mu_7= np.nanmean(FG_183_7V_norm)
sigma_7 = np.nanstd(FG_183_7V_norm)
median_7= np.nanmedian(FG_183_7V_norm)

#----------------------------------#
fig= plt.figure(figsize=(27,12))
ax0= plt.subplot(2,4,1)
ax0.set(yscale="log")
ax0.set_title("(a) 89V, absolute")
bins_FG = np.linspace(-100,100,1000)
histogram, bins = np.histogram(FG_89V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax0.plot(bin_centers, histogram, 'k')
ax0.plot(bins, 1/(sigma_0 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_0)**2 / (2 * sigma_0**2) ), 'k--')
ax0.set_xlim([-100,100])
ax0.set_ylim([1e-3, 1])
ax0.set_xlabel("FG Departure [K]")
ax0.set_ylabel("PDF")
skew_89V = skew(FG_89V)
kurtosis_89V = kurtosis(FG_89V)

ax1= plt.subplot(2,4,5)
ax1.set(yscale="log")
ax1.set_title('(e) 89V, normalized')
bins = np.linspace(-3, 3,300)
histogram, bins = np.histogram(FG_89V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax1.plot(bin_centers, histogram, 'k')
ax1.plot(bins, 1/(sigma_1 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_1)**2 / (2 * sigma_1**2) ), 'k--')
ax1.set_xlim([-3, 3])
ax1.set_ylim([1e-3, 1])
ax1.set_xlabel("Normalized FG Departure [K]")
ax1.set_ylabel("PDF")
skew_89V_norm = skew(FG_89V_norm)
kurtosis_89V_norm = kurtosis(FG_89V_norm)

ax2= plt.subplot(2,4,2)
ax2.set(yscale="log")
ax2.set_title("(b) 166V, absolute")
bins_FG = np.linspace(-100,100,1000)
histogram, bins = np.histogram(FG_166V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax2.plot(bin_centers, histogram, 'k')
ax2.plot(bins, 1/(sigma_2 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_2)**2 / (2 * sigma_2**2) ), 'k--')
ax2.set_xlim([-100,100])
ax2.set_ylim([1e-3, 1])
ax2.set_xlabel("FG Departure [K]")
ax2.set_ylabel("PDF")
FG_166V = FG_166V[~np.isnan(FG_166V)]
skew_166V = skew(FG_166V)
kurtosis_166V = kurtosis(FG_166V)

ax3= plt.subplot(2,4,6)
ax3.set(yscale="log")
ax3.set_title('(f) 166V, normalized')
bins = np.linspace(-3,3,300)
histogram, bins = np.histogram(FG_166V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax3.plot(bin_centers, histogram, 'k')
ax3.plot(bins, 1/(sigma_3 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_3)**2 / (2 * sigma_3**2) ), 'k--')
ax3.set_xlim([-3, 3])
ax3.set_ylim([1e-3, 1])
ax3.set_xlabel("Normalized FG Departure [K]")
ax3.set_ylabel("PDF")
skew_166V_norm = skew(FG_166V_norm)
kurtosis_166V_norm = kurtosis(FG_166V_norm)

ax4= plt.subplot(2,4,3)
ax4.set(yscale="log")
ax4.set_title("(c) $183\pm\ 3V$, absolute")
bins_FG = np.linspace(-100,100,1000)
histogram, bins = np.histogram(FG_183_3V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax4.plot(bin_centers, histogram, 'k')
ax4.plot(bins, 1/(sigma_4 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_4)**2 / (2 * sigma_4**2) ), 'k--')
ax4.set_xlim([-100,100])
ax4.set_ylim([1e-3, 1])
ax4.set_xlabel("FG Departure [K]")
ax4.set_ylabel("PDF")
FG_183_3V = FG_183_3V[~np.isnan(FG_183_3V)]
skew_183_3V = skew(FG_183_3V)
kurtosis_183_3V = kurtosis(FG_183_3V)

ax5= plt.subplot(2,4,7)
ax5.set(yscale="log")
ax5.set_title('(g) $183\pm\ 3V$, normalized')
bins = np.linspace(-3,3,300)
histogram, bins = np.histogram(FG_183_3V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax5.plot(bin_centers, histogram, 'k')
ax5.plot(bins, 1/(sigma_5 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_5)**2 / (2 * sigma_5**2) ), 'k--')
ax5.set_xlim([-3,3])
ax5.set_ylim([1e-3, 1])
ax5.set_xlabel("Normalized FG Departure [K]")
ax5.set_ylabel("PDF")
skew_183_3V_norm = skew(FG_183_3V_norm)
kurtosis_183_3V_norm = kurtosis(FG_183_3V_norm)

ax6= plt.subplot(2,4,4)
ax6.set(yscale="log")
ax6.set_title("(d) $183\pm\ 7V$, absolute")
bins_FG = np.linspace(-100,100,1000)
histogram, bins = np.histogram(FG_183_7V, bins=bins_FG, normed=True)
bin_centers = 0.5*(bins_FG[1:] + bins_FG[:-1])
pdf = norm.pdf(bin_centers)
ax6.plot(bin_centers, histogram, 'k')
ax6.plot(bins, 1/(sigma_6 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_6)**2 / (2 * sigma_6**2) ), 'k--')
ax6.set_xlim([-100,100])
ax6.set_ylim([1e-3, 1])
ax6.set_xlabel("FG Departure [K]")
ax6.set_ylabel("PDF")
FG_183_7V = FG_183_7V[~np.isnan(FG_183_7V)]
skew_183_7V = skew(FG_183_7V)
kurtosis_183_7V = kurtosis(FG_183_7V)

ax7= plt.subplot(2,4,8)
ax7.set(yscale="log")
ax7.set_title('(h) $183\pm\ 7V$, normalized')
bins = np.linspace(-3,3,300)
histogram, bins = np.histogram(FG_183_7V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax7.plot(bin_centers, histogram, 'k')
ax7.plot(bins, 1/(sigma_7 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_7)**2 / (2 * sigma_7**2) ), 'k--')
ax7.set_xlim([-3,3])
ax7.set_ylim([1e-3, 1])
ax7.set_xlabel("Normalized FG Departure [K]")
ax7.set_ylabel("PDF")
skew_183_7V_norm = skew(FG_183_7V_norm)
kurtosis_183_7V_norm = kurtosis(FG_183_7V_norm)

#plt.savefig('D:/Python_processing/cyclone/results/error_model/PDF_high_SI.png', dpi=300, bbox_inches='tight')
plt.show()

















