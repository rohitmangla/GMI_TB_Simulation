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

cloud_amount_std_19V  = np.zeros((36,2), dtype=np.float64)
cloud_amount_std_23V  = np.zeros((36,2), dtype=np.float64)
cloud_amount_std_37V  = np.zeros((36,2), dtype=np.float64)
dkel = 2.5
# for high frequency
FG_19V = GMI_low_final[:, 2]-Sim_low_final[:,2,4]
FG_23V = GMI_low_final[:, 4]-Sim_low_final[:,4,4]
FG_37V = GMI_low_final[:, 5]-Sim_low_final[:,4,4]
cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_19V[ind]
    cloud_amount_std_19V[ibin, 0] = cmin
    cloud_amount_std_19V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel
cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_23V[ind]
    cloud_amount_std_23V[ibin, 0] = cmin
    cloud_amount_std_23V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel
cmin = -10
for ibin in range (36):
    ind = np.where((Csym>cmin)*(Csym<=cmin+dkel))
    FG_sub = FG_37V[ind]
    cloud_amount_std_37V[ibin, 0] = cmin
    cloud_amount_std_37V[ibin, 1] = np.nanstd(FG_sub)
    cmin = cmin+dkel

# Standard deviation diagram
fig = plt.figure(figsize=(8, 5))
ax1= plt.subplot(1,1,1)
ax1.plot(cloud_amount_std_19V[:,0], cloud_amount_std_19V[:,1], dashes=[3, 10, 1, 10], color= 'black', label='19V')
ax1.plot(cloud_amount_std_23V[:,0], cloud_amount_std_23V[:,1], 'k--', label= '23V')
ax1.plot(cloud_amount_std_37V[:,0], cloud_amount_std_37V[:,1], 'k-.', label= '37V')
ax1.set_xlabel('$ C_{Sym} $')
ax1.set_ylabel('Std (K)')
ax1.legend(loc= 'upper right', ncol=2)
plt.savefig('D:/Python_processing/cyclone/results/Sdv_curve_low.png', dpi=300, bbox_inches= 'tight')

#-----------------------------------#
clr =0.0
cld_19V = 20.0
gclr_19V= 3.76
gcld_19V = 19.22

cld_23V = 20.0
gclr_23V= 2.28
gcld_23V = 8.96

cld_37V = 20.0
gclr_37V= 4.97
gcld_37V = 14.61

#------------------------------------------------------#
ind_19V_1 = np.where(Csym<clr)
FG_19V_1 = FG_19V[ind_19V_1]
FG_19V_1 = FG_19V_1[~np.isnan(FG_19V_1)]
FG_19V_norm_1 = FG_19V_1/gclr_19V

ind_23V_1 = np.where(Csym<clr)
FG_23V_1 = FG_23V[ind_23V_1]
FG_23V_1 = FG_23V_1[~np.isnan(FG_23V_1)]
FG_23V_norm_1 = FG_23V_1/gclr_23V

ind_37V_1 = np.where(Csym<clr)
FG_37V_1 = FG_37V[ind_37V_1]
FG_37V_1 = FG_37V_1[~np.isnan(FG_37V_1)]
FG_37V_norm_1 = FG_37V_1/gclr_37V

#-------------------------#
ind_19V_2 = np.where(Csym> cld_19V)
FG_19V_2 = FG_19V[ind_19V_2]
FG_19V_2 = FG_19V_2[~np.isnan(FG_19V_2)]
FG_19V_norm_2 = FG_19V_2/gcld_19V

ind_23V_2 = np.where(Csym> cld_23V)
FG_23V_2 = FG_23V[ind_23V_2]
FG_23V_2 = FG_23V_2[~np.isnan(FG_23V_2)]
FG_23V_norm_2 = FG_23V_2/gcld_23V

ind_37V_2 = np.where(Csym> cld_37V)
FG_37V_2 = FG_37V[ind_37V_2]
FG_37V_2 = FG_37V_2[~np.isnan(FG_37V_2)]
FG_37V_norm_2 = FG_37V_2/gcld_37V

#-------------------------#
ind_19V_3 = np.where((Csym>=clr)*(Csym<=cld_19V))
FG_19V_3= FG_19V[ind_19V_3]
CA_19V_3 = Csym[ind_19V_3]
SD_19V = gclr_19V+ (gcld_19V-gclr_19V)*(((CA_19V_3-clr)/(cld_19V-clr))**2)
FG_19V_norm_3 = FG_19V_3/SD_19V
FG_19V_norm_3 = FG_19V_norm_3[~np.isnan(FG_19V_norm_3)]

ind_23V_3 = np.where((Csym>=clr)*(Csym<=cld_23V))
FG_23V_3= FG_23V[ind_23V_3]
CA_23V_3 = Csym[ind_23V_3]
SD_23V = gclr_23V+ (gcld_23V-gclr_23V)*(((CA_23V_3-clr)/(cld_23V-clr))**2)
FG_23V_norm_3 = FG_23V_3/SD_23V
FG_23V_norm_3 = FG_23V_norm_3[~np.isnan(FG_23V_norm_3)]

ind_37V_3 = np.where((Csym>=clr)*(Csym<=cld_37V))
FG_37V_3= FG_37V[ind_37V_3]
CA_37V_3 = Csym[ind_37V_3]
SD_37V = gclr_37V+ (gcld_37V-gclr_37V)*(((CA_37V_3-clr)/(cld_37V-clr))**2)
FG_37V_norm_3 = FG_37V_3/SD_37V
FG_37V_norm_3 = FG_37V_norm_3[~np.isnan(FG_37V_norm_3)]

#---------------------#
FG_19V_norm = np.concatenate((FG_19V_norm_1, FG_19V_norm_2, FG_19V_norm_3), 0)
FG_23V_norm = np.concatenate((FG_23V_norm_1, FG_23V_norm_2, FG_23V_norm_3), 0)
FG_37V_norm = np.concatenate((FG_37V_norm_1, FG_37V_norm_2, FG_37V_norm_3), 0)

mu_1= np.mean(FG_19V_norm)
sigma_1 = np.std(FG_19V_norm)
median_1 = np.median(FG_19V_norm)

mu_2= np.mean(FG_23V_norm)
sigma_2 = np.std(FG_23V_norm)
median_2 = np.median(FG_23V_norm)

mu_3= np.mean(FG_37V_norm)
sigma_3 = np.std(FG_37V_norm)
median_3 = np.median(FG_37V_norm)

#----------------------------------#

fig= plt.figure(figsize=(20,4.5))
ax0= plt.subplot(1,3,1)
ax0.set(yscale="log")
ax0.set_title('(a) 19V, normalized')
bins = np.linspace(-6,6,30)
histogram, bins = np.histogram(FG_19V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax0.plot(bin_centers, histogram, 'k')
ax0.plot(bins, 1/(sigma_1 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_1)**2 / (2 * sigma_1**2) ), 'k--')
ax0.set_xlim([-8, 8])
ax0.set_ylim([1e-3, 1])
ax0.set_xlabel("Normalized FG Departure [K]")
ax0.set_ylabel("PDF")

ax1= plt.subplot(1,3,2)
ax1.set(yscale="log")
ax1.set_title('(b) 23V, normalized')
bins = np.linspace(-6,6,30)
histogram, bins = np.histogram(FG_23V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax1.plot(bin_centers, histogram, 'k')
ax1.plot(bins, 1/(sigma_2 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_2)**2 / (2 * sigma_2**2) ), 'k--')
ax1.set_xlim([-8, 8])
ax1.set_ylim([1e-3, 1])
ax1.set_xlabel("Normalized FG Departure [K]")
ax1.set_ylabel("PDF")

ax2= plt.subplot(1,3,3)
ax2.set(yscale="log")
ax2.set_title('(c) 37V, normalized')
bins = np.linspace(-6,6,30)
histogram, bins = np.histogram(FG_37V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers)
ax2.plot(bin_centers, histogram, 'k')
ax2.plot(bins, 1/(sigma_3 * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu_3)**2 / (2 * sigma_3**2) ), 'k--')
ax2.set_xlim([-8, 8])
ax2.set_ylim([1e-3, 1])
ax2.set_xlabel("Normalized FG Departure [K]")
ax2.set_ylabel("PDF")

plt.savefig('D:/Python_processing/cyclone/results/error_model/PDF_high_SI_low.png', dpi=300, bbox_inches='tight')
plt.show()

















