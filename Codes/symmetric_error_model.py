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
matplotlib.rcParams['figure.subplot.hspace'] = 0.25
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
gcld_23V = 9.87

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
cld_37V  = 0.45
gclr_37V = 1.979
gcld_37V = 15.03

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

# for 89 V
cld_89V = 0.5
gclr_89V = 2.37124962
gcld_89V = 35.20003476

ind_89V_1 = np.where(CA_avg<clr)
FG_89V_1 = FG_low_final[ind_89V_1, 7]
FG_89V_1 = FG_89V_1[~np.isnan(FG_89V_1)]
FG_89V_norm_1 = FG_89V_1/gclr_89V

ind_89V_2 = np.where(CA_avg>cld_89V)
FG_89V_2 = FG_low_final[ind_89V_2, 7]
FG_89V_2 = FG_89V_2[~np.isnan(FG_89V_2)]
FG_89V_norm_2 = FG_89V_2/gcld_89V

ind_89V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_89V))
FG_89V_3 = FG_low_final[ind_89V_3, 7]
CA_89V_3 = CA_avg[ind_89V_3]
SD_89V = gclr_89V+ ((gcld_89V-gclr_89V)/(cld_89V-clr))*(CA_89V_3-clr)
FG_89V_norm_3 = FG_89V_3/SD_89V
FG_89V_norm_3 = FG_89V_norm_3[~np.isnan(FG_89V_norm_3)]

FG_89V_norm = np.concatenate((FG_89V_norm_1, FG_89V_norm_2, FG_89V_norm_3), 0)

# for 166 V
cld_166V = 0.55
gclr_166V= 5.76937365
gcld_166V = 55.51169061

ind_166V_1 = np.where(CA_avg<clr)
FG_166V_1 = FG_high_final[ind_166V_1, 0]
FG_166V_1 = FG_166V_1[~np.isnan(FG_166V_1)]
FG_166V_norm_1 = FG_166V_1/gclr_166V

ind_166V_2 = np.where(CA_avg>cld_166V)
FG_166V_2 = FG_high_final[ind_166V_2, 0]
FG_166V_2 = FG_166V_2[~np.isnan(FG_166V_2)]
FG_166V_norm_2 = FG_166V_2/gcld_166V

ind_166V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_166V))
FG_166V_3 = FG_high_final[ind_166V_3, 0]
CA_166V_3 = CA_avg[ind_166V_3]
SD_166V = gclr_166V+ ((gcld_166V-gclr_166V)/(cld_166V-clr))*(CA_166V_3-clr)
FG_166V_norm_3 = FG_166V_3/SD_166V
FG_166V_norm_3 = FG_166V_norm_3[~np.isnan(FG_166V_norm_3)]

FG_166V_norm = np.concatenate((FG_166V_norm_1, FG_166V_norm_2, FG_166V_norm_3), 0)


# for 183+3 V
cld_183_3V = 0.5
gclr_183_3V= 3.47142742
gcld_183_3V = 34.02921401

ind_183_3V_1 = np.where(CA_avg<clr)
FG_183_3V_1 = FG_high_final[ind_183_3V_1, 2]
FG_183_3V_1 = FG_183_3V_1[~np.isnan(FG_183_3V_1)]
FG_183_3V_norm_1 = FG_183_3V_1/gclr_183_3V

ind_183_3V_2 = np.where(CA_avg> cld_183_3V)
FG_183_3V_2 = FG_high_final[ind_183_3V_2, 2]
FG_183_3V_2 = FG_183_3V_2[~np.isnan(FG_183_3V_2)]
FG_183_3V_norm_2 = FG_183_3V_2/gcld_183_3V

ind_183_3V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_183_3V))
FG_183_3V_3= FG_high_final[ind_183_3V_3, 2]
CA_183_3V_3 = CA_avg[ind_183_3V_3]
SD_183_3V = gclr_183_3V+ ((gcld_183_3V-gclr_183_3V)/(cld_183_3V-clr))*(CA_183_3V_3-clr)
FG_183_3V_norm_3 = FG_183_3V_3/SD_183_3V
FG_183_3V_norm_3 = FG_183_3V_norm_3[~np.isnan(FG_183_3V_norm_3)]

FG_183_3V_norm = np.concatenate((FG_183_3V_norm_1, FG_183_3V_norm_2, FG_183_3V_norm_3), 0)

# for 183+7 V
cld_183_7V = 0.55
gclr_183_7V= 5.06431253
gcld_183_7V = 48.22394227

ind_183_7V_1 = np.where(CA_avg<clr)
FG_183_7V_1 = FG_high_final[ind_183_7V_1, 3]
FG_183_7V_1 = FG_183_7V_1[~np.isnan(FG_183_7V_1)]
FG_183_7V_norm_1 = FG_183_7V_1/gclr_183_7V

ind_183_7V_2 = np.where(CA_avg> cld_183_7V)
FG_183_7V_2 = FG_high_final[ind_183_7V_2, 3]
FG_183_7V_2 = FG_183_7V_2[~np.isnan(FG_183_7V_2)]
FG_183_7V_norm_2 = FG_183_7V_2/gcld_183_7V

ind_183_7V_3 = np.where((CA_avg>=clr)*(CA_avg<=cld_183_7V))
FG_183_7V_3= FG_high_final[ind_183_7V_3, 3]
CA_183_7V_3 = CA_avg[ind_183_7V_3]
SD_183_7V = gclr_183_7V+ ((gcld_183_7V-gclr_183_7V)/(cld_183_7V-clr))*(CA_183_7V_3-clr)
FG_183_7V_norm_3 = FG_183_7V_3/SD_183_7V
FG_183_7V_norm_3 = FG_183_7V_norm_3[~np.isnan(FG_183_7V_norm_3)]

FG_183_7V_norm = np.concatenate((FG_183_7V_norm_1, FG_183_7V_norm_2, FG_183_7V_norm_3), 0)

fig = plt.figure(figsize=(27, 14))
ax1= plt.subplot(2,4,1)
ax1.set_title("19 V")
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax1.set(yscale="log")
histogram, bins = np.histogram(FG_19V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax1.plot(bin_centers, histogram, 'k')
ax1.plot(bin_centers, pdf, 'k--')
ax1.set_xlim([-6, 6])
ax1.set_ylim([1e-3, 1e0])
ax1.set_xlabel("Normalized FG Departures")
ax1.set_ylabel("PDF")

ax2= plt.subplot(2,4,2)
ax2.set_title("23 V")
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax2.set(yscale="log")
histogram, bins = np.histogram(FG_23V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax2.plot(bin_centers, histogram, 'k')
ax2.plot(bin_centers, pdf, 'k--')
ax2.set_xlim([-8, 8])
ax2.set_ylim([1e-3, 1e0])
ax2.set_xlabel("Normalized FG Departures")
ax2.set_ylabel("PDF")

ax3= plt.subplot(2,4,3)
ax3.set_title("37 V")
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax3.set(yscale="log")
histogram, bins = np.histogram(FG_37V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax3.plot(bin_centers, histogram, 'k')
ax3.plot(bin_centers, pdf, 'k--')
ax3.set_xlim([-8, 8])
ax3.set_ylim([1e-3, 1e0])
ax3.set_xlabel("Normalized FG Departures")
ax3.set_ylabel("PDF")

ax4= plt.subplot(2,4,4)
ax4.set_title("89 V")
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax4.set(yscale="log")
histogram, bins = np.histogram(FG_89V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax4.plot(bin_centers, histogram, 'k')
ax4.plot(bin_centers, pdf, 'k--')
ax4.set_xlim([-8, 8])
ax4.set_ylim([1e-3, 1e0])
ax4.set_xlabel("Normalized FG Departures")
ax4.set_ylabel("PDF")

ax5= plt.subplot(2,4,5)
ax5.set_title("166 V")
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax5.set(yscale="log")
histogram, bins = np.histogram(FG_166V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax5.plot(bin_centers, histogram, 'k')
ax5.plot(bin_centers, pdf, 'k--')
ax5.set_xlim([-8, 8])
ax5.set_ylim([1e-3, 1.5])
ax5.set_xlabel("Normalized FG Departures")
ax5.set_ylabel("PDF")

ax6= plt.subplot(2,4,6)
ax6.set_title('$183\pm\ 3V$')
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax6.set(yscale="log")
histogram, bins = np.histogram(FG_183_3V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax6.plot(bin_centers, histogram, 'k')
ax6.plot(bin_centers, pdf, 'k--')
ax6.set_xlim([-8, 8])
ax6.set_ylim([1e-3, 1.5])
ax6.set_xlabel("Normalized FG Departures")
ax6.set_ylabel("PDF")

ax7= plt.subplot(2,4,7)
ax7.set_title('$183\pm\ 7V$')
# Compute a histogram of the sample
bins = np.linspace(-6, 6, 240 )
ax7.set(yscale="log")
histogram, bins = np.histogram(FG_183_7V_norm, bins=bins, normed=True)
bin_centers = 0.5*(bins[1:] + bins[:-1])
# Compute the PDF on the bin centers from scipy distribution object
pdf = norm.pdf(bin_centers)
ax7.plot(bin_centers, histogram, 'k')
ax7.plot(bin_centers, pdf, 'k--')
ax7.set_xlim([-8, 8])
ax7.set_ylim([1e-3, 1.5])
ax7.set_xlabel("Normalized FG Departures")
ax7.set_ylabel("PDF")
plt.savefig('D:/Python_processing/cyclone/results/error_model/norm_pdf.png', dpi=300, bbox_inches='tight')

