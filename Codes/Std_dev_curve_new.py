import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 24
matplotlib.rcParams['legend.fontsize'] =15
matplotlib.rcParams['font.size']= 24
matplotlib.rcParams['axes.labelsize']= 20
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.35
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
    Simulated_Tb_low[ind_land, :, :]= np.nan
    Simulated_Tb_high[ind_land, :, :]= np.nan
    GMI_Tb_low[ind_land, :]      = np.nan
    GMI_Tb_high[ind_land, :]     = np.nan

    ind_high = np.where(GMI_Tb_low[:, 0] > 1)
    lat_high = lat[ind_high]
    lon_high = lon[ind_high]

    GMI_low   = np.squeeze(GMI_Tb_low[ind_high, :])              # (profiles, channels_low)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    GMI_high  = np.squeeze(GMI_Tb_high[ind_high, :])             # (profiles, channels_high)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    Sim_low   = np.squeeze(Simulated_Tb_low[ind_high, :, :])     # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)
    Sim_high  = np.squeeze(Simulated_Tb_high[ind_high, :, :])    # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)

GMI_low_final[0, :] = np.nan
GMI_high_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
Sim_high_final[0, :, :] = np.nan
n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan
FG_19V =  GMI_low_final[:, 2]-Sim_low_final[:,2,4]
FG_19V_qc =  GMI_low_final[:, 3]-Sim_low_final[:,2,4]
std_19V = np.nanstd(FG_19V)
FG_19V_qc[FG_19V_qc>2*std_19V]= np.nan
FG_19V_qc[FG_19V_qc<-2*std_19V]= np.nan


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
    ind_high = np.where(GMI_Tb_low[:, 0]>1)
    Sim_low_clr   = np.squeeze(Simulated_Tb_low_clr[ind_high, :])     # (profiles, channels)
    Sim_low_clr_final = np.concatenate((Sim_low_clr_final, Sim_low_clr), 0)
    Sim_high_clr  = np.squeeze(Simulated_Tb_high_clr[ind_high, :])    # (profiles, channels)
    Sim_high_clr_final = np.concatenate((Sim_high_clr_final, Sim_high_clr), 0)

Sim_low_clr_final[0, :]  = np.nan
Sim_high_clr_final[0, :] = np.nan

# polarization differences at 183 Ghz
PD_sim = (Sim_low_final[:, 5, 4]- Sim_low_final[:, 6, 4])/(Sim_low_clr_final[:, 5]-Sim_low_clr_final[:, 6])
CA_sim = 1-PD_sim
PD_obs = (GMI_low_final[:, 5]- GMI_low_final[:, 6])/(Sim_low_clr_final[:, 5]-Sim_low_clr_final[:, 6])
CA_obs = 1-PD_obs
CA_avg = (CA_obs+CA_sim)/2

#cloud_amount_std = np.zeros((22, 4, 13), dtype=np.float64)
cloud_amount_std_19V = np.zeros((22, 4), dtype=np.float64)
cmin = -0.05
nbins= 22
dbin = 0.05
for ibin in range (nbins):
    ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+dbin))
    FG_sub = FG_19V[ind]
    FG_sub_qc = FG_19V_qc[ind]
    cloud_amount_std_19V[ibin, 0] = cmin
    cloud_amount_std_19V[ibin, 1] = np.nanstd(FG_sub)
    cloud_amount_std_19V[ibin, 3] = np.nanstd(FG_sub_qc)
    cmin = cmin+dbin
#for i in range (9):
    # for low frequency
#    FG = GMI_low_final[:, i]-Sim_low_final[:,i,4]
#    stdk = np.nanstd(FG)
#    FG_qc = FG
#    FG_qc[FG_qc>2*stdk] = np.nan
#    FG_qc[FG_qc<-2*stdk] = np.nan
#    cmin = -0.075
#    for ibin in range (22):
#        ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+0.05))
#        FG_sub = FG[ind]
#        FG_sub_qc = FG_qc[ind]
#        count = np.count_nonzero(FG_sub)
#        cloud_amount_std[ibin, 0, i] = (cmin+cmin+0.05)/2
#        cloud_amount_std[ibin, 1, i] = np.nanstd(FG_sub)
#        cloud_amount_std[ibin, 2, i] = count
#        cloud_amount_std[ibin, 3, i] = np.nanstd(FG_sub_qc)
#        cmin = cmin+0.05
#for i in range (4):
    # for high frequency
#    FG = GMI_high_final[:, i]-Sim_high_final[:,i+9,4]
#    stdk = np.nanstd(FG)
#    FG_qc = FG
#    FG_qc[FG_qc > 2*stdk] = np.nan
#   FG_qc[FG_qc < -2*stdk] = np.nan
#    cmin = -0.075
#    for ibin in range (22):
#        ind = np.where((CA_avg>cmin)*(CA_avg<=cmin+0.05))
#        FG_sub = FG[ind]
#        FG_sub_qc = FG_qc[ind]
#        count = np.count_nonzero(FG_sub)
#        cloud_amount_std[ibin, 0, i+9] = (cmin+cmin+0.05)/2
#        cloud_amount_std[ibin, 1, i+9] = np.nanstd(FG_sub)
#        cloud_amount_std[ibin, 2, i+9] = count
#        cloud_amount_std[ibin, 3, i+9] = np.nanstd(FG_sub_qc)
#        cmin = cmin+0.05
# calculate percentage
#total_samples = np.sum(cloud_amount_std[:, 2,2])
#cloud_per= (cloud_amount_std[:,2,2]*100)/total_samples
fig = plt.figure(figsize=(27, 13))
ax1= plt.subplot(2,4,1)
ax1.set_title('(a) 19V')
ax1.plot(cloud_amount_std_19V[:,0], cloud_amount_std_19V[:,1], 'k')
ax1.plot(cloud_amount_std_19V[:,0], cloud_amount_std_19V[:,3], 'k--')
ax1.set_xlabel('$ C_{37avg} $')
ax1.set_ylabel('Std (K)')
ax1.set_ylim([0, 40])
ax1.set_xlim([0, 1.0])
ax1.set_xticks(np.arange(0, 1.2, step=0.2))
plt.show()
error



ax2= plt.subplot(2,4,2)
ax2.set_title('(b) 23V')
ax2.plot(cloud_amount_std[:,0, 4], cloud_amount_std[:,1, 4], 'k')
ax2.plot(cloud_amount_std[:,0, 4], cloud_amount_std[:,3, 4], 'k--')
ax2.set_xlabel('$ C_{37avg} $')
ax2.set_ylabel('Std (K)')
ax2.set_ylim([0, 40])
ax2.set_xlim([0, 1.0])
ax2.set_xticks(np.arange(0, 1.2, step=0.2))

ax3= plt.subplot(2, 4,3)
ax3.set_title('(c) 37V')
ax3.plot(cloud_amount_std[:,0, 5], cloud_amount_std[:,1, 5], 'k')
ax3.plot(cloud_amount_std[:,0, 5], cloud_amount_std[:,3, 5], 'k--')
ax3.set_xlabel('$ C_{37avg} $')
ax3.set_ylabel('Std (K)')
ax3.set_ylim([0, 40])
ax3.set_xlim([0, 1.0])
ax3.set_xticks(np.arange(0, 1.2, step=0.2))

ax4= plt.subplot(2, 4,4)
ax4.set_title('(d) 89V')
ax4.plot(cloud_amount_std[:,0, 7], cloud_amount_std[:,1, 7], 'k')
ax4.plot(cloud_amount_std[:,0, 7], cloud_amount_std[:,3, 7], 'k--')
ax4.set_xlabel('$ C_{37avg} $')
ax4.set_ylabel('Std (K)')
ax4.set_ylim([0, 70])
ax4.set_xlim([0, 1.0])
ax4.set_xticks(np.arange(0, 1.2, step=0.2))

ax5= plt.subplot(2,4,5)
ax5.set_title('(e) 166V')
ax5.plot(cloud_amount_std[:,0, 9], cloud_amount_std[:,1, 9], 'k')
ax5.plot(cloud_amount_std[:,0, 9], cloud_amount_std[:,3, 9], 'k--')
ax5.set_xlabel('$ C_{37avg} $')
ax5.set_ylabel('Std (K)')
ax5.set_ylim([0, 80])
ax5.set_xlim([0, 1.0])
ax5.set_xticks(np.arange(0, 1.2, step=0.2))

ax6= plt.subplot(2, 4,6)
ax6.set_title( ' (f) $183\pm\ 3V$')
ax6.plot(cloud_amount_std[:,0, 11], cloud_amount_std[:,1, 11], 'k')
ax6.plot(cloud_amount_std[:,0, 11], cloud_amount_std[:,3, 11], 'k--')
ax6.set_xlabel('$ C_{37avg} $')
ax6.set_ylabel('Std (K)')
ax6.set_ylim([0, 70])
ax6.set_xlim([0, 1.0])
ax6.set_xticks(np.arange(0, 1.2, step=0.2))

ax7= plt.subplot(2,4,7)
ax7.set_title('(g) $183\pm\ 7V$')
ax7.plot(cloud_amount_std[:,0, 12], cloud_amount_std[:,1, 12], 'k')
ax7.plot(cloud_amount_std[:,0, 12], cloud_amount_std[:,3, 12], 'k--')
ax7.set_xlabel('$ C_{37avg} $')
ax7.set_ylabel('Std (K)')
ax7.set_ylim([0, 80])
ax7.set_xlim([0, 1.0])
ax7.set_xticks(np.arange(0, 1.2, step=0.2))

ax8= plt.subplot(2,4,8)
ax8.set_title('(h) Number of Samples')
ax8.plot(cloud_amount_std[:,0, 2], cloud_per, 'k')
ax8.set_xlabel('$ C_{37avg} $')
ax8.set_ylabel('Number of Samples (%)')
ax8.set_ylim([0, 60])
ax8.set_xlim([-0.1, 1.0])
ax8.set_xticks(np.arange(0, 1.2, step=0.2))

plt.savefig('D:/Python_processing/cyclone/results/error_model/std_dev_new.png', dpi=300, bbox_inches='tight')
plt.show()

error
