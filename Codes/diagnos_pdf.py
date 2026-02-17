import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib
matplotlib.rcParams['axes.titlesize'] = 18
matplotlib.rcParams['legend.fontsize'] = 14.
matplotlib.rcParams['font.size'] = 22.
#matplotlib.rcParams['axes.labelsize'] = 16.
matplotlib.rcParams['xtick.labelsize'] = 17.
matplotlib.rcParams['ytick.labelsize'] = 17.
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.wspace'] = 0.3

# import datasets
PATH = "D:/Python_processing/cyclone/forecast_12hrs/"
FILES= glob.glob(PATH+"*.nc")

channels_low =9
channels_high= 4
channels = 13
mie_tables= 26
GMI_high_final = np.zeros((1, channels_high), dtype=np.float64)
GMI_low_final  = np.zeros((1, channels_low), dtype=np.float64)
Sim_high_final = np.zeros((1, channels, mie_tables), dtype=np.float64)
Sim_low_final  = np.zeros((1, channels, mie_tables), dtype=np.float64)

for ifile in FILES:
    print(ifile)
    ncfile = Dataset(ifile,'r')
    lat = ncfile.variables['lat'][:]    # latitude
    lon = ncfile.variables['lon'][:]    # longitude
    Simulated_Tb_high = ncfile.variables['Simulated_Tb_high'][:]  # high frequency Simulation
    Simulated_Tb_low  = ncfile.variables['Simulated_Tb_low'][:]  # low frequency Simulation
    GMI_Tb_high       = ncfile.variables['GMI_gridded_high'][:]  # high frequency GMI Tb
    GMI_Tb_low        = ncfile.variables['GMI_gridded_low'][:]  # high frequency GMI Tb
    lsm               = ncfile.variables['lsm'][:]  # land surface mask 0 for ocean, 1 for land
    # Land masking
    ind_land = np.where(lsm == 1)
    Simulated_Tb_high[ind_land, :, :] = np.nan
    Simulated_Tb_low[ind_land, :, :] = np.nan
    GMI_Tb_high[ind_land, :] = np.nan
    GMI_Tb_low[ind_land, :] = np.nan
    ind_high = np.where(GMI_Tb_high[:, 0] > 1)
    lat_high = lat[ind_high]
    lon_high = lon[ind_high]
    ind_low = np.where(GMI_Tb_low[:, 0] > 1)
    lat_low = lat[ind_low]
    lon_low = lon[ind_low]

    GMI_high = np.squeeze(GMI_Tb_high[ind_high, :])  # (profiles, channels_high)
    GMI_low = np.squeeze(GMI_Tb_low[ind_low, :])  # (profiles, channels_low)
    GMI_high_final = np.concatenate((GMI_high_final, GMI_high), 0)
    GMI_low_final = np.concatenate((GMI_low_final, GMI_low), 0)
    Sim_high = np.squeeze(Simulated_Tb_high[ind_high, :, :])  # (profiles, channels, mietables)
    Sim_high_final = np.concatenate((Sim_high_final, Sim_high), 0)
    Sim_low  = np.squeeze(Simulated_Tb_low[ind_low, :, :])  # (profiles, channels, mietables)
    Sim_low_final = np.concatenate((Sim_low_final, Sim_low), 0)

GMI_high_final[0, :] = np.nan
Sim_high_final[0, :, :] = np.nan
GMI_low_final[0, :] = np.nan
Sim_low_final[0, :, :] = np.nan
n= np.where((GMI_low_final[:, 5]>160)*(GMI_low_final[:, 5]<200))
GMI_low_final[n,5]= np.nan
Sim_low_final[n,5,:]= np.nan

nbins = 80
delta_bin= 2.5
obs_high = np.zeros((nbins,channels_high), dtype= np.float64)
sim_high = np.zeros((nbins,channels_high, mie_tables), dtype= np.float64)
obs_low = np.zeros((nbins,channels_low), dtype= np.float64)
sim_low = np.zeros((nbins,channels_low, mie_tables), dtype= np.float64)

# for high range frequencies:
bin_start = 100
for ibin in range (nbins):
    for ichannel in range (channels_high):
        for itable in range (mie_tables):
            ind_obs= np.where((GMI_high_final[:, ichannel]>= bin_start)*(GMI_high_final[:, ichannel]<bin_start+delta_bin))
            ind_sim= np.where((Sim_high_final[:, ichannel+9, itable]>= bin_start)*(Sim_high_final[:, ichannel+9, itable]<bin_start+delta_bin))
            obs = GMI_high_final[:, ichannel][ind_obs]
            obs_count = np.count_nonzero(obs)
            obs_high[ibin, ichannel] = obs_count
            sim = Sim_high_final[:, ichannel+9, itable][ind_sim]
            sim_count = np.count_nonzero(sim)
            sim_high[ibin,ichannel,itable] = sim_count
    bin_start= bin_start+delta_bin

# for low range frequencies:
bin_start = 100
for ibin in range (nbins):
    for ichannel in range (channels_low):
        for itable in range (mie_tables):
            ind_obs= np.where((GMI_low_final[:, ichannel]>= bin_start)*(GMI_low_final[:, ichannel]<bin_start+delta_bin))
            ind_sim= np.where((Sim_low_final[:, ichannel, itable]>= bin_start)*(Sim_low_final[:, ichannel, itable]<bin_start+delta_bin))
            obs = GMI_low_final[:, ichannel][ind_obs]
            obs_count = np.count_nonzero(obs)
            obs_low[ibin, ichannel] = obs_count
            sim = Sim_low_final[:, ichannel, itable][ind_sim]
            sim_count = np.count_nonzero(sim)
            sim_low[ibin,ichannel,itable] = sim_count
    bin_start= bin_start+delta_bin

# condider one location [20,8]:
k=100.0
l= 0.0
for ichannel in range (1):   # channels
    ichannel= 3
    for ibin in range (1):  # brightness temp
        ibin=72
        z_sim = sim_high[ibin, ichannel, :]
        z_obs = obs_high[ibin,ichannel   ]
        #if (max(z_sim)<4000) | (max(z_sim)>20000):continue
        q1 =  np.percentile(z_sim, 5)
        q2 =  np.percentile(z_sim, 25)
        q3 =  np.percentile(z_sim, 50)
        q4 =  np.percentile(z_sim, 75)
        q5 =  np.percentile(z_sim, 95)
        fig  = plt.figure(figsize=(12,4))
        #plt.hist (z_sim,20, range = (4000, 20000), histtype= 'step', facecolor= 'b', alpha= 0.5, linewidth = 2)
        plt.hist (z_sim, 20, range= (5500, 7000), histtype= 'step', facecolor= 'b', alpha= 0.5, linewidth = 2)
        plt.axvline(x=q1,ymin = 0, ymax= 3,  linewidth = 2, color= 'y',  linestyle=  'dashed')
        plt.axvline(x=q2,ymin = 0, ymax= 3,  linewidth = 2, color= 'c',  linestyle = 'dashed')
        plt.axvline(x=q3, ymin = 0, ymax= 3,  linewidth = 2, color= 'r', linestyle = 'dashed')
        plt.axvline(x=q4, ymin = 0, ymax= 3, linewidth = 2, color= 'm',  linestyle = 'dashed')
        plt.axvline(x=q5, ymin= 0 , ymax= 3, linewidth = 2, color= 'g',  linestyle = 'dashed')
        plt.annotate('[OBS]', (z_obs, 0), xytext=(6500, 8),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=14, horizontalalignment='right', verticalalignment='bottom')
        plt.annotate('(5%)',  (q1,3), fontsize=13, horizontalalignment='right')
        plt.annotate('(25%)', (q2,5), fontsize=13, horizontalalignment='right')
        plt.annotate('(50%)', (q3,7), fontsize=13, horizontalalignment='right')
        plt.annotate('(75%)', (q4,5), fontsize=13, horizontalalignment='right')
        plt.annotate('(95%)', (q5,9), fontsize=13, horizontalalignment='right')
        plt.ylim(0,23)
        plt.xlim(6300, 6900)
        plt.yticks(np.arange(0, 24, step=4))
        #plt.xticks(np.arange(0, 6000, step=1000))
        #plt.set_xticklabels(['', '0', '200', '400', '600', '800', '1000'])
        plt.xlabel('Number of WRF/OBS samples', fontsize = 15)
        plt.ylabel ('Number of RTTOVSCATT\n Configurations', fontsize = 15)
        plt.tight_layout()
        plt.savefig('D:/Python_processing/cyclone/results/pdf/pdf'+str(ibin)+'_Tb'+str(89)+'_ch'+'.png', dpi=300, bbox_inches = 'tight')
        k= k+delta_bin
        plt.close()
    l= l+1
    k= 100.0

