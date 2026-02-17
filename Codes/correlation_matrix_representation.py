import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib

matplotlib.rcParams['axes.titlesize'] = 16
matplotlib.rcParams['legend.fontsize'] =12
matplotlib.rcParams['font.size']= 18
matplotlib.rcParams['axes.labelsize']= 15
matplotlib.rcParams['xtick.labelsize'] = 7
matplotlib.rcParams['ytick.labelsize'] = 7
matplotlib.rcParams['figure.subplot.bottom'] = 0.1
matplotlib.rcParams['figure.subplot.top'] = 0.9
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
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

cor_mat_low  = np.zeros((channels_low, mie_tables), dtype = np.float64)
cor_mat_high = np.zeros((channels_high, mie_tables), dtype = np.float64)

for ichannel in range (channels_low):
    for itable in range (mie_tables):
        obs = GMI_low_final[:, ichannel]
        obs = obs[~np.isnan(obs)]
        sim = Sim_low_final[:, ichannel, itable]
        sim = sim[~np.isnan(sim)]
        coef = pearsonr(obs, sim)
        cor_mat_low[ichannel, itable] = coef[0]

for ichannel in range (channels_high):
    for itable in range (mie_tables):
        obs = GMI_high_final[:, ichannel]
        ind = np.argwhere(np.isnan(obs))
        obs = obs[~np.isnan(obs)]
        sim = Sim_high_final[:, ichannel+9, itable]
        sim = np.delete(sim, ind, 0)
        coef = pearsonr(obs, sim)
        cor_mat_high[ichannel, itable] = coef[0]
cor_mat_total = np.concatenate((cor_mat_low, cor_mat_high), 0)
cor_mat_filter = np.zeros((6, 26), dtype= np.float64)
cor_mat_filter[0, :] = cor_mat_total[2,:]
cor_mat_filter[1, :] = cor_mat_total[4,:]
cor_mat_filter[2, :] = cor_mat_total[5,:]
cor_mat_filter[3, :] = cor_mat_total[7,:]
cor_mat_filter[4, :] = cor_mat_total[9,:]
cor_mat_filter[5, :] = cor_mat_total[12,:]

x_axis_labels = ['long column', 'short column', 'block column','thick plate', 'thin plate', '3-bullet rosette', \
                 '4-bullet rosette', '5-bullet rosette', '6-bullet rosette', 'sector snowflake', 'dendrite snowflake',\
                 'PlateType1', 'ColumnType1', '6-BulletRosette', 'Perpendicular 4-Bullet Rosette', 'Flat3-Bulletrosette',\
                 'IconCloudIce', 'Sector Snowflake', 'EvansSnowAggregate', '8-columnAggregate','Large-PlateAggregate', \
                 'Large Column Aggregate', 'Large Block Aggregate', 'IconSnow', 'IconHail','GemGraupel'] # labels for x-axis

#y_axis_labels= ['10V', '10H', '19V', '19H', '24V', '37V', '37H', '89V', '89H', '166V', '166H', '183+3', '183+7']
y_axis_labels= ['19V', '23V', '37V', '89V', '166V', '$183\pm\ 7V$']

fig= plt.figure(figsize=(8,6))
ax = sns.heatmap(cor_mat_filter, vmin= 0.2, vmax=0.8, cmap=plt.cm.jet, \
                 xticklabels=x_axis_labels, yticklabels=y_axis_labels, square=True, cbar_kws={"shrink": .5})
ax.set_xticklabels(ax.get_xticklabels(),rotation=45,horizontalalignment='right')
plt.savefig('D:/Python_processing/cyclone/results/correlation_matrix.png', dpi=300, bbox_inches= 'tight')


#cmap=sns.diverging_coolwarm(220, 10, as_cmap=True)