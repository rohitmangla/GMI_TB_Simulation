import numpy as np
from pyhdf.SD import SD, SDC
import os
import glob
os.environ['PROJ_LIB'] = r'C:\Users\Rohit\Anaconda3\pkgs\proj4-5.2.0-h6538335_1006\Library\share'
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['axes.titlesize']  = 25
matplotlib.rcParams['legend.fontsize'] = 15.
matplotlib.rcParams['font.size']       = 26.
matplotlib.rcParams['axes.labelsize']  = 26.
matplotlib.rcParams['xtick.labelsize'] = 28.
matplotlib.rcParams['ytick.labelsize'] = 28.

PATH = 'D:/PhD/First_Paper/GPM_PF/'
FILES = glob.glob(PATH+"*.hdf")

Maxref_final_IGP= np.zeros((1,40), dtype = np.float64)
Maxref_final_WG= np.zeros((1,40), dtype = np.float64)

for ifile in FILES:
    print (ifile)
    file= SD (ifile, SDC.READ)
    # Maximum Reflectivtity Profile
    sds_dBZ = file.select('MAXDBZ')
    Maxref  = sds_dBZ.get()
    # Latitude
    lat    = file.select('LAT')
    LAT    = lat.get()
    # Longitude
    lon    = file.select('LON')
    LON    = lon.get()
    # Number of pixels
    Npixel = file.select('NPIXELS')
    NPIXEL = Npixel.get()

    # IGP Subset
    ind_IGP    = np.where((LAT>=23)*(LAT<=30)*(LON>=73)*(LON<=85))
    Lat_IGP    = LAT[ind_IGP]
    Lon_IGP    = LON[ind_IGP]
    Maxref_IGP = np.squeeze(Maxref[ind_IGP, :])
    Maxref_IGP[Maxref_IGP<0] = 0
    NPIXEL_IGP = NPIXEL[ind_IGP]

    # WG Subset
    ind_WG  = np.where((LAT >= 10) * (LAT <= 18) * (LON >= 72) * (LON <= 77))
    Lat_WG = LAT[ind_WG]
    Lon_WG = LON[ind_WG]
    Maxref_WG = np.squeeze(Maxref[ind_WG, :])
    Maxref_WG[Maxref_WG < 0] = 0
    NPIXEL_WG = NPIXEL[ind_WG]
    # Area greater than 100 km2(pixel greater than 4)
    ind_npixel_IGP = np.where(NPIXEL_IGP >= 4)
    ind_npixel_WG  = np.where(NPIXEL_WG >= 4)

    Lat_IGP        = Lat_IGP[ind_npixel_IGP]
    Lon_IGP        = Lon_IGP[ind_npixel_IGP]
    Maxref_IGP     =  np.squeeze(Maxref_IGP[ind_npixel_IGP, :])

    Lat_WG = Lat_WG[ind_npixel_WG]
    Lon_WG = Lon_WG[ind_npixel_WG]
    Maxref_WG = np.squeeze(Maxref_WG[ind_npixel_WG, :])
    # Extend Arrays

    Maxref_final_IGP = np.append(Maxref_final_IGP, Maxref_IGP, axis=0)
    Maxref_final_WG = np.append(Maxref_final_WG, Maxref_WG, axis=0)
    error

# Remove first row
Maxref_final_IGP = Maxref_final_IGP[1:, :]
Maxref_final_WG  = Maxref_final_WG[1:, :]

# Condition for CbT Clouds [Keep vertical layer only those with h>12 km and Z>20 dBZ]
ind_CbT_IGP = np.where(Maxref_final_IGP[:, 24]>=20)
ind_CbT_WG  = np.where(Maxref_final_WG[:, 24]>=20)

Maxref_IGP_CbT = np.squeeze(Maxref_final_IGP[ind_CbT_IGP, :])
Maxref_WG_CbT  = np.squeeze(Maxref_final_WG[ind_CbT_WG, :])

# Condition for ICC8 Clouds [Keep vertical layer only those with h>8 km and Z>32 dBZ]
ind_ICC8_IGP = np.where(Maxref_final_IGP[:, 16]>=32)
ind_ICC8_WG  = np.where(Maxref_final_WG[:, 16]>=31)

Maxref_IGP_ICC8 = np.squeeze(Maxref_final_IGP[ind_ICC8_IGP, :])
Maxref_WG_ICC8 = np.squeeze(Maxref_final_WG[ind_ICC8_WG, :])

# Condition for ICC3 clouds
ind_ICC3_IGP = np.where(Maxref_final_IGP[:, 6]>=41)
ind_ICC3_WG = np.where(Maxref_final_WG[:, 6]>=40)

Maxref_IGP_ICC3 = np.squeeze(Maxref_final_IGP[ind_ICC3_IGP, :])
Maxref_WG_ICC3 = np.squeeze(Maxref_final_WG[ind_ICC3_WG, :])

# Probability Diagram
count_IGP_CbT  = np.zeros((70, 34), dtype = np.float64)
count_IGP_ICC8 = np.zeros((70, 34), dtype = np.float64)
count_IGP_ICC3 = np.zeros((70, 34), dtype = np.float64)

count_WG_CbT  = np.zeros((70, 34), dtype = np.float64)
count_WG_ICC8 = np.zeros((70, 34), dtype = np.float64)
count_WG_ICC3 = np.zeros((70, 34), dtype = np.float64)

Ref_IGP_CbT_pr  = np.zeros((70, 34), dtype= np.float64)
Ref_IGP_ICC8_pr = np.zeros((70, 34), dtype= np.float64)
Ref_IGP_ICC3_pr = np.zeros((70, 34), dtype= np.float64)

Ref_WG_CbT_pr  = np.zeros((70, 34), dtype= np.float64)
Ref_WG_ICC8_pr = np.zeros((70, 34), dtype= np.float64)
Ref_WG_ICC3_pr = np.zeros((70, 34), dtype= np.float64)

ini_ref = 20

for i in range (70):
    for j in range (34):
        x_IGP_CbT  = Maxref_IGP_CbT[:, j]
        x_IGP_ICC8 = Maxref_IGP_ICC8[:, j]
        x_IGP_ICC3 = Maxref_IGP_ICC3[:, j]

        x_WG_CbT = Maxref_WG_CbT[:, j]
        x_WG_ICC8 = Maxref_WG_ICC8[:, j]
        x_WG_ICC3 = Maxref_WG_ICC3[:, j]

        x_IGP_CbT[x_IGP_CbT<18]   = np.nan
        x_IGP_ICC8[x_IGP_ICC8<18] = np.nan
        x_IGP_ICC3[x_IGP_ICC3<18] = np.nan

        x_WG_CbT[x_WG_CbT < 18]   = np.nan
        x_WG_ICC8[x_WG_ICC8 < 18] = np.nan
        x_WG_ICC3[x_WG_ICC3 < 18] = np.nan

        x_IGP_CbT = x_IGP_CbT[~np.isnan(x_IGP_CbT)]
        x_IGP_ICC8 = x_IGP_ICC8[~np.isnan(x_IGP_ICC8)]
        x_IGP_ICC3 = x_IGP_ICC3[~np.isnan(x_IGP_ICC3)]

        x_WG_CbT = x_WG_CbT[~np.isnan(x_WG_CbT)]
        x_WG_ICC8 = x_WG_ICC8[~np.isnan(x_WG_ICC8)]
        x_WG_ICC3 = x_WG_ICC3[~np.isnan(x_WG_ICC3)]

        y_IGP_CbT = np.where((x_IGP_CbT >= ini_ref)*(x_IGP_CbT <ini_ref + 1))
        y_IGP_ICC8 = np.where((x_IGP_ICC8 >= ini_ref)*(x_IGP_ICC8 <ini_ref + 1))
        y_IGP_ICC3 = np.where((x_IGP_ICC3 >= ini_ref)*(x_IGP_ICC3 <ini_ref + 1))

        y_WG_CbT = np.where((x_WG_CbT >= ini_ref) * (x_WG_CbT < ini_ref + 1))
        y_WG_ICC8 = np.where((x_WG_ICC8 >= ini_ref) * (x_WG_ICC8 < ini_ref + 1))
        y_WG_ICC3 = np.where((x_WG_ICC3 >= ini_ref) * (x_WG_ICC3 < ini_ref + 1))

        count_IGP_CbT[i, j] = np.size(y_IGP_CbT)
        count_IGP_ICC8[i, j] = np.size(y_IGP_ICC8)
        count_IGP_ICC3[i, j] = np.size(y_IGP_ICC3)

        count_WG_CbT[i, j] = np.size(y_WG_CbT)
        count_WG_ICC8[i, j] = np.size(y_WG_ICC8)
        count_WG_ICC3[i, j] = np.size(y_WG_ICC3)

    ini_ref= ini_ref+1

total_count_IGP_CbT= np.max(count_IGP_CbT)
total_count_IGP_ICC8= np.max(count_IGP_ICC8)
total_count_IGP_ICC3= np.max(count_IGP_ICC3)

total_count_WG_CbT = np.max(count_WG_CbT)
total_count_WG_ICC8= np.max(count_WG_ICC8)
total_count_WG_ICC3= np.max(count_WG_ICC3)

# Normalized the occurrences
for k in range (70):
    for l in range (34):
        Ref_IGP_CbT_pr[k,l] = count_IGP_CbT[k,l]/total_count_IGP_CbT
        Ref_IGP_ICC8_pr[k,l] = count_IGP_ICC8[k,l]/total_count_IGP_ICC8
        Ref_IGP_ICC3_pr[k,l] = count_IGP_ICC3[k,l]/total_count_IGP_ICC3

        Ref_WG_CbT_pr[k, l]   = count_WG_CbT[k, l] / total_count_WG_CbT
        Ref_WG_ICC8_pr[k, l]  = count_WG_ICC8[k, l] / total_count_WG_ICC8
        Ref_WG_ICC3_pr[k, l] = count_WG_ICC3[k, l] / total_count_WG_ICC3


Ref_IGP_CbT_pr[Ref_IGP_CbT_pr<0.1] =np.nan
Ref_IGP_ICC8_pr[Ref_IGP_ICC8_pr<0.1] =np.nan
Ref_IGP_ICC3_pr[Ref_IGP_ICC3_pr<0.1] =np.nan

Ref_WG_CbT_pr[Ref_WG_CbT_pr<0.1] =np.nan
Ref_WG_ICC8_pr[Ref_WG_ICC8_pr<0.1] =np.nan
Ref_WG_ICC3_pr[Ref_WG_ICC3_pr<0.1] =np.nan

# Spatial Distribution plot
cmap  = plt.cm.jet
pr_min = 0
pr_max = 1
res_colorbar =  0.2
clevs = np.arange(0, (pr_max-pr_min)/res_colorbar+1)*res_colorbar+pr_min
# Create Spatial Grids
xi = np.linspace(21, 90, 70)
yi = np.linspace(0.5, 17, 34)
xx, yy = np.meshgrid(yi, xi)
# Contour Plot
fig = plt.figure(figsize=(25, 16))
ax0  = plt.subplot(2, 3, 1)
cs1 = plt.contourf(yy, xx, Ref_IGP_CbT_pr, 15, levels = clevs, cmap= cmap)
ax0.set_title(' (a) PF-CbT_IGP')
ax0.set_xlabel('Maximum Reflectivity (dBZ)')
ax0.set_ylabel('Height (km)')
ax0.set_ylim(0, 18)
ax0.set_xlim(21, 60)
plt.tight_layout()

ax1  = plt.subplot(2, 3, 2)
cs2 = plt.contourf(yy, xx, Ref_IGP_ICC8_pr, 15, levels = clevs, cmap= cmap)
ax1.set_title('(b) PF-ICC8_IGP')
ax1.set_xlabel('Maximum Reflectivity (dBZ)')
ax1.set_ylabel('Height (km)')
ax1.set_ylim(0, 18)
ax1.set_xlim(21, 60)
plt.tight_layout()

ax2  = plt.subplot(2, 3, 3)
cs3 = plt.contourf(yy, xx, Ref_IGP_ICC3_pr, 15, levels = clevs, cmap= cmap)
ax2.set_title(' (c) PF-ICC3_IGP')
ax2.set_xlabel('Maximum Reflectivity (dBZ)')
ax2.set_ylabel('Height (km)')
ax2.set_ylim(0, 18)
ax2.set_xlim(21, 60)
plt.tight_layout()

ax3  = plt.subplot(2, 3, 4)
cs4 = plt.contourf(yy, xx, Ref_WG_CbT_pr, 15, levels = clevs, cmap= cmap)
ax3.set_title('(d) PF-CbT_WG')
ax3.set_xlabel('Maximum Reflectivity (dBZ)')
ax3.set_ylabel('Height (km)')
ax3.set_ylim(0, 18)
ax3.set_xlim(21, 60)
plt.tight_layout()

ax4  = plt.subplot(2, 3, 5)
cs5 = plt.contourf(yy, xx, Ref_WG_ICC8_pr, 15, levels = clevs, cmap= cmap)
ax4.set_title('(e) PF-ICC8_WG')
ax4.set_xlabel('Maximum Reflectivity (dBZ)')
ax4.set_ylabel('Height (km)')
ax4.set_ylim(0, 18)
ax4.set_xlim(21, 60)
plt.tight_layout()

ax5  = plt.subplot(2, 3, 6)
cs6 = plt.contourf(yy, xx, Ref_WG_ICC3_pr, 15, levels = clevs, cmap= cmap)
ax5.set_title('(f) PF-ICC3_WG')
ax5.set_xlabel('Maximum Reflectivity (dBZ)')
ax5.set_ylabel('Height (km)')
ax5.set_ylim(0, 18)
ax5.set_xlim(21, 60)
plt.tight_layout()

cb_ax = fig.add_axes([1, 0.1, 0.02, 0.8])
cbar = plt.colorbar(cs1,  cax=cb_ax)
plt.savefig('D:/PhD/First_Paper/GPM_PF/Figure_2.tiff', dpi= 300, bbox_inches = 'tight')
plt.show()












