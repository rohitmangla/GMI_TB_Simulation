import numpy as np
from pyhdf.SD import SD, SDC
import glob
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['axes.titlesize']  = 26
matplotlib.rcParams['legend.fontsize'] = 15.
matplotlib.rcParams['font.size']       = 28.
matplotlib.rcParams['axes.labelsize']  = 24.
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

# Median Vertical Profile Diagram
Med_ref_IGP_CbT  = np.zeros((36, 1), dtype = np.float64)
Med_ref_IGP_ICC8 = np.zeros((36, 1), dtype = np.float64)
Med_ref_IGP_ICC3   = np.zeros((36, 1), dtype = np.float64)

Med_ref_WG_CbT  = np.zeros((36, 1), dtype = np.float64)
Med_ref_WG_ICC8 = np.zeros((36, 1), dtype = np.float64)
Med_ref_WG_ICC3 = np.zeros((36, 1), dtype = np.float64)

ini_ref = 20

for i in range (36):

      x_IGP_CbT  = Maxref_IGP_CbT[:, i]
      x_IGP_ICC8 = Maxref_IGP_ICC8[:, i]
      x_IGP_ICC3 = Maxref_IGP_ICC3[:, i]

      x_WG_CbT = Maxref_WG_CbT[:, i]
      x_WG_ICC8 = Maxref_WG_ICC8[:, i]
      x_WG_ICC3 = Maxref_WG_ICC3[:, i]

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

      Med_ref_IGP_CbT[i] = np.median(x_IGP_CbT)
      Med_ref_IGP_ICC8[i] = np.median(x_IGP_ICC8)
      Med_ref_IGP_ICC3[i] = np.median(x_IGP_ICC3)

      Med_ref_WG_CbT[i]  = np.median(x_WG_CbT)
      Med_ref_WG_ICC8[i] = np.median(x_WG_ICC8)
      Med_ref_WG_ICC3[i] = np.median(x_WG_ICC3)


# Line Plots
height = np.linspace(0.5, 18, 36)

# Contour Plot
fig = plt.figure(figsize=(25, 6))
ax0  = plt.subplot(1, 3, 1)
ax0.plot(Med_ref_IGP_CbT, height, 'r',  linewidth=2, label='IGP')
ax0.plot(Med_ref_WG_CbT, height, 'b',  linewidth=2, label='WG')
ax0.set_title(' (a) PF-CbT')
ax0.set_xlabel('Median Maximum DPR Reflectivity (dBZ)')
ax0.set_ylabel('Height (km)')
ax0.set_xlim(20, 60)
ax0.set_ylim(1, 18)
ax0.set_xticks(np.arange(20, 60, step=10))
ax0.set_yticks(np.arange(1, 18, step=4))

ax0.legend()
plt.tight_layout()

ax1  = plt.subplot(1, 3, 2)
ax1.plot(Med_ref_IGP_ICC8, height, 'r', linewidth=2)
ax1.plot(Med_ref_WG_ICC8, height, 'b', linewidth=2)
ax1.set_title('(b) PF-ICC8')
ax1.set_xlabel('Median Maximum DPR Reflectivity (dBZ)')
ax1.set_ylabel('Height (km)')
ax1.set_ylim(1, 18)
ax1.set_xlim(20, 60)
ax1.set_yticks(np.arange(1, 18, step= 4))
ax1.set_xticks(np.arange(20, 60, step=10))

plt.tight_layout()

ax2  = plt.subplot(1, 3, 3)
ax2.plot(Med_ref_IGP_ICC3, height, 'r', linewidth=2)
ax2.plot(Med_ref_WG_ICC3, height, 'b', linewidth=2)
ax2.set_title('(c) PF-ICC3')
ax2.set_xlabel('Median Maximum DPR Reflectivity (dBZ)')
ax2.set_ylabel('Height (km)')
ax2.set_ylim(1, 18)
ax2.set_xlim(20, 60)
ax2.set_yticks(np.arange(1, 18, step=4))
ax2.set_xticks(np.arange(20, 60, step=10))
plt.tight_layout()

plt.savefig('D:/PhD/First_Paper/GPM_PF/Figure_3.tiff', dpi= 300, bbox_inches = 'tight')
plt.show()












