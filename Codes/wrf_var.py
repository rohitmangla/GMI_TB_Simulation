import numpy as np



# wrf model variables
def wrf_var(ncfile):
   # extract variables
   lon       =  ncfile.variables['XLONG'   ][:]   # longitude
   lat       =  ncfile.variables['XLAT'    ][:]   # latitude
   lon_west  =  ncfile.variables['XLONG_U' ][:]   # longitude west
   lat_south =  ncfile.variables['XLAT_V'  ][:]   # latitude south
   T2        =  ncfile.variables['T2'      ][:]   # temperature at 2m
   tsk       =  ncfile.variables['TSK'     ][:]   # skin tempearture
   Q2        =  ncfile.variables['Q2'      ][:]   # specific humidity at 2m
   u10       =  ncfile.variables['U10'     ][:]   # u wind component at 10m
   v10       =  ncfile.variables['V10'     ][:]   # v wind component at 10m
   sfc       =  ncfile.variables['PSFC'    ][:]   # surface pressure
   lsm       =  ncfile.variables['LANDMASK'][:]   # landmask

   Q         =  ncfile.variables['QVAPOR'  ][:]   # water vapour or specific humiidty profile
   Clw       =  ncfile.variables['QCLOUD'  ][:]   # cloud liquid water  profile
   Ciw       =  ncfile.variables['QICE'    ][:]   # cloud ice  profile
   Snow      =  ncfile.variables['QSNOW'   ][:]   # snow  profile
   Rain      =  ncfile.variables['QRAIN'   ][:]   # RAIN  profile
   QGRP      =  ncfile.variables['QGRAUP'  ][:]   # Grapual profile
   Cc        =  ncfile.variables['CLDFRA'  ][:]   # Cloud Fraction
   times     =  ncfile.variables['Times'   ][:]   # time information
   hgt       =  ncfile.variables['HGT'     ][:]   # terrain height
   P         =  ncfile.variables['P'       ][:]   # perturbation pressure
   PB        =  ncfile.variables['PB'      ][:]   # Base state pressure
   # convert 3 to 2 dimensions

   lon        =  np.squeeze(lon)
   lat        =  np.squeeze(lat)
   lon_west   =  np.squeeze(lon_west)
   lat_south  =  np.squeeze(lat_south)
   T2         =  np.squeeze(T2)
   tsk        =  np.squeeze(tsk)
   Q2         =  np.squeeze(Q2)
   u10        =  np.squeeze(u10)
   v10        =  np.squeeze(v10)
   sfc        =  np.squeeze(sfc)
   lsm        =  np.squeeze(lsm)
   hgt        =  np.squeeze(hgt)
   
   Q          =  np.squeeze(Q)
   QGRP       =  np.squeeze(QGRP)
   Clw        =  np.squeeze(Clw)
   Ciw        =  np.squeeze(Ciw)
   Snow       =  np.squeeze(Snow)
   Rain       =  np.squeeze(Rain)
   Cc         =  np.squeeze(Cc)
   times      =  np.squeeze(times)
   P          =  np.squeeze(P)
   PB         =  np.squeeze(PB)
   
   nlevels   =   len(Rain)
   nrows     =   len(Rain[0])
   ncols     =   len(Rain[0][0])
   nprofiles =   nrows*ncols

   ## WRF variables  in RTTOV format

   P        =   np.transpose(np.reshape(P,    [nlevels, nprofiles]))
   P        =   np.fliplr(P)
   
   PB       =   np.transpose(np.reshape(PB,    [nlevels, nprofiles]))
   PB       =   np.fliplr(PB)
   
   Q        =   np.reshape(Q,    [nlevels, nprofiles])
   Q        =   np.transpose(Q)
   Q        =   np.fliplr(Q)

   QGRP     =   np.reshape(QGRP,    [nlevels, nprofiles])
   QGRP     =   np.transpose(QGRP)
   QGRP     =   np.fliplr(QGRP)
   
   Clw      =    np.reshape(Clw,  [nlevels, nprofiles])
   Clw      =    np.transpose(Clw)
   Clw      =   np.fliplr(Clw)

   Ciw      =    np.reshape(Ciw,  [nlevels, nprofiles])
   Ciw      =    np.transpose(Ciw)
   Ciw      =    np.fliplr(Ciw)
   
   Snow     =    np.reshape(Snow, [nlevels, nprofiles])
   Snow     =    np.transpose(Snow)
   Snow     =    np.fliplr(Snow)
   
   Rain     =    np.reshape(Rain, [nlevels, nprofiles])
   Rain     =    np.transpose(Rain)
   Rain     =    np.fliplr(Rain)
    
   Cc       =    np.reshape(Cc,   [nlevels, nprofiles])
   Cc       =    np.transpose(Cc)
   Cc       =    np.fliplr(Cc)
   
   lat      =    np.transpose(np.reshape(lat,  nprofiles))
   lon      =    np.transpose(np.reshape(lon,  nprofiles))
   sfc      =    np.transpose(np.reshape(sfc,  nprofiles))
   sfc      =    sfc/100
   T2       =    np.transpose(np.reshape(T2,   nprofiles))          
   Q2       =    np.transpose(np.reshape(Q2,   nprofiles))              
   u10      =    np.transpose(np.reshape(u10,  nprofiles))          
   v10      =    np.transpose(np.reshape(v10,  nprofiles))          
   tsk      =    np.transpose(np.reshape(tsk,  nprofiles))          
   lsm      =    np.transpose(np.reshape(lsm,  nprofiles))
   hgt      =    np.transpose(np.reshape(hgt,  nprofiles))
   hgt      =    hgt/1000
   
   return lon,lon_west, lat,lat_south, nprofiles, nlevels,Q, Clw, Ciw, Snow, Rain,QGRP, Cc, \
          sfc, T2, Q2, u10, v10,tsk,lsm, hgt, P, PB, nrows, ncols


  
