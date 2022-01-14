#Evaluo como simulan los modelos los modos de oscilacion del pacifico tropical
#imports
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature
from cartopy.util import add_cyclic_point
import matplotlib.path as mpath
import os
import glob
import pandas as pd
import xarray as xr
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth
import cartopy.util as cutil
import logging
import os, fnmatch
from scipy import signal 
import statsmodels.api as sm
from scipy import stats
import funciones 
import diccionarios
import clases

#Evaluo observaciones primero
#Open datasets
path_obs_sst = '/home/julia.mindlin'
ssts_kaplan = xr.open_dataset(path_obs_sst+'/Kaplan_SST/sst.mon.anom_remap.nc').sel(time=slice('1856-01','2020-12'))
ssts_ersstv = xr.open_dataset(path_obs_sst+'/ERSSTv5/ersstv5.sst.mnmean_remap.nc').sel(time=slice('1856-01','2020-12'))

anom_ssts_kaplan = funciones.anomalias(ssts_kaplan,ssts_kaplan).sel(time=slice('1900-01','2014-12'))  # Evaluo anomalias
anom_ssts_kaplan_djf = funciones.seasonal_data(anom_ssts_kaplan,'DJF')			# Selecciono la estacion de verano
dic_kaplan = diccionarios.indices(anom_ssts_kaplan_djf)					# Genero un diccionario de indices - cajas - eofs

#Armo la clase regresion para calcular mapas de regresion con los EOFs
reg = clases.regression()
regressors = pd.DataFrame({'mode1':funciones.standardize(dic_kaplan['sst_mode1']),'mode2':funciones.standardize(dic_kaplan['sst_mode2']),
                          'mode3':funciones.standardize(dic_kaplan['sst_mode3']),'mode4':funciones.standardize(dic_kaplan['sst_mode4'])})
reg.regressors = regressors                                                                #asigno el dataframe como regresores de la clase
out_reg_kaplan = reg.perform_regression(anom_ssts_kaplan_djf.sst)      #hago regresion
levels = np.arange(-.5,.55,.05)
fv = [round(dic_kaplan['fv'].sel(mode=0).values*100,2),round(dic_kaplan['fv'].sel(mode=1).values*100,2),
     round(dic_kaplan['fv'].sel(mode=2).values*100,2),round(dic_kaplan['fv'].sel(mode=3).values*100,2)]
dic_tit_mapa = ['Kaplan SSTs EOF1 '+str(fv[0])+'% 1900 - 2014 (DJF)','Kaplan SSTs EOF2 '+str(fv[1])+'% 1900 - 2014 (DJF)',
               'Kaplan SSTs EOF3 '+str(fv[2])+'% 1900 - 2014 (DJF)','Kaplan SSTs EOF4 '+str(fv[3])+'% 1900 - 2014 (DJF)']
dic_tit_serie = ['PC1','PC2','PC3','PC4']
datos = [-out_reg_kaplan['mode1']['coef'],out_reg_kaplan['mode2']['coef'],
        out_reg_kaplan['mode3']['coef'],out_reg_kaplan['mode4']['coef']]
datos_r2 = datos #[out_reg_kaplan['mode2']['r2'],out_reg_kaplan['mode2']['r2']]
t = [dic_kaplan['time'],dic_kaplan['time'],dic_kaplan['time'],dic_kaplan['time']]
series = [dic_kaplan['sst_mode1'].values,dic_kaplan['sst_mode2'].values,
          dic_kaplan['sst_mode3'].values,dic_kaplan['sst_mode4'].values]
levels = [np.arange(-1,1,.1),np.arange(-1,1,.1),np.arange(-1,1,.1),np.arange(-1,1,.1)]
fig = funciones.fig_sst_multiple(datos,datos_r2,t,series,dic_tit_mapa,dic_tit_serie,levels)
fig.savefig('/home/julia.mindlin/Tesis/Capitulo3/scripts/EOF_CMIP6_historical/figuras/Kaplan_EOF_analysis_1900_2014_DJF')

