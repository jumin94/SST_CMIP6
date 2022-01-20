#Evaluo como simulan los modelos los modos de oscilacion del pacifico tropical
#imports
import cartopy.crs as ccrs
import netCDF4
import numpy as np
import pandas as pd
import numpy as np
import cartopy.feature
from cartopy.util import add_cyclic_point
import matplotlib.path as mpath
import os
import glob
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

#Evaluo simulaciones historicas
#Tres apporaches
#1. 
#Open datasets


#Funciones--------------------------------------------------------
def cargo_todo_crudos_remap(scenarios,models,ruta,var):
    os.chdir(ruta)
    os.getcwd()
    dic = {}
    dic['historical'] = {}
    dic['ssp585'] = {}
    for scenario in dic.keys():
        listOfFiles = os.listdir(ruta+'/'+scenario+'/'+var)
        for model in models:
            dic[scenario][model] = []
            pattern = "*"+model+"*"+scenario+"*remap*"
            for entry in listOfFiles:
                if fnmatch.fnmatch(entry,pattern):
                    print(pattern)
                    dato = xr.open_dataset(ruta+'/'+scenario+'/'+var+'/'+entry)
                    if scenario == "historical":
                        dato = dato#.sel(time=slice('1900','1999'))
                        dic[scenario][model].append(dato)
                    else:
                        dato = dato#.sel(time=slice('2014','2099'))
                        dic[scenario][model].append(dato)
    return dic

#Abro datos---------------------------
models = [
            'ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0',
            'CNRM-CM6-1', 'CESM2_','CESM2-WACCM', 'CanESM5', 'CMCC-CM2-SR5',
            'CNRM-ESM2-1','EC-Earth3', 'FGOALS-g3', 'GFDL-ESM4','HadGEM3-GC31-LL','HadGEM3-GC31-MM',
            'IITM-ESM','INM-CM4-8','IPSL-CM6A-LR','INM-CM5-0','KACE-1-0-G',
            'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR',
            'MRI-ESM2-0', 'NESM3', 'NorESM2-LM', 'NorESM2-MM','TaiESM1','UKESM1-0-LL'
            ]

scenarios = ['historical','ssp585']
variables = ['tos']

path = '/pikachu/datos/CMIP6_backup/Data_used/'
ruta = '/datos/julia.mindlin/CMIP6_remap'

var = ['/pr','/tos','/psl']
scenarios = ['historical','ssp585']
dato_sst = cargo_todo_crudos_remap(scenarios,models,ruta,var[1])



from eofs.xarray import Eof
#Hago analisis de componentes principales para DJF y guardo analisis en un .csv
path = '/home/julia.mindlin/Tesis/Capitulo3/scripts/EOF_SST_evaluation/figuras'

dic = {}
for model in models[:]:
    dic[model] = {}
    full_period = xr.merge([dato_sst[scenarios[0]][model][0].sel(time=slice('1950','2014')).tos,dato_sst[scenarios[1]][model][0].sel(time=slice('2015','2099')).tos])
    full_period = funciones.anomalias(full_period,full_period.sel(time=slice('1950','2000')))
    gw = full_period.mean(dim='lon').mean(dim='lat')
    reg = clases.regression()
    regressors = pd.DataFrame({'gw':gw.tos})
    reg.regressors = regressors
    aux = full_period.tos.rename({'time':'year'})
    out_reg1 = reg.perform_regression(aux)
    ssts_wo_gw = full_period - out_reg1['gw']['coef']*gw.tos
    
    dat = ssts_wo_gw.sel(lat=slice(15,-15)).sel(lon=slice(150,280))#.sel(time=slice('1950-01','1999-12'))
    lat = dat.lat; lon = dat.lon; time = dat.time
    solver = Eof(dat.tos)
    pcs = solver.pcs()
    eofs = solver.eofs()
    a = 0; b = 0; c = 0; d = 0
    if eofs.sel(mode=0).sel(lat=slice(5,-5)).sel(lon=slice(190,240)).mean(dim='lat').mean(dim='lon') > 0.:
        dic[model]['sst_mode1'] = pcs.sel(mode=0)
    else:
        dic[model]['sst_mode1'] = -pcs.sel(mode=0)
        a = 1
    if eofs.sel(mode=1).sel(lat=slice(5,-5)).sel(lon=slice(190,240)).mean(dim='lat').mean(dim='lon') > 0.:
        dic[model]['sst_mode2'] = pcs.sel(mode=1)
    else:
        dic[model]['sst_mode2'] = -pcs.sel(mode=1)
        b = 1
    if eofs.sel(mode=2).sel(lat=slice(5,-5)).sel(lon=slice(210,270)).mean(dim='lat').mean(dim='lon') > 0.:
        dic[model]['sst_mode3'] = pcs.sel(mode=2)
    else:
        dic[model]['sst_mode3'] = -pcs.sel(mode=2)
        c = 1
    if eofs.sel(mode=3).sel(lat=slice(5,-5)).sel(lon=slice(190,240)).mean(dim='lat').mean(dim='lon') > 0.:
        dic[model]['sst_mode4'] = pcs.sel(mode=3)
    else:
        dic[model]['sst_mode4'] = -pcs.sel(mode=3)
        d = 1

    dic[model]['E_index'] = (dic[model]['sst_mode1'] - dic[model]['sst_mode2'])/np.sqrt(2)
    dic[model]['C_index'] = (dic[model]['sst_mode1'] + dic[model]['sst_mode2'])/np.sqrt(2)
    dic[model]['nino34'] = dat.sel(lat=slice(5,-5)).sel(lon=slice(190,240)).mean(dim='lat').mean(dim='lon')
    dic[model]['nino3'] = dat.sel(lat=slice(5,-5)).sel(lon=slice(210,270)).mean(dim='lat').mean(dim='lon')
    dic[model]['nino4'] = dat.sel(lat=slice(5,-5)).sel(lon=slice(160,210)).mean(dim='lat').mean(dim='lon')
    dic[model]['nino12'] = dat.sel(lat=slice(0,-10)).sel(lon=slice(270,280)).mean(dim='lat').mean(dim='lon')
    dic[model]['west'] = dat.sel(lat=slice(5,-5)).sel(lon=slice(160,190)).mean(dim='lat').mean(dim='lon')
    dic[model]['fv'] = solver.varianceFraction()
    reg = clases.regression()
    regressors = pd.DataFrame({'mode1':funciones.standardize(dic[model]['sst_mode1']),'mode2':funciones.standardize(dic[model]['sst_mode2']),
                              'mode3':funciones.standardize(dic[model]['sst_mode3']),'mode4':funciones.standardize(dic[model]['sst_mode4'])})
    reg.regressors = regressors
    out_reg = reg.perform_regression(aux)
    levels = np.arange(-.5,.55,.05)

    #aux = out_reg_kaplan['mode2']['coef'].sel(lat=slice(15,-15)).sel(lon=slice(150,280)).fillna(0)
    #dic[model]['corr_kaplan'] = pattern_corr(eofs.sel(mode=1).fillna(0),aux)
    datos = []
    if a == 0:
        datos.append(out_reg1['gw']['coef'])
    else:
        datos.append(out_reg1['gw']['coef'])

    if b == 0:
        datos.append(out_reg['mode1']['coef'])
    else:
        datos.append(out_reg['mode1']['coef'])

    if c == 0:
        datos.append(out_reg['mode2']['coef'])
    else:
        datos.append(out_reg['mode2']['coef'])

    if d == 0:
        datos.append(out_reg['mode3']['coef'])
    else:
        datos.append(out_reg['mode3']['coef'])

    fv = [round(dic[model]['fv'].sel(mode=0).values*100,2),round(dic[model]['fv'].sel(mode=1).values*100,2),
          round(dic[model]['fv'].sel(mode=2).values*100,2),round(dic[model]['fv'].sel(mode=3).values*100,2)]
    dic_tit_mapa = [model+' GW 1950 - 2099 (DJF)',
                    model+' EOF1 '+str(fv[0])+'% 1950 - 2099 (DJF)',
                    model+' EOF2 '+str(fv[1])+'% 1950 - 2099 (DJF)',
                    model+' EOF3 '+str(fv[2])+'% 1950 - 2099 (DJF)']
    dic_tit_serie = ['GW','PC1','PC2','PC3']
    datos_r2 = datos
    #datos = [out_reg_ersstv5['mode1']['coef'],out_reg_ersstv5['mode2']['coef']]
    #datos_r2 = [out_reg_ersstv5['mode1']['r2'],out_reg_ersstv5['mode1']['r2']]
    t = [time,time,time,time]
    series = [gw.tos,dic[model]['sst_mode1'],dic[model]['sst_mode2'],dic[model]['sst_mode3']]
    levels = [np.arange(-3,3,.5),np.arange(-2,2,.5),np.arange(-.5,.6,.1),np.arange(-.5,.6,.1)]
    fig = funciones.fig_sst_multiple2(datos,datos_r2,t,series,dic_tit_mapa,dic_tit_serie,levels)
    fig.savefig(path+'/EOF_regress_out_GW_'+model+'_1950_2099.png') 


