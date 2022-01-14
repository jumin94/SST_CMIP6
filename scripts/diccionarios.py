#Imports
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

from eofs.xarray import Eof
def indices(dato):
    dic = {}
    if 'year' in dato:
        dic['time'] = dato.year
    else:
        dic['time'] = dato.time        
    dic['nino34'] = dato.sel(lat=slice(5,-5)).sel(lon=slice(190,240)).mean(dim='lat').mean(dim='lon')
    dic['nino3'] = dato.sel(lat=slice(5,-5)).sel(lon=slice(210,270)).mean(dim='lat').mean(dim='lon')
    dic['nino4'] = dato.sel(lat=slice(5,-5)).sel(lon=slice(160,210)).mean(dim='lat').mean(dim='lon')
    dic['nino12'] = dato.sel(lat=slice(0,-10)).sel(lon=slice(270,280)).mean(dim='lat').mean(dim='lon')
    dic['west'] = dato.sel(lat=slice(5,-5)).sel(lon=slice(160,190)).mean(dim='lat').mean(dim='lon')
    dic['east'] = dato.sel(lat=slice(5,-5)).sel(lon=slice(240,270)).mean(dim='lat').mean(dim='lon')
    if 'year' in dato:
        a = dato.sst.rename({'year':'time'})
    else:
        a = dato.sst
    solver = Eof(a.sel(lat=slice(20,-20)).sel(lon=slice(150,280)))
    pcs = solver.pcs()
    dic['sst_mode1'] = pcs.sel(mode=0)
    dic['sst_mode2'] = pcs.sel(mode=1)
    dic['sst_mode3'] = pcs.sel(mode=2)
    dic['sst_mode4'] = pcs.sel(mode=3)
    dic['fv'] = solver.varianceFraction()
    dic['E_index'] = (dic['sst_mode1'] - dic['sst_mode2'])/np.sqrt(2)
    dic['C_index'] = (dic['sst_mode1'] + dic['sst_mode2'])/np.sqrt(2)
    return dic

