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

import statsmodels.api as sm
from sklearn import linear_model
class regression(object):
    from sklearn import linear_model
    def __init__(self):
        self.what_is_this = 'This performs a regression patterns for principal component analysis'
    
    def create_x(self,i,j,dato):
        x = dato[:,i,j].values
        return x
    
    def regressors(self,regressors):
        """Recibe un DataFrame con las series para hacer regresion"""
                #por ejemplo:
                #regressors = pd.DataFrame({'pc1':self.stand(pcs['pc1']),
                #                           'pc2':self.stand(pcs['pc2'])})
        self.regressors = regressors
        
    def perform_regression(self,dato):
        #Regresion lineal
        y = self.regressors.values
        #y = sm.add_constant(regressors.values)
        lat = dato.lat
        lon = dato.lon
        reg = linear_model.LinearRegression()
        
        campo = dato
        dic_out = {}
        for k in range(len(self.regressors.keys())):
            dic_out[self.regressors.keys()[k]] = {}
            dic_out[self.regressors.keys()[k]]['coef'] = campo.isel(year=0).copy()
            dic_out[self.regressors.keys()[k]]['pval'] = campo.isel(year=0).copy()
            dic_out[self.regressors.keys()[k]]['r2'] = campo.isel(year=0).copy()
            
        for i in range(len(lat)):
            for j in range(len(lon)):
                x = self.create_x(i,j,campo)
                if np.isnan(x).all():
                    res = 0.
                    for k in range(len(self.regressors.keys())):
                        dic_out[self.regressors.keys()[k]]['coef'][i,j] = res
                        dic_out[self.regressors.keys()[k]]['pval'][i,j] = res
                        dic_out[self.regressors.keys()[k]]['r2'][i,j] = res
                else:
                    res = sm.OLS(x,y).fit()
                    for k in range(len(self.regressors.keys())):
                        dic_out[self.regressors.keys()[k]]['coef'][i,j] = res.params[k]
                        dic_out[self.regressors.keys()[k]]['pval'][i,j] = res.pvalues[k]
                        dic_out[self.regressors.keys()[k]]['r2'][i,j] = res.rsquared
                    
        return dic_out
        

