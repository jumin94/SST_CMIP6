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

#Remove climatology
def anomalias(dato,reference):
    dato_anom = dato.groupby('time.month') - reference.groupby('time.month').mean('time')
    return dato_anom

#Select DJF months
def seasonal_data(data,season):
    # select DJF
    DA_DJF = data.sel(time = data.time.dt.season==season)

    # calculate mean per year
    DA_DJF = DA_DJF.groupby(DA_DJF.time.dt.year).mean("time")
    return DA_DJF



def standardize(dato):
    return (dato - np.mean(dato))/np.std(dato)

def add_box(box,texto,color='black',text_color='black'):
    x_start, x_end = box[0], box[1]
    y_start, y_end = box[2], box[3]
    margin = 0.0007
    margin_fractions = np.array([margin, 1.0 - margin])
    x_lower, x_upper = x_start + (x_end - x_start) * margin_fractions
    y_lower, y_upper = y_start + (y_end - y_start) * margin_fractions
    box_x_points = x_lower + (x_upper - x_lower) * np.array([0, 1, 1, 0, 0])
    box_y_points = y_lower + (y_upper - y_lower) * np.array([0, 0, 1, 1, 0])
    plt.plot(box_x_points, box_y_points, transform=ccrs.PlateCarree(),linewidth=1, color=color, linestyle='-')
    #plt.text(x_start + (x_end - x_start)*0.4, y_start + (y_end - y_start)*0.4, texto,transform=ccrs.PlateCarree( ),color=text_color)
    
def fig_sst_multiple(dato,dato_r2,t,series,titulo,titulo_serie,levels = [np.arange(-1,1.1,.1)]):
    n = len(dato)
    lon = dato[0].lon; lat = dato[0].lat
    fig = plt.figure(figsize=(n*5, n*3),dpi=300,constrained_layout=True)
    data_crs = ccrs.PlateCarree()
    proj = ccrs.PlateCarree(180)
    for ii in range(n):
        #print('mapa',int(ii*2+1))
        level=levels[ii]
        #print(ii+1,np.max(dato[ii].values),np.max(dato[ii].values)/10)
        ax1 = plt.subplot(n,2,int(ii*2+1),projection=proj)
        im1=ax1.contourf(lon, lat,dato[ii].values,level,transform=data_crs,cmap='RdBu_r',extend='both')
        clevels = np.arange(-1,1,0.2)
        r2 =  ax1.contour(lon, lat,dato_r2[ii].values,level,colors='k',linewidths=0.5,transform=data_crs)
        ax1.clabel(r2, fontsize=10, inline=True)
        ax1.set_title(titulo[ii],fontsize=10)
        ax1.add_feature(cartopy.feature.COASTLINE,alpha=.5)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax1.gridlines(crs=data_crs, linewidth=0.3, linestyle='-')
        ax1.set_extent([-60, 120, -40, 40], ccrs.PlateCarree(central_longitude=180))
        #Add gridlines
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.3, color='gray', alpha=0.01, linestyle='-')
        #Add niño boxes
        # Niño 4
        box = [160,210,-5,5]
        add_box(box,'Niño 4','green','green')
        box = [190,240,-5,5]
        add_box(box,'Niño 3.4','red','red')
        box = [210,270,-5,5]
        add_box(box,'Niño 3')
        box = [270,280,-10,0]
        add_box(box,'Niño 1+2')
        box = [140,170,-5,5]
        add_box(box,'West Pacific')
        
        if ii == 1:
            plt1_ax = plt.gca()
            left, bottom, width, height = plt1_ax.get_position().bounds
            colorbar_axes = fig.add_axes([left + 0.6, bottom,0.02, height])
            cbar = fig.colorbar(im1, colorbar_axes, orientation='vertical')
            cbar.set_label(r'K',fontsize=10) #rotation= radianes
            cbar.ax.tick_params(axis='both',labelsize=5)
            
            
        [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(series[ii])),
                                                   series[ii])
        #plot time serie
        ax1 = plt.subplot(n,2,int(ii*2)+2)
        #print('serie',int(ii*2)+2)
        try:
          ax1.plot(t[ii],series[ii], 'r', linewidth=1.5,label='slope: '+str(round(slope,3)))
        except:
           ax1.plot(series[ii], 'r', linewidth=1.5,label='slope: '+str(round(slope,3)))
        plt.axhline(y=0,xmin=0,linestyle='--',linewidth=0.8,color='k',alpha=0.5)
        ax1.set_ylim((np.min(series[ii])-1, np.max(series[ii])+1))
        ax1.legend()
        x = np.arange(0, len(series[ii]))
        try:
          ax1.plot(t[ii],interc+x*slope,'--k')
        except:
          ax1.plot(interc+x*slope,'--k')
        ax1.set_title(titulo_serie[ii])
        fig.tight_layout()

    return fig


def fig_sst_multiple2(dato,dato_r2,t,series,titulo,titulo_serie,levels = [np.arange(-1,1.1,.1)]):
    n = len(dato)
    lon = dato[0].lon; lat = dato[0].lat
    fig = plt.figure(figsize=(n*6, n*3),dpi=200,constrained_layout=True)
    data_crs = ccrs.PlateCarree()
    proj = ccrs.PlateCarree(180)
    for ii in range(n):
        #print('mapa',int(ii*2+1))
        level=levels[ii]
        #print(ii+1,np.max(dato[ii].values),np.max(dato[ii].values)/10)
        ax1 = plt.subplot(n,2,int(ii*2+1),projection=proj)
        im1=ax1.contourf(lon, lat,dato[ii].values,level,transform=data_crs,cmap='RdBu_r',extend='both')
        clevels = np.arange(-1,1,0.2)
        r2 =  ax1.contour(lon, lat,dato_r2[ii].values,level,colors='k',linewidths=0.5,transform=data_crs)
        ax1.clabel(r2, fontsize=10, inline=True)
        ax1.set_title(titulo[ii],fontsize=10)
        ax1.add_feature(cartopy.feature.COASTLINE,alpha=.5)
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax1.gridlines(crs=data_crs, linewidth=0.3, linestyle='-')
        #ax1.set_extent([-60, 120, -40, 40], ccrs.PlateCarree(central_longitude=180))
        ax1.set_extent([-180, 180, -90, 90], ccrs.PlateCarree(central_longitude=180))
        #Add gridlines
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.3, color='gray', alpha=0.01, linestyle='-')
        #Add niño boxes
        # Niño 4
        box = [160,210,-5,5]
        add_box(box,'Niño 4','green','green')
        box = [190,240,-5,5]
        add_box(box,'Niño 3.4','red','red')
        box = [210,270,-5,5]
        add_box(box,'Niño 3')
        box = [270,280,-10,0]
        add_box(box,'Niño 1+2')
        box = [140,170,-5,5]
        add_box(box,'West Pacific')
        
        if ii == 1:
            plt1_ax = plt.gca()
            left, bottom, width, height = plt1_ax.get_position().bounds
            colorbar_axes = fig.add_axes([left + 0.6, bottom,0.02, height])
            cbar = fig.colorbar(im1, colorbar_axes, orientation='vertical')
            cbar.set_label(r'K',fontsize=10) #rotation= radianes
            cbar.ax.tick_params(axis='both',labelsize=5)
            
            
        [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(series[ii])),
                                                   series[ii])
        #plot time serie
        ax1 = plt.subplot(n,2,int(ii*2)+2)
        #print('serie',int(ii*2)+2)
        try:
          ax1.plot(t[ii],series[ii], 'r', linewidth=1.5,label='slope: '+str(round(slope,3)))
        except:
           ax1.plot(series[ii], 'r', linewidth=1.5,label='slope: '+str(round(slope,3)))
        plt.axhline(y=0,xmin=0,linestyle='--',linewidth=0.8,color='k',alpha=0.5)
        ax1.set_ylim((np.min(series[ii])-1, np.max(series[ii])+1))
        ax1.legend()
        x = np.arange(0, len(series[ii]))
        try:
          ax1.plot(t[ii],interc+x*slope,'--k')
        except:
          ax1.plot(interc+x*slope,'--k')
        ax1.set_title(titulo_serie[ii])
        fig.tight_layout()

    return fig

