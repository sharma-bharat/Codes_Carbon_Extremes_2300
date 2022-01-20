#Bharat Sharma
#python 2.7.13

from __future__ import division                                         
from scipy import stats
from scipy import ndimage
import glob
import sys
import netCDF4 as nc4
import numpy as np
#import pandas as pd
import datetime as dt
from calendar import monthrange
import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
#importing my functions
from functions import time_dim_dates, index_and_dates_slicing, geo_idx
from timeit import default_timer as timer
import argparse


#Just enter the flux that you want to covert to gC
parser  = argparse.ArgumentParser()
parser.add_argument('--variable',       '-var'  ,    help = "Variable to decompose" , type= str, default= 'n')
parser.add_argument('--model_config',   '-mod_conf', help = "Model Configuration"   , type= str, default= 'w_lulcc')
args = parser.parse_args()

variable= args.variable
# import datasets
#new : 
if args.model_config == 'w_lulcc': 					#with landuse and land cover change
	in_path     = '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'
	nc_var      = nc4.Dataset(in_path+'cesm1bgc_with_lulcc_'+variable+'.nc',mode='r')
elif args.model_config == 'wo_lulcc':   			#without landuse and landcover change
	in_path = '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
	nc_var      = nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'.nc',mode='r')

nc_gpp      = nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/'+'cesm1bgc_pftcon_gpp.nc',mode='r')
####################################################################################
ga_names    = nc_var.ncattrs()
#Reading data from nc_var/anomalies nc file
var     	= nc_var.variables[variable]
#var     	= nc_var.variables[variable.upper()]
lat     	= nc_var.variables['lat']
lon    	 	= nc_var.variables['lon']
time        = nc_var.variables['time']
time_bounds = nc_gpp.variables['time_bounds']
#Reading data from nc_gpp/original nc file
area        = nc_gpp.variables['area']
lf      	= nc_gpp.variables['landfrac']
area_act    = area[...] * lf[...]           #area_act is the effective or actual area of that pixels
time_days   = [int(time_bounds[i][1]-time_bounds[i][0]) for i in range(time_bounds.shape[0])]   #time_bounds from the variables of timebound
time_sec    = np.array(time_days)*24*3600   
vol_m2s     = time_sec[:,np.newaxis,np.newaxis] * area_act *(10**6)   #units vol: m^2*s
var_gC      = var * vol_m2s                                 #units var: gC, masked array with fill_value = 1e+36

if args.model_config == 'wo_lulcc':
	out_filename = in_path+'cesm1bgc_pftcon_'+variable+'_gC.nc'
elif args.model_config == 'w_lulcc':
	out_filename = in_path+'cesm1bgc_with_lulcc_'+variable+'_gC.nc'

with nc4.Dataset(out_filename,mode = 'w') as dset:
	for ga_name in  ga_names:
		gval=   nc_var.__dict__[ga_name]
		dset    .setncattr(ga_name,gval)
	dset        .createDimension( "time",size = time.size)
	dset        .createDimension( "lat" ,size = lat.size)
	dset        .createDimension( "lon" ,size = lon.size)
	t   =   dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time"), fill_value = 1e+36)
	x   =   dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon") , fill_value = 1e+36)
	y   =   dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat") , fill_value = 1e+36)
	z   =   dset.createVariable(varname = variable  ,datatype = float, dimensions = ("time","lat","lon"),fill_value = 1e+36)
	t.axis  =   "T"
	x.axis  =   "X"
	y.axis  =   "Y"
	t[...]  =   time  [...]
	x[...]	=	lon   [...]
	y[...]	=	lat   [...]
	z[...]	=	var_gC[...]
	z.missing_value	= 1e+36
	z.stardard_name = variable+" in gC" + var.long_name
	z.setncattr		  ("cell_methods",var.cell_methods)
	z.units			= "gC/month"
	x.units         =   lat.units
	x.missing_value =   1e+36
	x.setncattr         ("long_name",lat.long_name)
	y.units         =   lon.units
	y.missing_value =   1e+36
	y.setncattr         ("long_name",lon.long_name)
	t.units         =   time.units
	t.setncattr         ("calendar",time.calendar)
	t.setncattr         ("long_name",time.long_name)
	t.missing_value =	1e+36
	
