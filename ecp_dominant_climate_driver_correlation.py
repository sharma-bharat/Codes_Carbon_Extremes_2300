#python 2.7
"""
	This code will find the most dominant climate driver based on the correlation coefficients and pvalues using
	the results of the individual correlation coefficints
"""

from __future__ import division                                                    
from scipy import stats
from scipy import ndimage
import glob
import sys
import netCDF4 as nc4
import numpy as np
import datetime as dt
from calendar import monthrange
import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pyplot as plt
#importing my functions
from functions import time_dim_dates, index_and_dates_slicing, geo_idx, mpi_local_and_global_index, create_seq_mat, cumsum_lagged,patch_with_gaps_and_eventsize, norm
from timeit import default_timer as timer
from scipy.stats.stats import pearsonr
from mpi4py import MPI
import pandas as pd
import argparse
import  collections

parser  = argparse.ArgumentParser()
parser.add_argument('--model_config', '-config' 	, help = "Model Configuration"   					, type= str, default= 'w_lulcc') # w_lulcc, wo_lulcc, all
parser.add_argument('--thres_type'  ,   '-th_typ'   , help = "Type of thresholds (independent or misc)" , type= str, default= 'misc'   ) # this is mainly for the misc case study filter
parser.add_argument('--ext_type'    ,   '-ext_typ'  , help = "Type of extreme analysis (pos or neg)"	, type= str, default= 'neg'    )
args = parser.parse_args()
#Running: run ecp_dominant_climate_driver_correlation.py -config wo_lulcc -th_typ misc -ext_typ pos
th_type     = args.thres_type
ext_type    = args.ext_type
conf		= args.model_config

print args
# Misc case (th_type_ for this code is basically based on the 'keys_wrt_wo_lulcc_pos' 
keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
case_name   = ext_type+'_'+conf+'_based_pos_wo_lulcc'

driver_consider = 7
#drivers						 = np.array(['prcp','sm','pme','tmax','tsa','tmin','fclosscol','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider] 
drivers						 = np.array(['prcp','sm','pme','tmax','tsa','tmin','col_fire_closs','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider] 
drivers_names				 = np.array(['Prcp','Soilmoist',r'$P-ET$',r'$T_{max}$',r'$T_{sa}$',r'$T_{min}$','Fire','SPI_gamma_06','SPI_pearson_06','SPEI_gamma_06','SPEI_gamma_12','SPEI_pearson_06'])[:driver_consider]
drivers_code				 = np.array([  10,   20,   30,   40,       50,		60,			70,			80,		90,			100, 110,120])[:driver_consider]
features 			 		 = {}
#features['abr']				 = np.array(['prcp','sm','pme','tmax','tsa', 'tmin','fclosscol','gpp','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider]
features['abr']				 = np.array(['prcp','sm','pme','tmax','tsa', 'tmin','col_fire_closs','gpp','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider]
features['filenames']		 = {}


features['filenames']['gpp'		] = '/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_anomalies_gC.nc'   #just for extraction of (lat,lon,time)

if th_type == 'independent':
	if conf == 'wo_lulcc':
		for dri in drivers:
			features['filenames'][dri] 	= '/home/ud4/CESM1-BGC/'+dri+'/without_land_use_change/corr_triggers/corr_coff_cumlag.nc'

if th_type == 'misc':
	for dri in drivers:
		features['filenames'][dri]	= '/home/ud4/CESM1-BGC/'+dri+'/comparison/pftcon_lulcc/misc/corr_trigger/'+case_name+'/corr_coff_cumlag.nc'

#features['filenames']['prcp'	] = '/home/ud4/CESM1-BGC/prcp/without_land_use_change/corr_triggers/corr_coff_cumlag.nc'

# The name with which the variables are stored in the nc files:
features['Names']			 	= {}
features['Names']['prcp']	 	= 'prcp'
features['Names']['sm']	 	 	= 'sm'
features['Names']['tmax']	 	= 'tmax'
features['Names']['tmin']	 	= 'tmin'
features['Names']['pme']	 	= 'pme'
features['Names']['col_fire_closs']	= 'col_fire_closs'
features['Names']['gpp']	 	= 'gpp'
features['Names']['spi_g06'] 	= 'spi_gamma_06'

#reading time, lat, lon
time 	= nc4.Dataset(features['filenames']['gpp'],mode = 'r')['time']
lat		= nc4.Dataset(features['filenames']['gpp'],mode = 'r')['lat']
lon		= nc4.Dataset(features['filenames']['gpp'],mode = 'r')['lon']

#reading the NC files
features['ncfiles']			 = {}
for dr in features['abr']: 
	print dr
	features['ncfiles'][dr] = nc4.Dataset(features['filenames'][dr],mode = 'r')

#Checks:
p_value_th				= .1 #or 90% confidence
#creating the empty arrays
(nwin,nlag,nlat,nlon)	= features['ncfiles']['prcp']['corr_coeff'].shape

dom_dri_corr_coeff		= np.ma.masked_all((driver_consider,nwin,nlag,nlat,nlon))
dom_dri_ids 			= np.ma.copy(dom_dri_corr_coeff,True)

corr_dri_mask			= np.ma.ones((nwin,nlag,nlat,nlon))
for dr in features['abr']: 
	corr_dri_mask.mask 	= features['ncfiles'][dr]['corr_coeff'][...].mask +corr_dri_mask.mask

# the array below will add all the values of different time-periods and lags to a 2-d array thus giving us the mask of a 2-d array where calculations should be performed.
corr_dri_mask_2d		= np.ma.copy(corr_dri_mask, True)
corr_dri_mask_2d[corr_dri_mask_2d.mask == True] = 0
corr_dri_mask_2d		= np.ma.sum(corr_dri_mask_2d,axis = 0) #4d to 3d
corr_dri_mask_2d		= np.ma.sum(corr_dri_mask_2d,axis = 0) #3d to 2d
corr_dri_mask_2d		= np.ma.masked_equal(corr_dri_mask_2d, 0)

lt_ln_mat   			= create_seq_mat()
corr_dri_mask_1d     	= lt_ln_mat[~corr_dri_mask_2d[...].mask]

wins 	= np.array(range(nwin))
lags 	= np.array(range(nlag))

for win in wins:
	for lg in lags:
		for pixel in corr_dri_mask_1d:
			lt,ln=np.argwhere(lt_ln_mat == pixel)[0]
			dri_id = []
			dri_cc = []
			dri_pv = []
			for idx,dr in enumerate (features['abr']):
				dri_cc.append(features['ncfiles'][dr]['corr_coeff'][win,lg,lt,ln])
				dri_id.append(drivers_code[idx])
				dri_pv.append(features['ncfiles'][dr]['p_value'][win,lg,lt,ln])
			dri_id	= np.array(dri_id)
			dri_cc 	= np.array(dri_cc)
			dri_pv	= np.array(dri_pv)
			if (np.unique(dri_pv < p_value_th, return_counts= 1)[0][0] == False and np.unique(dri_pv < p_value_th, return_counts= 1)[1][0] == driver_consider):
				dri_dom_cc = np.ma.masked_all(driver_consider)
				dri_dom_id = np.ma.masked_all(driver_consider)
			else:
				dri_cp 		= dri_cc[dri_pv < p_value_th]				# driver's corr coeff with confidence
				dri_idp 	= dri_id[dri_pv < p_value_th]				# driver's id codes with confidence
				dri_cp_abs 	= np.abs(dri_cp)							# absolute driver's corr coeff with confidence
				dri_cp_des 	= dri_cp_abs[np.argsort(-dri_cp_abs)]	# absolute driver's corr coeff with confidence in desending order
				dri_dom_cc 	= []
				dri_dom_id	= []
				for c in  dri_cp_des:
					dri_dom_cc.	append(dri_cp[np.argwhere(dri_cp_abs == c)[0][0]])
					dri_dom_id.	append(dri_idp[np.argwhere(dri_cp_abs == c)[0][0]])

			dri_dom_cc_ma 		= np.ma.masked_all(driver_consider)
			dri_dom_id_ma 		= np.ma.masked_all(driver_consider)

			dri_dom_cc_ma[:len(dri_dom_cc)]	= dri_dom_cc
			dri_dom_id_ma[:len(dri_dom_id)]	= dri_dom_id
			dom_dri_corr_coeff	[:,win,lg,lt,ln]= dri_dom_cc_ma
			dom_dri_ids			[:,win,lg,lt,ln]= dri_dom_id_ma

if th_type == 'independent':
	out_path = '/home/ud4/CESM1-BGC/corr_triggers/'
if th_type == 'misc':
	out_path = '/home/ud4/CESM1-BGC/corr_triggers/misc/' + case_name +'/'

with nc4.Dataset(out_path+'dominant_driver_correlation.nc',mode="w") as dset: 
	dset        .createDimension( "rank" ,size = drivers.size) 															##################
	dset        .createDimension( "win" ,size = nwin) 															##################
	dset        .createDimension( "lag" ,size = nlag) 															##################
	dset        .createDimension( "lat" ,size = lat.size) 															##################
	dset        .createDimension( "lon" ,size = lon.size)
	t   =   dset.createVariable(varname = "rank"  		,datatype = float	, dimensions = ("rank") 			,fill_value = np.nan)
	v   =   dset.createVariable(varname = "win"  		,datatype = float	, dimensions = ("win") 				,fill_value = np.nan)
	w   =   dset.createVariable(varname = "lag"  		,datatype = float	, dimensions = ("lag") 				,fill_value = np.nan)
	y   =   dset.createVariable(varname = "lat"  		,datatype = float	, dimensions = ("lat") 				,fill_value = np.nan)
	x   =   dset.createVariable(varname = "lon"  		,datatype = float	, dimensions = ("lon") 				,fill_value = np.nan)
	z   =   dset.createVariable(varname = "dri_id"  	,datatype = float	, dimensions = ("rank","win","lag","lat","lon")	,fill_value = np.nan)
	zz  =   dset.createVariable(varname = "dri_coeff"  	,datatype = float	, dimensions = ("rank","win","lag","lat","lon")	,fill_value = np.nan)
	t.axis  =   "T"
	v.axis  =   "V"
	w.axis  =   "W"
	x.axis  =   "X"
	y.axis  =   "Y"
	t[...]	=	range(len(drivers_code))
	v[...]	=	np.array([format(i,'002') for i in range(nwin)]) 
	w[...]  =   np.array([format(i,'002') for i in range(nlag)])
	x[...]  =   lon[...]
	y[...]  =   lat[...]
	z[...]	= 	dom_dri_ids
	zz[...]	= 	dom_dri_corr_coeff
	t.units			=   "rank of dominant drivers"
	t.missing_value = 	np.nan
	v.units			=   "time periods(wins) 0:1850-74, 1:1875-99...."
	v.missing_value = 	np.nan
	x.units         =   lat.units
	x.missing_value =   np.nan
	x.setncattr         ("long_name",lat.long_name)
	y.units         =   lon.units
	y.missing_value =   np.nan
	y.setncattr         ("long_name",lon.long_name)
	z.missing_value =   np.nan
	z.stardard_name =   "Dominant driver ID "
	z.setncattr         ("method",'Multiple linear regression')
	text = ""
	for i in range(driver_consider): text = text+str(drivers_code[i]) + ":" + str(drivers[i]) +","
	z.units         =   text
	zz.missing_value =   np.nan
	zz.stardard_name =   "Dominant driver Coefficents "
	zz.setncattr         ("method",'Linear Regression')
	zz.units         =   "coefficients"

