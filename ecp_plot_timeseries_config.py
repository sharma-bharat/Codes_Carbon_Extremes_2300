""" 
Bharat Sharma
python 2.7

Running mean timeseries and the variable values at 1899,1999,2099,2199,2299
to plot the timeseries of the variables of the CESM1-BGC
Running mean(not center) : 5 years meaning the mean of 1,2,3,4,5th values will be saved at 5th location
"""
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
mpl.use('Agg')

import matplotlib.pyplot as plt
#importing my functions
from functions import time_dim_dates, index_and_dates_slicing, geo_idx, mpi_local_and_global_index
from timeit import default_timer as timer
from scipy.stats.stats import pearsonr
from mpi4py import MPI
import pandas as pd

#Inputs:
import argparse
parser  = argparse.ArgumentParser()
parser.add_argument('--variable'	, '-var'		, help = "variable whose timeseries is to be plotted?"	, type= str, default= 'n'		) # gpp, nep,...
parser.add_argument('--model_config', '-config' 	, help = "Model Configuration"   						, type= str, default= 'w_lulcc'	) # w_lulcc, wo_lulcc, all
parser.add_argument('--var_in_gC'	, '-var_gC'		, help = "If the variable is gC" 						, type= str, default= 'n'		) # gC, ori
parser.add_argument('--plt_global'	, '-globe'		, help = "plot_global time series of the variable" 		, type= str, default= 'pgts'	) # n, pgts
parser.add_argument('--lat_range'	, '-lt_range'	, help = "lat range e.g. -30to-60 (Maxrange: -90to90)"	, type= str, default= 'n'		) # n, -90to90
parser.add_argument('--lon_range'	, '-ln_range'	, help = "lon range e.g. -30to-60 (Maxrange: 0to360)"	, type= str, default= 'n'		) # n, 0to360
parser.add_argument('--plt_pft'		, '-pft'		, help = "plot different pft plots?000"					, type= str, default= 'n'		) # n, y
args = parser.parse_args()

# Running : python ecp_plot_timeseries_config.py -var wood_harvestc -var_gC gC -config all 
#temp
print args
variable 		= args.variable #climate variable original

paths						= {}
paths['in' ]				= {}
paths['in' ]['wo_lulcc'] 	= '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
paths['in' ]['w_lulcc' ] 	= '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'
paths['out']				= {}
paths['out']['wo_lulcc'] 	= '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/results/'
paths['out']['w_lulcc' ] 	= '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/results/'
paths['out']['comp'    ] 	= '/home/ud4/CESM1-BGC/'+variable+'/comparison/pftcon_lulcc/'


nc_data						= {}
nc_data['wo_lulcc']			= {}
nc_data['w_lulcc' ]			= {}
nc_data['wo_lulcc']['ori']	= nc4.Dataset(paths['in']['wo_lulcc']	+'cesm1bgc_pftcon_'		+variable+'.nc'		)
nc_data['w_lulcc' ]['ori']	= nc4.Dataset(paths['in']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+variable+'.nc'		)
if args.var_in_gC == "gC":
	nc_data['wo_lulcc']['gC' ]	= nc4.Dataset(paths['in']['wo_lulcc']	+'cesm1bgc_pftcon_'		+variable+'_gC.nc'	)
	nc_data['w_lulcc' ]['gC' ]	= nc4.Dataset(paths['in']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+variable+'_gC.nc'	)


if args.model_config == 'wo_lulcc': 	#without landuse and landcover change
	configs = ['wo_lulcc']
	if args.var_in_gC =='n':
		var_units = ['ori']
	else: 
		var_units = ['gC']
elif args.model_config == 'w_lulcc': 		#with landuse and land cover change
	configs = ['w_lulcc']
	if args.var_in_gC =='n':
		var_units = ['ori']
	else: 
		var_units = ['gC']
elif args.model_config == 'all':
	configs = ['wo_lulcc','w_lulcc']
	if args.var_in_gC =='n':
		var_units = ['ori']
	else: 
		var_units = ['gC']
elif args.model_config == 'all':
	configs = ['wo_lulcc','w_lulcc']
	if args.var_in_gC =='all':
		var_units = ['ori','gC']
else:
	print "Error: Enter the model configuration and the unit types correctly"



#...........save_path 		= in_path+'results/'

#.........nc_flux_ori		= nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'.nc',mode='r')
#.........nc_flux_gC		= nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'_gC.nc',mode='r')

"""	Timeseries plot types:
	1.	Monthly gC						'ts_var_mon_gC'
	2. 	Monthly area weighted mean		'ts_var_mon_awm_gC'
	3.	Annual gC						'ts_var_yr_gC'
	4. 	Annual area weighted mean		'ts_var_yr_awm_gC'
	5.  5-yr rolling mean				'ts_var_5yr_ma_gC'
"""

area	= nc_data['wo_lulcc']['ori'].variables['area']
for conf in configs:
	for unit in var_units:
		lat		= nc_data[conf][unit].variables['lat']
		lon		= nc_data[conf][unit].variables['lon']
		time	= nc_data[conf][unit].variables['time']
		#area	= nc_data[conf][unit].variables['area']
		break
	break

#...lat				= nc_flux_gC.variables['lat']
#...lon				= nc_flux_gC.variables['lon']
#...time			= nc_flux_gC.variables['time']
if (args.lat_range != 'n' and args.lon_range != 'n'):
	lt_range_start	= float((args.lat_range).split('to')[0][1:])
	lt_range_end	= float((args.lat_range).split('to')[1])
	ln_range_start	= float((args.lon_range).split('to')[0])
	ln_range_end	= float((args.lon_range).split('to')[1])
	print type(lt_range_start),lt_range_start
	lat_idx_start	= geo_idx(lt_range_start,lat[...])
	lat_idx_end		= geo_idx(lt_range_end,lat[...])
	lon_idx_start   = geo_idx(ln_range_start,lon[...])
	lon_idx_end   	= geo_idx(ln_range_end,lon[...])

yrs                 = np.array(range(1850,2301))
yrs_selection		= np.array([1899,1999,2099,2199,2299])
yrs_sel_idx			= np.array([np.where(yrs==yr)[0][0] for yr in yrs_selection])	
results				= {}
for conf in configs:
	if (conf in results.keys()) is False:
		results[conf] = {}
	for unit in var_units:
		if (unit in results[conf].keys()) is False:
			results[conf][unit] = {}
		if args.plt_global == 'pgts':
			results[conf][unit]['pgts'] = {}

			var					= nc_data[conf][unit].variables[variable]
			results[conf][unit]['pgts']["var_units"			]			= var.units				#actual units 
			results[conf][unit]['pgts']['mon_global'	]		= np.array([np.sum(var[i,:,:]) for i in range(var.shape[0])]) 								#gram	
			results[conf][unit]['pgts']['mon_awm_global']		= np.array([np.average(var[i,:,:],weights = area[...]) for i in range(var.shape[0])]) 		#gram
			results[conf][unit]['pgts']['yr_global_av'		]	= np.array( [np.mean(results[conf][unit]['pgts']['mon_global'][i*12:(i*12)+12]) for i in range(len(range(1850,2301)))]) #gram
			results[conf][unit]['pgts']['yr_awm_global_av' ]	= np.array([np.average(np.mean(var[i*12:(i*12+12),:,:],axis = 0),weights = area[...]) for i in range(len(range(1850,2301)))])
			results[conf][unit]['pgts']['yr_global_tot'	]		= np.array( [np.sum(results[conf][unit]['pgts']['mon_global'][i*12:(i*12)+12]) for i in range(len(range(1850,2301)))]) #gram
			results[conf][unit]['pgts']['yr_awm_global_tot']	= np.array([np.average(np.sum(var[i*12:(i*12+12),:,:],axis = 0),weights = area[...]) for i in range(len(range(1850,2301)))])
			
			var_pd_global_tot									= pd.Series(results[conf][unit]['pgts']['yr_global_tot']) #g/y
			results[conf][unit]['pgts']['rm_5yr_global_tot'	]	= var_pd_global_tot.rolling(window=5,center = False).mean()# 5 year rolling mean
			var_pd_global_av									= pd.Series(results[conf][unit]['pgts']['yr_awm_global_av']) #g/y
			results[conf][unit]['pgts']['rm_5yr_global_av'	]	= var_pd_global_av.rolling(window=5,center = False).mean()# AWM 5 year rolling mean
			results[conf][unit]['pgts'][variable+'_values']		= np.array([results[conf][unit]['pgts']['rm_5yr_global_av'][idx] for idx in yrs_sel_idx])

# Updating the Units of the variable
# ----------------------------------
if nc_data[conf][unit].variables[variable].units.split('/')[-1] == 's':
	var_unit = nc_data[conf][unit].variables[variable].units [:-2] + '/year'
elif nc_data[conf][unit].variables[variable].units.split('/')[-1] == 'month':
	var_unit = nc_data[conf][unit].variables[variable].units [:-6] + '/year'
else:
	var_unit = nc_data[conf][unit].variables[variable].units + '/year'

# Saving the value of the variable at selected years:
# -------------------------------------------------
for conf in configs:
	val_dict  = {}
	val_dict [var_unit] = results[conf][unit]['pgts'][variable+'_values']
	data = pd.DataFrame(val_dict)
	data
	data.index = yrs_selection
	data.to_csv(paths['out']['comp'] + 'values_%s_%s_%s.csv'%(variable,conf,unit), sep = ',')

if args.plt_pft == 'y':
	pft_names		= [	'Bare Ground', 'NET Temperate','NET Boreal','NDT Boreal','BET Tropical','BET Temperate','BDT Tropical','BDT Temperate','BDT Boreal','BES Temperate','BDS Temperate','BDS Boreal','C 3 arctic grass','C 3 grass','C 4 grass','Crop1','Crop2']
	pft_idx			= [i for i in range(len(pft_names)-1)]


# Adding legend names:
# --------------------
legends = {}
legends['wo_lulcc'] = "without LULCC"
legends['w_lulcc'] = "with LULCC"


"""
Definitions:
------------
1. var_units:
	gC:		 			for fluxes converted to mass
	n:					for original flux units
2. mon_global:			monthly timeseries of the mean of that variable over all lats and lons (5412)
3. mon_awn_global:		monthly timeseries of the area weighted average of the variables (5412)
4. yr_global_av:		yearly timeseries of the yearly mean of global monthly ts (def.2) (5412/12 = 451)
5. yr_awm_global_av:	yearly timeseries of the yearly mean of the area weighted monthly ts (451)
6. yr_global_tot:		yearly timeseries of the variable summed yearly and globally  (451)_
7. yr_awm_global_tot:	yearly timeseries of the variable summed yearly but area weighted averaged (451)
8. rm_5yr_global_tot:	5-yr rolling mean timeseries of yearly ts of variable summed over time and space
9. rm_5yr_global_av:	5-yr rolling mean timeseries of yearly ts of varible summed over time and area weighted averaged

"""
years_ticks = np.arange(1850,2301,50)
# if you want to plot the fluxes which has same or similar units are the variables below :
# ---------------------------------------------------------------------------------------
if variable in ['gpp','npp','nep','nbp','col_fire_closs','wood_harvestc']:
	print "Running the analysis of %s ...."%variable 
	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for conf in configs:
			for unit in var_units:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot']/10**15, label = legends[conf])
				plt.title ("Global 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel("P"+results[conf][unit]['pgts']["var_units"].split('/')[0]+'/year',fontsize =14)
				plt.xlabel("Years",fontsize=14)
				plt.xticks(years_ticks)
				plt.grid  (True, linestyle='--',linewidth = .5)
				plt.legend()
				#setting y axis limits
				if variable == 'gpp': plt.ylim(110,230)
				if variable == 'nep': plt.ylim(1,8)
				if variable == 'npp': plt.ylim(40,65)
				if variable == 'nbp': plt.ylim(-2,3)
				if variable == 'col_fire_closs': plt.ylim(1.75,3.75)

				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))

	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for unit in var_units:
			fig_rm5yr_global = plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
			plt.style.use("classic")

			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot']/10**15, label = legends[conf], lw = 1.5,
							 color='navy', ls= ('--' if conf =='wo_lulcc' else '-'))
			plt.title ("Global 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel(" Global Integrated "+variable.upper() + " (P"+results[conf][unit]['pgts']["var_units"].split('/')[0]+'/year)', 
						fontsize=14)
			plt.xlabel("Time", fontsize =14)
			plt.xticks(years_ticks)
			plt.legend(loc=0)
			plt.grid  (which='both', linestyle='--', linewidth='0.3', color='gray')
			#setting y axis limits
			if variable == 'gpp': plt.ylim(110,230)
			if variable == 'nep': plt.ylim(1,8)
			if variable == 'npp': plt.ylim(40,65)

			if variable == 'col_fire_closs': plt.ylim(1.75,3.75)

			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))


elif variable in ['pme','prcp','sm']:
	print "Running the analysis of %s ...."%variable 
	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for conf in configs:
			for unit in var_units:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot'], label = legends[conf])
				plt.title ("Global 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel(results[conf][unit]['pgts']["var_units"])
				plt.xlabel("Years")
				plt.grid  (True, linestyle='--',linewidth = .5)
				plt.legend()
				#setting y axis limits
				if variable == 'gpp': plt.ylim(110,230)
				
				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))
							
				print " plotting the global average yearly rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_av_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_av'], label = legends[conf])
				plt.title ("Global AWM 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel(results[conf][unit]['pgts']["var_units"])
				plt.xlabel("Years")
				plt.grid  (True, linestyle='--',linewidth = .5)
				plt.legend()
				#setting y axis limits
				if variable == 'gpp': plt.ylim(110,230)

				fig_rm5yr_av_global.savefig(paths['out'][conf] + 'ts_global_awm_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_av_global.savefig(paths['out'][conf] + 'ts_global_awm_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))

	#Timeseries : Rolling Mean Global on the same plot
	if args.plt_global == 'pgts':
		for unit in var_units:
			fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)

			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot'], label = legends[conf])
			plt.title ("Global 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel(results[conf][unit]['pgts']["var_units"])
			plt.xlabel("Years")
			plt.legend()
			plt.grid  (True, linestyle='--',linewidth = .5)
			#setting y axis limits
			if variable == 'gpp': plt.ylim(110,230)

			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))
		
	# Time series : Rolling Mean of Mean yearly global on the same plot
			fig_rm5yr_av_global = plt.figure(tight_layout = True, dpi = 400)

			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_av'], label = legends[conf])
			plt.title ("Global AWM 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel(results[conf][unit]['pgts']["var_units"])
			plt.xlabel("Years")
			plt.legend()
			plt.grid  (True, linestyle='--',linewidth = .5)
			#setting y axis limits
			if variable == 'gpp': plt.ylim(110,230)

			fig_rm5yr_av_global.savefig(paths['out']['comp'] + 'ts_global_awm_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_av_global.savefig(paths['out']['comp'] + 'ts_global_awm_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))
	
elif variable in ['tmax','tmin','tsa']:
	print "Running the analysis of %s ...."%variable 
	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for conf in configs:
			for unit in var_units:
				print " plotting the global AWM rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_awm_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_av']-273.16, label = legends[conf])
				plt.title ("Global AWM 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel(r'$^\circ$ Celsius')
				plt.xlabel("Years")
				plt.grid  (True, linestyle='--',linewidth = .5)
				plt.legend()
				#setting y axis limits
				if variable == 'gpp': plt.ylim(110,230)

				fig_rm5yr_awm_global.savefig(paths['out'][conf] + 'ts_global_AWM_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_awm_global.savefig(paths['out'][conf] + 'ts_global_AWM_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))

	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for unit in var_units:
			# Time series : Rolling Mean of Mean yearly global on the same plot
			fig_rm5yr_awm_global = plt.figure(tight_layout = True, dpi = 400)

			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_av']-273.16, label = legends[conf])
			plt.title ("Global AWM 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel(r'$^\circ$ Celsius')
			plt.xlabel("Years")
			plt.legend()
			plt.grid  (True, linestyle='--',linewidth = .5)
			#setting y axis limits
			if variable == 'gpp': plt.ylim(110,230)

			fig_rm5yr_awm_global.savefig(paths['out']['comp'] + 'ts_global_awm_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_awm_global.savefig(paths['out']['comp'] + 'ts_global_awm_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))

# If the variable is a ratio and has no units
# ==========================================
elif variable in ['tlai']:
	print "Running the analysis of %s ...."%variable 
	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for conf in configs:
			for unit in var_units:
				print " plotting the global AWM rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_awm_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_av'], label = legends[conf])
				plt.title ("Global AWM 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel(variable.upper())
				plt.xlabel("Years")
				plt.grid  (True, linestyle='--',linewidth = .5)
				plt.legend()
				#setting y axis limits
				if variable == 'gpp': plt.ylim(110,230)

				fig_rm5yr_awm_global.savefig(paths['out'][conf] + 'ts_global_AWM_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_awm_global.savefig(paths['out'][conf] + 'ts_global_AWM_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))

	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for unit in var_units:
			# Time series : Rolling Mean of Mean yearly global on the same plot
			fig_rm5yr_awm_global = plt.figure(tight_layout = True, dpi = 400)
			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_av'], label = legends[conf])
			plt.title ("Global AWM 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel(variable.upper())
			plt.xlabel("Years")
			plt.legend()
			plt.grid  (True, linestyle='--',linewidth = .5)
			#setting y axis limits
			if variable == 'gpp': plt.ylim(110,230)

			fig_rm5yr_awm_global.savefig(paths['out']['comp'] + 'ts_global_awm_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_awm_global.savefig(paths['out']['comp'] + 'ts_global_awm_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))

# If the variable is a mass i.e. the units are gC (in carbon pools)
# ================================================================
elif variable in ['totecosysc']:
	print "Running the analysis of %s ...."%variable 
	print nc_data[conf][unit].variables[variable].units
	print nc_data[conf][unit]
	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for conf in configs:
			for unit in var_units:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot']/10**15, label = legends[conf])
				plt.title ("Global 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel("P"+nc_data[conf][unit].variables[variable].units)
				plt.xlabel("Years")
				plt.legend()
				plt.grid  (True, linestyle='--',linewidth = .5)
				
				#setting y axis limits
				if variable == 'gpp': plt.ylim(110,230)
				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))

	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for unit in var_units:
			fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)

			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot']/10**15, label = legends[conf])
			plt.title ("Global 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel("P"+ nc_data[conf][unit].variables[variable].units)
			plt.xlabel("Years")
			plt.legend()
			plt.grid  (True, linestyle='--',linewidth = .5)
			#setting y axis limits
			if variable == 'gpp': plt.ylim(110,230)

			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))

			
else:
	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for conf in configs:
			for unit in var_units:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot']/10**15, label = legends[conf])
				plt.title ("Global 5 yr Running Mean %s for %s (1850-2300)"%(variable.upper(),conf))
				plt.ylabel("P"+results[conf][unit]['pgts']["var_units"].split('/')[0]+'/year')
				plt.xlabel("Years")
				plt.grid  (True, linestyle='--',linewidth = .5)
				plt.legend()
				#setting y axis limits
				if variable == 'npp': plt.ylim(40,65)

				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.pdf'%(variable,conf,unit))
				fig_rm5yr_global.savefig(paths['out'][conf] + 'ts_global_%s_%s_5yr_rm_%s.png'%(variable,conf,unit))

	#Timeseries : Rolling Mean Global
	if args.plt_global == 'pgts':
		for unit in var_units:
			fig_rm5yr_global = plt.figure(tight_layout = True, dpi = 400)

			for conf in configs:
				print " plotting the global rolling mean for config %s and for %s flux"%(conf, unit)
				plt.plot  (yrs, results[conf][unit]['pgts']['rm_5yr_global_tot']/10**15, label = legends[conf])
			plt.title ("Global 5 yr Running Mean %s for with vs without LULCC (1850-2300)"%(variable.upper()))
			plt.ylabel("P"+results[conf][unit]['pgts']["var_units"].split('/')[0]+'/year')
			plt.xlabel("Years")
			plt.legend()
			plt.grid  (True, linestyle='--',linewidth = .5)
			#setting y axis limits
			if variable == 'npp': plt.ylim(40,65)

			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.pdf'%(variable,unit))
			fig_rm5yr_global.savefig(paths['out']['comp'] + 'ts_global_%s_pftcon_lulcc_5yr_rm_%s.png'%(variable,unit))

		
# plotting the rolling mean of CO2 concentration
# ----------------------------------------------
path_co2	= "/home/ud4/CESM1-BGC/rcp85_co2/"
filename_co2= "%srcp85_co2.csv"%path_co2
#read csv
df_co2 = pd.read_csv(filename_co2, header=None)
s_idx 	= 85
l_idx	= 85+450
time_co2 	= df_co2.iloc[:,0] [s_idx:l_idx+1]
ppm_co2 	= df_co2.iloc[:,1] [s_idx:l_idx+1]


#Timeseries : yearly Global
print ("Plotting the atmospheric CO2 concentration")
fig_co2 = plt.figure(tight_layout = True, dpi = 400, figsize=(9,5))
plt.plot  (time_co2, ppm_co2, label = "CO2 (ppm)")
plt.title (r'Atmospheric CO$_2$ Concentration')
plt.ylabel(r'CO$_2$ (ppm)',fontsize =14)
plt.xlabel("Years",fontsize=14)
plt.xticks(years_ticks)
plt.grid  (True, linestyle='--',linewidth = .5)
#plt.legend()
#setting y axis limits
plt.ylim(200,2000)
fig_co2.savefig(path_co2 + 'ts_CO2_ppm.pdf')
fig_co2.savefig(path_co2 + 'ts_CO2_ppm.png')