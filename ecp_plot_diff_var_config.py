"""
Changes in GPP varmalies for 21st, 22nd and 23rd century
-The benchmark time window is 1975-1999
-The sum of GPP for 2075-2099, 2175-2199 and 2275-2299 will be calculated based on defined percentile (1.0)
-The color bar might have to be changed to not to show the tail events.

"""
from    __future__ import division                                                    
from    scipy import stats
from    scipy import ndimage
import  glob
import  sys
import  netCDF4 as nc4
import  numpy as np
import  pandas as pd
import  datetime as dt
from    calendar import monthrange
import  matplotlib as mpl
mpl.use('Agg')
import  matplotlib.pyplot as plt
#importing my functions
from    functions import time_dim_dates, index_and_dates_slicing, geo_idx, mpi_local_and_global_index,ts_lagged,percentile_colorbar
from    timeit import default_timer as timer
from    scipy.stats.stats import pearsonr
from    mpi4py import MPI
import  statsmodels.formula.api as smf
import  collections

import 	argparse

parser	= argparse.ArgumentParser()
parser	. add_argument('--var_diff'		, '-vdiff'	, help = "plot difference of var"						, type= str, default= 'y'		)		
parser	. add_argument('--neg_extremes'	, '-negext'	, help = "plot for negative extremes"					, type= str, default= 'n'		)
parser	. add_argument('--variable'		, '-var'	, help = "variable whose timeseries is to be plotted?"	, type= str, default= 'n'		) # gpp, nep,...
parser	. add_argument('--model_config'	, '-config'	, help = "Model Configuration"   						, type= str, default= 'w_lulcc'	) # w_lulcc, wo_lulcc, all
parser	. add_argument('--var_in_gC'	, '-var_gC'	, help = "If the variable is gC" 						, type= str, default= 'n'		) # gC, ori
parser  . add_argument('--thres_type'   , '-th_typ' , help = "Type of thresholds (independent or misc)" 	, type= str, default= 'misc'    )
args 	= parser.parse_args()

print args
# running: python ecp_plot_diff_var_config.py -vdiff y	-negext y -var gpp -config wo_lulcc -var_gC gC -th_typ misc
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
if args.neg_extremes == 'y':
	nc_data['wo_lulcc']['ano' ]	= nc4.Dataset(paths['in']['wo_lulcc']	+'cesm1bgc_pftcon_'		+variable+'_anomalies_gC.nc'	)
	nc_data['w_lulcc' ]['ano' ]	= nc4.Dataset(paths['in']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+variable+'_anomalies_gC.nc'	)

win_start_years = np.arange(1850,2300,25)

if args.thres_type	== 'misc':
	nc_data['wo_lulcc']['bin_misc']	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_ext_neg_wo_lulcc_based_pos_wo_lulcc.nc')
	nc_data['w_lulcc'] ['bin_misc']	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_ext_neg_w_lulcc_based_pos_wo_lulcc.nc')

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


for conf in configs:
	for unit in var_units:
		lat		= nc_data[conf][unit].variables['lat']
		lon		= nc_data[conf][unit].variables['lon']
		time	= nc_data[conf][unit].variables['time']
		area	= nc_data[conf][unit].variables['area']
		break
	break

#Constraints
per 		= 1.0 	#percentile  
win_len 	= 12*25 #25 year windows
start_yrs 	= np.array([1975])

var_diff	= args.var_diff
calc_negext = args.neg_extremes


#RdGn colormap
import colorsys as cs
val 		= 0.8
Rd			= cs.rgb_to_hsv(1,0,0)
Rd			= cs.hsv_to_rgb(Rd[0],Rd[1],val)
Gn			= cs.rgb_to_hsv(0,1,0)
Gn			= cs.hsv_to_rgb(Gn[0],Gn[0],val)
RdGn		= {'red'  :	((0.0,	0.0,	Rd[0]),
						 (0.5,	1.0,	1.0	 ),
						 (1.0,	Gn[0],	0.0	 )),
			   'green':	((0.0,	0.0, 	Rd[1]),
					     (0.5,	1.0,	1.0	 ),
						 (1.0,	Gn[1],	0.0  )),
			   'blue' : ((0.0,	0.0,	Rd[2]),
						 (0.5,	1.0,	1.0	 ),
						 (1.0,	Gn[2],	0.0	 ))}
plt.register_cmap(name	= 'RdGn',data = RdGn)

#importing dataset

dates_ar    = time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates = np.array([dates_ar[i*win_len] for i in range(int(time.size/win_len))]) 			#list of start dates of 25 year window
end_dates   = np.array([dates_ar[i*win_len+win_len -1] for i in range(int(time.size/win_len))])	#list of end dates of the 25 year window

#to find the argument or index of the benchmark window 1975-99
for i,date in enumerate(start_dates):
	if date.year in start_yrs:
		start_yr_id = i

per			= 1 #percentile for extreme event definition
yrs			= np.array(range(1850,2301))

# Preparing the dictonary for results storage
# -------------------------------------------
results		= {}
var_sum_win	= []
var_negext_sum_win	= []
for conf in configs:
	if (conf in results.keys()) is False:
		results[conf] = {}
	for unit in var_units:
		if (unit in results[conf].keys()) is False:
			results[conf][unit] = {}
		if args.var_diff == 'y':
			results[conf][unit]['var_diff'] = {}
		if args.neg_extremes == 'y':
			results[conf][unit]['var_negext_diff'] = {}

if (args.neg_extremes == 'y' and args.thres_type != 'misc'):
	var	= nc_data[conf]['ano'].variables[variable]
	for i in range(len(start_dates)):
		print "Calculating negexts for the win %d"%i 
		idx_loc, dates_loc  = index_and_dates_slicing(dates_ar,start_dates[i],end_dates[i])
		var_loc             = var[...][idx_loc[0]:idx_loc[-1]+1,:,:]
		threshold           = np.percentile(var_loc[var_loc.mask == False],per)
		bin_ano_neg         = var_loc < threshold
		ano_loc_ext         = bin_ano_neg * var_loc  #these are basically the anomalies
		var_negext_sum_win	. append (np.ma.sum(ano_loc_ext,axis =0))

if (args.neg_extremes == 'y' and args.thres_type == 'misc'):
	var	= nc_data[conf]['ano'].variables[variable]
	bin_ext = nc_data[conf]['bin_misc'].variables['gpp_bin_ext']
	for i in range(len(start_dates)):
		print "Calculating negexts for the win %d"%i 
		idx_loc, dates_loc  = index_and_dates_slicing(dates_ar,start_dates[i],end_dates[i])
		var_loc             = var[...][idx_loc[0]:idx_loc[-1]+1,:,:]
		bin_ext_loc         = bin_ext[...][idx_loc[0]:idx_loc[-1]+1,:,:]
		ano_loc_ext         = bin_ext_loc * var_loc  #these are basically the anomalies
		var_negext_sum_win	. append (np.ma.sum(ano_loc_ext,axis =0))


if var_diff == 'y':
	var             = nc_data[conf][unit].variables[variable]
	for i in range(len(start_dates)):
		print "Calculating sum ori fluxes for the win %d"%i 
		idx_loc, dates_loc  = index_and_dates_slicing(dates_ar,start_dates[i],end_dates[i])
		var_loc             = var[...][idx_loc[0]:idx_loc[-1]+1,:,:]
		var_sum_win			. append (np.ma.sum(var_loc, axis = 0))
	


if var_diff == 'y':
	var_sum_diff   = [] #calculate the difference of var from the base win 1975-1999
	for i in range(len(var_sum_win)):
		diff        	= np.array(var_sum_win[i]-var_sum_win[start_yr_id])
		var_sum_diff   . append(diff	/ 10**12)
	var_sum_diff		= np.ma.array(var_sum_diff)
	results[conf][unit]['var_diff'] = var_sum_diff 
if calc_negext == 'y':
	var_negext_sum_diff   = [] #calculate the difference of var from the base win 1975-1999
	for i in range(len(var_sum_win)):
		diff        	= np.array(var_negext_sum_win[i]-var_negext_sum_win[start_yr_id])
		var_negext_sum_diff   . append(diff	/ 10**12)
	var_negext_sum_diff		= np.ma.array(var_negext_sum_diff)
	results[conf][unit]['var_negext_diff'] = var_negext_sum_diff 


from    mpl_toolkits.basemap import Basemap
from    matplotlib import cm 
import matplotlib.patches as patches

text = []
for i in start_dates:
	a = i.year
	b = a+24
	c = str(a)+'-'+str(b)[2:]
	text.append(c)

#########################################################################################################################################################

if var_diff == 'y':
	#-------------------*
	#Hardcoding y limits
	if variable == 'gpp':
		ymin = -500
		ymax = 500
	#-------------------*
	for i,data in enumerate (results[conf][unit]['var_diff']):
		print "Plotting ori flux for the win %d"%i 
		fig,ax  = plt.subplots(figsize = (6,2.8),tight_layout=True,dpi=500)
		bmap    = Basemap(  projection  =   'eck4',
    		                lon_0       =   0.,
	    		            resolution  =   'c')
		LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij')
		ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_invalid(data),latlon=True,cmap= 'RdGn',vmax = ymax ,vmin= ymin)
		cbar    = plt.colorbar(ax)
		cbar    .ax.set_ylabel('T'+var.units) #refer to line 97 "diff/10**12"
		#cbar    .ax.set_ylabel(r'kgC/$m^{2}$/yr')
		bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
		bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
		bmap    .drawcoastlines(linewidth = .25,color='lightgrey')
		plt.title(" %s minus %s" %(text[i], text[start_yr_id]))
		if args.thres_type == 'misc':
			fig.savefig('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/diff_plots/th_misc/tot_diff_%s_sce_%s_win_%s.pdf'%(variable,conf, format(i,'02')))    
		else:
			fig.savefig(paths['out'][conf]+'tot_diff_%s_sce_%s_win_%s.pdf'%(variable,conf, format(i,'02')))    
		plt.close(fig)

if calc_negext == 'y':
	#-------------------*
	#Hardcoding y limits
	if variable == 'gpp':
		ymin = -14
		ymax = 14
	#-------------------*
	for i,data in enumerate (results[conf][unit]['var_negext_diff']):

		print "Plotting negexts for the win %d"%i 
		fig_ne,ax  = plt.subplots(figsize = (6,2.8),tight_layout=True,dpi=500)
		bmap    = Basemap(  projection  =   'eck4',
    		                lon_0       =   0.,
	    		            resolution  =   'c')
		LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij')
		ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_invalid(data),latlon=True,cmap= 'RdGn',vmax = ymax ,vmin= ymin)
		cbar    = plt.colorbar(ax)
		cbar    .ax.set_ylabel('T'+var.units) #refer to line 97 "diff/10**12"
		#cbar    .ax.set_ylabel(r'kgC/$m^{2}$/yr')
		bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
		bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
		bmap    .drawcoastlines(linewidth = .25,color='lightgrey')
		plt.title(" %s minus %s" %(text[i], text[start_yr_id]))
		if args.thres_type == 'misc':
			fig_ne.savefig('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/diff_plots/th_misc/tot_diff_negext_%s_sce_%s_win_%s.pdf'%(variable,conf, format(i,'02')))
		else:
			fig_ne.savefig(paths['out'][conf]+'th_misc/tot_diff_negext_%s_sce_%s_win_%s.pdf'%(variable,conf, format(i,'02')))
		plt.close(fig_ne)


