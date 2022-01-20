# Bharat Sharma
# python 2.7
# Different correlation approaches for all the pixels with atleast one neg and one pos TCE in time period 2000-24

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
from functions import time_dim_dates, index_and_dates_slicing, geo_idx, mpi_local_and_global_index, create_seq_mat, cumsum_lagged,patch_with_gaps_and_eventsize, norm, cum_av_lagged
from timeit import default_timer as timer
from scipy.stats.stats import pearsonr
from mpi4py import MPI
import pandas as pd
import argparse
import  collections

import timeit

start_time = timeit.default_timer()

#Inputs:
#------
parser  = argparse.ArgumentParser()
parser.add_argument('--drivers'		, '-dri'  	, help = "Type all the variables you want to plot( separated by '-' minus sign" , type= str, default= "prcp-sm-pme-tmax-tsa-tmin-col_fire_closs")
parser.add_argument('--variable'   	, '-var'    , help = "Anomalies of carbon cycle variable" 	, type= str		, default= 'gpp'	)
parser.add_argument('--var_in_gC'	, '-var_gC'	, help = "Unit:If the variable is gC" 			, type= str		, default= 'gC'		) # gC, ori
parser.add_argument('--model_config', '-config' , help = "Model Configuration"   				, type= str		, default= 'w_lulcc') # w_lulcc, wo_lulcc, all
parser.add_argument('--thres_type'	, '-th_typ' , help = "Type of thresholds (independent or misc)", type= str	, default= 'misc'	) # this is mainly for the misc case study filter
parser.add_argument('--ext_type'	, '-ext_typ', help = "Type of extreme analysis (pos or neg)", type= str		, default= 'neg'	) #pos/neg
#parser.add_argument('--lat_idx'		, '-lat_idx', help = ' Lat coordinate of the location'		, type = int	, default = 102	    ) 
#parser.add_argument('--lon_idx'		, '-lon_idx', help = ' Lon coordinate of the location'		, type = int	, default = 26	    ) 
parser.add_argument('--plot_win'    , '-pwin'	, help = "which time period to plot? 2000-24:win 06", type= int	, default=  6   	)
parser.add_argument('--tce_len_p'   , '-tce_p'	, help = "percent len of TCE to look at (trigger)", type= float	, default=  .25   	)

args = parser.parse_args()
print args
#running run ecp_correlation_methods_pixel_win.py -dri prcp-sm-pme-tmax-tsa-tmin-col_fire_closs -th_typ misc -pwin 6 -tce_p 0.25 -config wo_lulcc -ext_typ neg 

drivers_list=  (args.drivers).split("-")
variable	= args.variable 	# variable original
th_type		= args.thres_type
ext_type	= args.ext_type
pwin    	= args.plot_win
#lat_in		= args.lat_idx
#lon_in		= args.lon_idx
conf		= args.model_config
tce_len_per = args.tce_len_p

#	import datasets
paths							= {}
paths['in' ]					= {}
paths['in' ]['var']				= {}
paths['in' ]['var']['wo_lulcc'] = '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
paths['in' ]['var']['w_lulcc' ] = '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'

paths['out']					= {}
#paths['out']['wo_lulcc'] 		= '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/results/'
#paths['out']['w_lulcc' ] 		= '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/results/'

nc_data								= {}
nc_data['var']						= {}
nc_data['var']['wo_lulcc']			= {}
nc_data['var']['w_lulcc' ]			= {}
nc_data['var']['wo_lulcc']['ori']	= nc4.Dataset(paths['in']['var']['wo_lulcc']	+'cesm1bgc_pftcon_'		+variable+'.nc'		)
nc_data['var']['w_lulcc' ]['ori']	= nc4.Dataset(paths['in']['var']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+variable+'.nc'		)
if args.var_in_gC == "gC": #Anomalies of gpp in gC
	nc_data['var']['wo_lulcc']['gC' ]	= nc4.Dataset(paths['in']['var']['wo_lulcc']	+'cesm1bgc_pftcon_'		+variable+'_anomalies_gC.nc'	)
	nc_data['var']['w_lulcc' ]['gC' ]	= nc4.Dataset(paths['in']['var']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+variable+'_anomalies_gC.nc'	)

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
		lat			= nc_data['var'][conf][unit].variables['lat' ]
		lon			= nc_data['var'][conf][unit].variables['lon' ]
		time		= nc_data['var'][conf][unit].variables['time']
		ano			= nc_data['var'][conf][unit].variables[variable]
		#time_bounds	= nc_data['var'][conf][unit].variables['time_bounds']
		break
	break

##############################################
window 		= 25 #years
win_len     = 12 * window            #number of months in window years
nwin		= int(time.size/win_len) #number of windows

wins    	= np.array([format(i,'02' ) for i in range(nwin)])
#number of months you want to use for lags : 0(during) to 12months)
#max_lag		=	12 #months
#lags    	= np.array([format(i,'02' ) for i in range(max_lag +1)])
dates_ar    = time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates	= [dates_ar[i*win_len] for i in range(nwin)]#list of start dates of 25 year window
end_dates   = [dates_ar[i*win_len+win_len -1] for i in range(nwin)]#list of end dates of the 25 year window

# Calculation of anomalies
idx_dates_win= []   #indicies of time in 30yr windows
dates_win   = []    #sel dates from time variables in win_len windows

# Reading driver filenames!
#--------------------------
paths['in' ]['driver_ano']				= {}
paths['in' ]['driver_ano']['wo_lulcc'] = {}
paths['in' ]['driver_ano']['w_lulcc' ] = {}
for dri in drivers_list:
	paths['in' ]['driver_ano']['wo_lulcc'][dri] = '/home/ud4/CESM1-BGC/'+dri+'/without_land_use_change/'
	paths['in' ]['driver_ano']['w_lulcc' ][dri] = '/home/ud4/CESM1-BGC/'+dri+'/with_land_use_change/'

# Reading the nc data:
# --------------------
nc_data['driver_ano']					= {}
nc_data['driver_ano']['wo_lulcc']		= {}
nc_data['driver_ano']['w_lulcc' ]		= {}
for dri in drivers_list:
	nc_data['driver_ano']['wo_lulcc'][dri]		= nc4.Dataset(paths['in']['driver_ano']['wo_lulcc'][dri] + 'cesm1bgc_pftcon_'+ dri + '_anomalies.nc'	 )
	nc_data['driver_ano']['w_lulcc' ][dri]		= nc4.Dataset(paths['in']['driver_ano']['w_lulcc' ][dri] + 'cesm1bgc_with_lulcc_'+ dri + '_anomalies.nc' )

# Getting the land frac data and ocean will be masked
lf		= nc_data['var'][conf]['gC']['landfrac'][...]
lf_mask	= np.ma.masked_equal(lf,0)

# Different Cases Under Considerations:
keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
case_name   = ext_type+'_'+conf+'_based_pos_wo_lulcc'

# Saving the dataframe
# --------------------
paths['out']['TEST'] = '/home/ud4/CESM1-BGC/attr_2019/misc/'+conf+'/'


# Correlation Testing
# -------------------
"""
	The correlations are tested for following cases:
	1. fts 		: full time series, the correlations are computed for the complete timeseries
	2. fts_rm	: full time series, the correlations are computed for the complete timeseries with 3 months of running mean
	3. tce_len 	: computing the correlations only at locations of TCE keeping the len constant even at lags
	4. tce_len_b: computing the correlations only at locations of TCE keeping the len constant even at lags (TCE locations for positive as well as negative should be considered)
	5. tce_cum_av_lag_b : when the mean cumulative lag of drivers are correlated with TCE for both positive and negative TCEs
	6. tce_cum_av_lag_trig_b: when the mean cumulative lag of drivers are correlated with the triggers of TCE for both positive and negative TCE ( constrained by what percent of length of TCE to be evaluated)
	7. tce_trig_b : computing the correlations only at locations of triggers of TCE keeping the len constant even at lags (TCE locations for positive as well as negative should be considered)

"""
#lags = [0,1,2,3,4,5]
lags = [0,1,2,3,4]

# For non-normalized correlation
# ------------------------------
corr_cases	= ['fts','fts_rm', 'tce_len','tce_len_b' ,'tce_cum_av_lag_b','tce_cum_av_lag_trig_b' ,'tce_trig_b']
corr		= {}
for cors in corr_cases:
	corr[cors] = {}
	for dri in drivers_list:
		corr[cors][dri] = {}
		for lag in lags:
			corr[cors][dri][lag] = {}

# For normalized correlations
# ---------------------------
corr_norm		= {}
for cors in corr_cases:
	corr_norm[cors] = {}
	for dri in drivers_list:
		corr_norm[cors][dri] = {}
		for lag in lags:
			corr_norm[cors][dri][lag] = {}


# Creating two dataframes (2-D) one to capture the len of the TCE and other for correlation coefficients and p value for different lags:
# -------------------------------------------------------------------------------------------------------------------------------------

df_tce = pd.DataFrame(columns = ['tce_neg', 'tce_pos', 'tce_len_neg', 'tce_len_pos', 'tce_len_tot'])

col_names = []

df_cc_pv = {}
for cors in corr_cases:
	df_cc_pv[cors]= {}

col_names = []
for dri in drivers_list:
	col_names .append(dri+'-cc')
	col_names .append(dri+'-pv')
for cors in corr_cases:
	df_cc_pv[cors] = pd.DataFrame(columns = col_names)
# to find the particular file name
# --------------------------------
win_start_years = np.arange(1850,2300,25)

from    scipy   import  ndimage
for lat_in in range(len(lat[...])):
	for lon_in in range(len(lon[...])):
		#binary of 1s for extreme events for 2000-24
		lag_idx_bin = -1 # which is 4 months, this is done so that the bin TS is same for correlations from lag 0 to 4
		bin_tce_1s	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_1s_'+case_name+'_%d.nc'%win_start_years[pwin])['gpp_TCE_1s'] [lag_idx_bin,:,lat_in,lon_in]
		bin_tce_01s	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_01s_'+case_name+'_%d.nc'%win_start_years[pwin])['gpp_TCE_01s'] [lag_idx_bin,:,lat_in,lon_in]
		if ext_type == 'neg':
			case_name_alt 	= 'pos'+'_'+conf+'_based_pos_wo_lulcc'
			bin_tce_1s_alt	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_1s_'+case_name_alt+'_%d.nc'%win_start_years[pwin])['gpp_TCE_1s'] [lag_idx_bin,:,lat_in,lon_in]
			bin_tce_01s_alt	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_01s_'+case_name_alt+'_%d.nc'%win_start_years[pwin])['gpp_TCE_01s'] [lag_idx_bin,:,lat_in,lon_in]
		if ext_type == 'pos':
			case_name_alt 	= 'neg'+'_'+conf+'_based_pos_wo_lulcc'
			bin_tce_1s_alt	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_1s_'+case_name_alt+'_%d.nc'%win_start_years[pwin])['gpp_TCE_1s'] [lag_idx_bin,:,lat_in,lon_in]
			bin_tce_01s_alt	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_01s_'+case_name_alt+'_%d.nc'%win_start_years[pwin])['gpp_TCE_01s'] [lag_idx_bin,:,lat_in,lon_in]
		
		#TCE
		larray,narray   = ndimage.label(bin_tce_1s,structure = np.ones(3))
		locations       = ndimage.find_objects(larray)
		#TCE _alt
		larray_alt,narray_alt   = ndimage.label(bin_tce_1s_alt,structure = np.ones(3))
		locations_alt       = ndimage.find_objects(larray_alt)
		

		if (narray > 1 and narray_alt > 1): #CHECK  great than 1 because this package return 1 event if all are masked or 0's and with only 1 tce and other 0's
			print "Check Passed"	
			df_tce_index = "%s_%s"%(format(lat_in,'003'),format(lon_in,'003'))
			df_tce = df_tce.reindex(df_tce.index.values.tolist()+[df_tce_index])
			
			# Normalized data
			data_norm = {}
			for dri in drivers_list: # This loop is just for the data frame of normalized data
				data_norm [dri] 	= norm (nc_data['driver_ano'][conf][dri][dri][pwin*win_len:(pwin+1)*win_len, lat_in, lon_in])
				data_norm ['TCE']	= larray
				data_norm ['TCE_alt']	= larray_alt
				data_norm [variable]= norm(nc_data['var'][conf]['gC'][variable][pwin*win_len:(pwin+1)*win_len, lat_in, lon_in])

			df_norm	= pd.DataFrame(data_norm, index = dates_ar[pwin*win_len:(pwin+1)*win_len])
			
			# Case 01: fts
			# ------------	
			cors = 'fts'
			
			df_index = []
			for lag in lags:
				df_cc_pv_index = "%s_%s_%s"%(format(lat_in,'003'),format(lon_in,'003'),format(lag,'02'))
				df_cc_pv[cors] = df_cc_pv[cors].reindex(df_cc_pv[cors].index.values.tolist()+[df_cc_pv_index])
				df_index.append(df_cc_pv_index)
			print df_cc_pv[cors].shape 
			for dri in drivers_list: # This loop is for cc pv calculations
				# Case 01: fts
				# ------------	
				cors = 'fts'
				for lag in lags:
					if lag == 0:
						cc,pv  = stats.linregress(df_norm[dri],df_norm[variable])[2:4]
					if lag > 0:
						cc,pv = stats.linregress(df_norm[dri][:-lag],df_norm[variable][lag:])[2:4]
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)
					print "dri: %s lag: %s"%(dri,lag),cc,pv

				# Case 02: fts_rm (normalized)
				# ----------------------------
				cors = 'fts_rm'
				tmp_pd_dri  = pd.Series(df_norm[dri])
				tmp_pd_var  = pd.Series(df_norm[variable])
				# Rolling mean of 3 months
				tmp_rm_var  = tmp_pd_var.rolling(window=3, center= True).mean()[1:-1]
				tmp_rm_dri  = tmp_pd_dri.rolling(window=3, center= True).mean()[1:-1]
				for lag in lags:
					tmp_pd_dri  = pd.Series(df_norm[dri])
					tmp_pd_var  = pd.Series(df_norm[variable])
					#Rolling mean of 3 months
					tmp_rm_var  = tmp_pd_var.rolling(window=3, center= True).mean()[1:-1]
					if lag == 0:
						cc,pv = stats.linregress(tmp_rm_dri,tmp_rm_var)[2:4]
					else:
						cc,pv = stats.linregress(tmp_rm_dri[:-lag] , tmp_rm_var[lag:])[2:4]
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)


				# Case 03: tce_len (normalized)
				# ---------------------------------
				cors = 'tce_len'
				ts_tce	= df_norm['TCE']
				idxs    = np.asarray(np.arange(len(ts_tce)), dtype = int)
				tce_idxs= np.array([])
				for loc in locations: 
					tce_idxs = np.append(tce_idxs, idxs[loc])

				tce_idxs = np.asarray(tce_idxs, dtype = int)

				for lag in lags:
					if lag == 0:
						ts_dri_tce  = np.array([np.array(df_norm[dri])[i] for i in tce_idxs])
						ts_var_tce  = np.array([np.array(df_norm[variable])[i] for i in tce_idxs])
						cc,pv 		= stats.linregress(ts_dri_tce, ts_var_tce)[2:4]
					else:
						if len(tce_idxs) <= len(lags):
							cc,pv = (np.nan,np.nan)
						elif lag<=tce_idxs[0]:
							ts_dri_tce  = np.array([np.array(df_norm[dri])[i] for i in (tce_idxs-lag)])
							ts_var_tce	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs])
							cc,pv = stats.linregress(ts_dri_tce, ts_var_tce)[2:4]
						elif lag <= tce_idxs[len(lags)]:
							ts_dri_tce  = np.array([np.array(df_norm[dri])[i] for i in (tce_idxs-lag)])
							ts_var_tce	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs])
							cc,pv = stats.linregress(ts_dri_tce, ts_var_tce)[2:4]
						else:
							cc,pv = (np.nan,np.nan)
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)


				# Case 04: tce_len_b (normalized)
				# ---------------------------------
				cors = 'tce_len_b'
				ts_tce	= df_norm['TCE']
				idxs    = np.asarray(np.arange(len(ts_tce)), dtype = int)
				tce_idxs= np.array([])
				for loc in locations: 
					tce_idxs = np.append(tce_idxs, idxs[loc])
				tce_idxs = np.asarray(tce_idxs, dtype = int)
				#alt:
				ts_tce_alt	= df_norm['TCE_alt']
				idxs_alt    = np.asarray(np.arange(len(ts_tce_alt)), dtype = int)
				tce_idxs_alt= np.array([])
				for loc in locations_alt: 
					tce_idxs_alt = np.append(tce_idxs_alt, idxs[loc])
				tce_idxs_alt = np.asarray(tce_idxs_alt, dtype = int)

				tce_idxs_concate	= np.concatenate((tce_idxs, tce_idxs_alt))
				tce_idxs_concate.sort()

				for lag in lags:
					if lag == 0:
						ts_dri_tce_b  	= np.array([np.array(df_norm[dri])[i] for i in tce_idxs_concate])
						ts_var_tce_b  	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
						cc,pv			= stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
					else:
						if len(tce_idxs_concate) <= len(lags):
							cc,pv = (np.nan,np.nan)
						elif lag<=tce_idxs_concate[0]:
							ts_dri_tce_b  	= np.array([np.array(df_norm[dri])[i] for i in (tce_idxs_concate-lag)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv 			= stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						elif lag <= tce_idxs_concate[len(lags)]: 
							ts_dri_tce_b  	= np.array([np.array(df_norm[dri])[i] for i in (tce_idxs_concate-lag)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv 			= stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						else:
							cc,pv = (np.nan,np.nan)
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)


				# Case 05: tce_cum_av_lag_b (normalized)
				# ---------------------------------
				cors = 'tce_cum_av_lag_b'
				ts_tce	= df_norm['TCE']
				idxs    = np.asarray(np.arange(len(ts_tce)), dtype = int)
				tce_idxs= np.array([])
				for loc in locations: 
					tce_idxs = np.append(tce_idxs, idxs[loc])
				tce_idxs = np.asarray(tce_idxs, dtype = int)
				#alt:
				ts_tce_alt	= df_norm['TCE_alt']
				idxs_alt    = np.asarray(np.arange(len(ts_tce_alt)), dtype = int)
				tce_idxs_alt= np.array([])
				for loc in locations_alt: 
					tce_idxs_alt = np.append(tce_idxs_alt, idxs[loc])
				tce_idxs_alt = np.asarray(tce_idxs_alt, dtype = int)
				tce_idxs_concate	= np.concatenate((tce_idxs, tce_idxs_alt))
				tce_idxs_concate.sort()
				if ext_type == 'neg':
					df_tce .loc[df_tce_index,'tce_neg'] 	=  len(locations)
					df_tce .loc[df_tce_index,'tce_pos'] 	=  len(locations_alt)
					df_tce .loc[df_tce_index,'tce_len_neg'] =  len(tce_idxs)
					df_tce .loc[df_tce_index,'tce_len_pos'] =  len(tce_idxs_alt)
					df_tce .loc[df_tce_index,'tce_len_tot'] =  len(tce_idxs_concate)

				for lag in lags:
					if lag == 0:
						cc,pv = (np.nan,np.nan)
					else:
						if len(tce_idxs_concate) <= len(lags):
							cc,pv = (np.nan,np.nan)
						elif lag<=tce_idxs_concate[0]:
							ts_dri_tce_b  = np.array([cum_av_lagged(df_norm[dri],ignore_t0 = True, lag = lag)[i] for i in (tce_idxs_concate)]) # "tce_idxs_concate-lag would be incorrect because lagged effects are considered by cum_av_lagged
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv = stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						elif lag <= tce_idxs_concate[len(lags)]:
							ts_dri_tce_b  = np.array([cum_av_lagged(df_norm[dri],ignore_t0 = True, lag = lag)[i] for i in (tce_idxs_concate)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv = stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						else:
							cc,pv = (np.nan,np.nan)
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)

				# Case 06: tce_cum_av_lag_trig_b (normalized)
				# ---------------------------------
				cors = 'tce_cum_av_lag_trig_b'
				ts_tce	= df_norm['TCE']
				idxs    = np.asarray(np.arange(len(ts_tce)), dtype = int)
				per_trig_len = tce_len_per
				per_tce_len  = []
				for loc in locations:
					per_tce_len . append( int(np.ceil(per_trig_len * bin_tce_1s[loc].size)))
				tce_idxs= np.asarray([],dtype = int)
				for i,loc in enumerate(locations): 
					tce_idxs = np.append(tce_idxs, idxs[loc][:per_tce_len[i]])
				#alt:
				ts_tce_alt	= df_norm['TCE_alt']
				per_tce_len_alt  = []
				for loc in locations_alt:
					per_tce_len_alt . append( int(np.ceil(per_trig_len * bin_tce_1s_alt[loc].size)))
				idxs_alt    = np.asarray(np.arange(len(ts_tce_alt)), dtype = int)
				tce_idxs_alt= np.array([],dtype = int)
				for i,loc in enumerate(locations_alt): 
					tce_idxs_alt = np.append(tce_idxs_alt, idxs[loc][:per_tce_len_alt[i]])
				idxs_concate	= np.concatenate((tce_idxs, tce_idxs_alt))
				tce_idxs_concate.sort()
				for lag in lags:
					if lag == 0:
						cc,pv = (np.nan,np.nan)
					else:
						if len(tce_idxs_concate) <= len(lags):
							cc,pv = (np.nan,np.nan)
						elif lag<=tce_idxs_concate[0]:
							ts_dri_tce_b  = np.array([cum_av_lagged(df_norm[dri],ignore_t0 = True, lag = lag)[i] for i in (tce_idxs_concate)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv = stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						elif lag <= tce_idxs_concate[int(np.ceil(len(lags)*per_trig_len))]:
							ts_dri_tce_b  = np.array([cum_av_lagged(df_norm[dri],ignore_t0 = True, lag = lag)[i] for i in (tce_idxs_concate)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv = stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						else:
							cc,pv = (np.nan,np.nan)					
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)

				# Case 07: tce_trig_b (normalized)
				# ---------------------------------
				cors = 'tce_trig_b'
				ts_tce	= df_norm['TCE']
				idxs    = np.asarray(np.arange(len(ts_tce)), dtype = int)
				per_trig_len = tce_len_per
				per_tce_len  = []
				for loc in locations:
					per_tce_len . append( int(np.ceil(per_trig_len * bin_tce_1s[loc].size)))
				tce_idxs= np.asarray([],dtype = int)
				for i,loc in enumerate(locations): 
					tce_idxs = np.append(tce_idxs, idxs[loc][:per_tce_len[i]])

				#alt:
				ts_tce_alt	= df_norm['TCE_alt']
				per_tce_len_alt  = []
				for loc in locations_alt:
					per_tce_len_alt . append( int(np.ceil(per_trig_len * bin_tce_1s_alt[loc].size)))

				idxs_alt    = np.asarray(np.arange(len(ts_tce_alt)), dtype = int)
				tce_idxs_alt= np.array([],dtype = int)
				for i,loc in enumerate(locations_alt): 
					tce_idxs_alt = np.append(tce_idxs_alt, idxs[loc][:per_tce_len_alt[i]])

				tce_idxs_concate	= np.concatenate((tce_idxs, tce_idxs_alt))
				tce_idxs_concate.sort()

					
				for lag in lags:
					print ">>>>>>>>>>>>>>>>>> .... ",cors,lag,lat_in,lon_in
					if lag == 0:
						cc,pv = (np.nan,np.nan)
					else:
						if len(tce_idxs_concate) <= len(lags):
							cc,pv = (np.nan,np.nan)
						elif lag<=tce_idxs_concate[0]:
							ts_dri_tce_b  = np.array([np.array(df_norm[dri])[i] for i in (tce_idxs_concate-lag)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv = stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						elif lag <= tce_idxs_concate[int(np.ceil(len(lags)*per_trig_len))]:
							ts_dri_tce_b  = np.array([np.array(df_norm[dri])[i] for i in (tce_idxs_concate-lag)])
							ts_var_tce_b	= np.array([np.array(df_norm[variable])[i] for i in tce_idxs_concate])
							cc,pv = stats.linregress(ts_dri_tce_b, ts_var_tce_b)[2:4]
						else:
							cc,pv = (np.nan,np.nan)					
					df_cc_pv [cors].loc [df_index[lag],dri+'-cc'] = round(cc,4)
					df_cc_pv [cors].loc [df_index[lag],dri+'-pv'] = round(pv,4)

for cors in corr_cases:
	df_cc_pv [cors].to_csv(paths['out']['TEST']+cors+'-cc_pv.csv',sep = ',')

tce_stats = {}
for typ in df_tce.columns:
	tce_stats[typ] = {}
	tce_stats[typ]['mean']  = round(np.array(df_tce[typ]).mean(),3)
	tce_stats[typ]['std']	= round(np.array(df_tce[typ]).std(),3)
	tce_stats[typ]['max']  	= np.array(df_tce[typ]).max()

for i in tce_stats[typ].keys():
	df_tce = df_tce.reindex(df_tce.index.values.tolist()+[i])
for x in tce_stats[typ].keys():
	for y in df_tce.columns:
		df_tce.loc[x,y] = tce_stats[y][x]

df_tce .to_csv(paths['out']['TEST']+conf+'_TCE_stats.csv',sep = ',')
print "Time taken to run this code: ",(timeit.default_timer() - start_time), "seconds"
