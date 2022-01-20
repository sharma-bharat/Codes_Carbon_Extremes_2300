from scipy import stats
from scipy import ndimage
import glob
import sys
import netCDF4 as nc4
import numpy as np
import pandas as pd
import datetime as dt
from calendar import monthrange
import matplotlib as mpl
mpl.use('Agg')
from    mpl_toolkits.basemap import Basemap
from    matplotlib import cm 
import matplotlib.patches as patches

import matplotlib.pyplot as plt
#importing my functions
from functions import time_dim_dates, index_and_dates_slicing, geo_idx,patch_with_gaps_and_eventsize
import argparse

""" Arguments to enter while running the python file
1st argument 	: percentile under consideration
				 looking at the negative/positive  tail of gpp events: {eg.1,5,10,90,95,99}
2nd argument	: 'b' if you want to look at both negative and positive percentiles
3rd argument	: 'pth' if you want to see the plot of the thresholds with time
4th argument	: 'pgts' if you want to see the plot of global timeseries
5th argument	: 'pf' if you want to plot the frequency of anomalies
6th argument	: 'pcr1999' if you want to plot the count of the extremes relative to the thresholds uptill 1999 (from 1850)
7th argument	: integer value of number of years you want the window be to calculate threshold
8th argument	: 'pwin' if you want the window size plots
9th argument	: lat value for which you want to see timeseries [-90, 90]
10th argument	: lon value for which you want to see timeseries [  0,360]
11th argument	: index of the pft type  you want to analyse [ 0 : 16]
				 + '99' if you want all the pft to be analyzed
12th argrument	: 'comp' if you want to compare the results with the overall gpp anomalies
13th argument   : 'compwin' if you want to see the consicutive window plots of pft contribution 
14th argument	: 'productivity type - gpp/ npp/ nep/ nbp 
15th argument 	: 'pct' if you want to use the actual values of the pft percentage and attributing anomalies to pft based on percentage
"""
#p_value/percentile under consideration
# 1, 5 or 10 for percentile looking at the negative tail of gpp events
parser  = argparse.ArgumentParser()
parser.add_argument('--percentile'		,'-per'  		, help = "Threshold Percentile?" 			, type= int, 	default= 1			)
parser.add_argument('--per_both'		,'-p_both'  	, help = "want to look at both percentiles?", type= str, 	default= 'b'		)
parser.add_argument('--plt_threshold'	,'-pl_th'  		, help = "plot threshold with time?" 		, type= str, 	default= 'pth'		)	
parser.add_argument('--plt_global_ts'	,'-pl_gts'  	, help = "plot global timeseries var?" 		, type= str, 	default= 'pgts'		)
parser.add_argument('--plt_freq'		,'-pl_hz'  		, help = "plot frequency ?" 				, type= str, 	default= 'pf'		)
parser.add_argument('--plt_count_r1999'	,'-pc_r1999'	, help = "ext rel to threshold (1850-1999)" , type= str, 	default= 'pcr1999'	)
parser.add_argument('--window'			,'-wsize'		, help = "window size (25 years)?" 			, type= int, 	default= 25			)	
parser.add_argument('--p_win'			,'-p_win'		, help = "window plot sizes?" 				, type= str, 	default= 'pwin'		)
parser.add_argument('--lat_in'			,'-lt_id'		, help = "lat value for ts [-90,90]" 		, type= float, 	default= 999		)
parser.add_argument('--lon_in'			,'-ln_id'		, help = "lon value for ts [0,360]" 		, type= float, 	default= 999		)
parser.add_argument('--pft_idx_in'		,'-pft_idx'		, help = "index of the pft type [0,16]" 	, type= int, 	default= 55			)
parser.add_argument('--pft_comp'		,'-pft_comp'	, help = "compare with the overall gpp ano" , type= str, 	default= 'na'		) # 'comp'
parser.add_argument('--pft_comp_win'	,'-pft_comp_win', help = "compare with gpp ano_cons_win" 	, type= str, 	default= 'na'		) # 'compwin'
parser.add_argument('--variable'		,'-var'			, help = "variable? gpp/npp/nep/nbp,,,,"	, type= str, 	default= 'gpp'		)
parser.add_argument('--model_config'	,'-config'		, help = "w_lulcc or wo_lulcc"				, type= str, 	default= 'misc'	)
parser.add_argument('--save_nc_bin'		,'-nc'			, help = "save nc binary matrix rel th+ wolulcc", type= str, 	default= 'n'	)
parser.add_argument('--save_nc_bin_TCE'	,'-nc_TCE'		, help = "save nc binary TCE matrix rel th+ wolulcc", type= str, 	default= 'n'	)
args = parser.parse_args()

per				= int	(args.percentile)
per_both		= str	(args.per_both)
plt_threshold	= str	(args.plt_threshold)
plt_global_ts	= str	(args.plt_global_ts)
plt_freq		= str	(args.plt_freq)
plt_count_r1999	= str	(args.plt_count_r1999)
window			= int	(args.window)
p_win 	        = str   (args.p_win)
lat_in			= float	(args.lat_in)
lon_in			= float	(args.lon_in)
pft_idx_in		= int	(args.pft_idx_in)
pft_comp		= str	(args.pft_comp)
pft_comp_win	= str	(args.pft_comp_win)
variable		= str 	(args.variable)
save_binary_matrix	= str 	(args.save_nc_bin)
save_TCE_binary	= str 	(args.save_nc_bin_TCE)
#pft_pct		= str	(sys.argv[14])
"""
per             = 1
per_both        = 'b'
plt_threshold   = 'pth'
plt_global_ts   = 'pgts'
plt_freq        = 'pf'
plt_count_r1999 = 'pcr1999'
window          = 25
p_win           = 'pwin'
lat_in          = 999
lon_in          = 999
pft_idx_in      = 4
pft_comp        = 'comp'
pft_comp_win	= 'compwin'
variable		= 'gpp'
#pft_pct		= 'pct'
"""
#running run ecp_extremes_analysis.py -per 1 -p_both b -pl_th pth -pl_gts pgts -pl_hz pf -pc_r1999 pcr1999 -wsize 25 -p_win p -lt_id 999 -ln_id 999 -pft_idx 88 -pft_comp c -pft_comp_win c -var gpp --model_config wo_lulcc --save_nc_bin n --save_nc_bin_TCE n
print (args)

#if 'per_both' ='b' then neg and positive both the percentiles will be analysed
win_len		=	12 * window 			#number of months in window years

# Paths:
# -------
paths                       = {}
paths['in' ]                = {}
paths['in' ]['wo_lulcc']    = '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
paths['in' ]['w_lulcc' ]    = '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'

paths['out' ]				= {}
paths['out' ]['wo_lulcc']	= '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
paths['out' ]['w_lulcc' ]	= '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'
paths['out' ]['comp_lulcc']	= '/home/ud4/CESM1-BGC/'+variable+'/comparison/pftcon_lulcc/'
# Data:
# ------
# 1. The conf are either with or without LULCC
# 2. The original flux (g/m^2/s) is saved in the sub-directory 'var'
# 3. The original flux and anomalies are also saved in units of 'gC'

nc_data 						= {}
nc_data['wo_lulcc']       		= {}
nc_data['w_lulcc' ]         	= {}
nc_data['wo_lulcc']['var']  	= nc4.Dataset(paths['in']['wo_lulcc']   +'cesm1bgc_pftcon_'     +variable+'.nc'     )
nc_data['w_lulcc' ]['var']  	= nc4.Dataset(paths['in']['w_lulcc' ]   +'cesm1bgc_with_lulcc_' +variable+'.nc'     )
nc_data['wo_lulcc']['var_gC' ]  = nc4.Dataset(paths['in']['wo_lulcc']   +'cesm1bgc_pftcon_'     +variable+'_gC.nc'  )
nc_data['w_lulcc' ]['var_gC' ]  = nc4.Dataset(paths['in']['w_lulcc' ]   +'cesm1bgc_with_lulcc_' +variable+'_gC.nc'  )
nc_data['wo_lulcc']['ano_gC' ] 	= nc4.Dataset(paths['in']['wo_lulcc']   +'cesm1bgc_pftcon_'     +variable+'_anomalies_gC.nc')
nc_data['w_lulcc' ]['ano_gC' ]	= nc4.Dataset(paths['in']['w_lulcc' ]   +'cesm1bgc_with_lulcc_' +variable+'_anomalies_gC.nc')

# Universal Variables
# -------------------
lat			= nc_data['wo_lulcc' ]['ano_gC' ].variables['lat']
lon			= nc_data['wo_lulcc' ]['ano_gC' ].variables['lon']
time		= nc_data['wo_lulcc' ]['ano_gC' ].variables['time']
area		= nc_data['wo_lulcc' ]['ano_gC' ].variables['area']
land_frac 	= nc_data['wo_lulcc' ]['var' ]	 .variables['landfrac']

# Arranging Time Array for plotting and calling
# --------------------------------------------
dates_ar	= time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates = np.array([dates_ar[i*win_len] for i in range(int(time.size/win_len))]) 	#list of start dates of 25 year window
end_dates	= np.array([dates_ar[i*win_len+win_len -1] for i in range(int(time.size/win_len))])	#list of end dates of the 25 year window

idx_yr_2100 = 3012 #some the year 1850 if the data is monthly
idx_yr_2300 = 5412 #some the year 1850 if the data is monthly
idx_yr_2299 = 5400 #some the year 1850 if the data is monthly

# Working with the pft types and masks
# -----------------------------------
#	1. The Names of different PFTs indexed from 0 to 16
#	2. ALL PFT is basically the sum of all pfts productivity or it is the original flux; indexed as 99
pft_names		= [	'Bare Ground', 'NET Temperate','NET Boreal',
					'NDT Boreal','BET Tropical','BET Temperate',
					'BDT Tropical','BDT Temperate','BDT Boreal',
					'BES Temperate','BDS Temperate','BDS Boreal',
					'C 3 arctic grass','C 3 grass','C 4 grass',
					'Crop1','Crop2', 'ALL PFT']
pft_idx			= [i for i in range(len(pft_names)-1)]
pft_idx			. append(99)

#later# if pft_idx_in in pft_idx:
#later# 	nc_pft			= nc4.Dataset(in_path+'surfdata_0.9x1.25_urb3den_simyr1850_c090702.nc') 
#later# 	pft				= nc_pft.variables['PCT_PFT']			#pft variable with the different pft types and their percentage in that particular cell
#later# 	sr_pft      	= np.array([i for i in range(pft.shape[0])])
#later# 	pct_pfts    	= []
#later# 	for i in range(pft.shape[0]):
#later# 		tmp     = np.zeros(pft[0,:,:].shape)
#later# 		tmp     = pft[i,:,:]
#later# 		pf      = np.ma.masked_equal(tmp,0)					#mask when the value is 0
#later# 		pct_pfts. append(pf)
#later# 	pft_types	= np.ma.array(pct_pfts)
#later# 
#later# 	pft_sum_pct     = np.sum(pft[1:-1,:,:],axis=0) 	# summing all the pfts expect bareground and crop2 	
#later# 	pft_rel_pct		= pft_types/pft_sum_pct	# rel contritution of pft : pct of pft cover for that pixel divided by total sum of pct for that pixel exlucing bareground and crop2

# Initiation:
# -----------
def TS_Dates_and_Index (dates_ar = dates_ar,start_dates = start_dates, end_dates=end_dates ):
	"""
	Returns the TS of the dates and index of consecutive windows of len 25 years
	
	Parameters:
	-----------
	dates_ar :  an array of dates in datetime.date format 
				the dates are chosen from this array
	start_dates: an array of start dates, the start date will decide the dates and index of the first entry for final time series for that window
	end_dates: similar to start_dates but for end date

	Returns:
	--------
	dates_win: a 2-d array with len of start dates/  total windows and each row containing the dates between start and end date
	idx_dates_win : a 2-d array with len of start dates/  total windows and each row containing the index of dates between start and end date
	"""
	idx_dates_win	= []	#indicies of time in 25yr windows
	dates_win		= []	#sel dates from time variables in win_len windows
	for i in range(len(start_dates)):
		idx_loc, dates_loc 	= index_and_dates_slicing(dates_ar,start_dates[i],end_dates[i]) # see functions.py
		idx_dates_win		. append	(idx_loc)
		dates_win			. append	(dates_loc)
	return np.array(dates_win), np.array(idx_dates_win)

# Calling the function "ts_dates_and_index"; Universal for rest of the code
dates_win, idx_dates_win = TS_Dates_and_Index ()

# The saving the results in a dictionary
# --------------------------------------
Results = {}
Results['wo_lulcc'] = {}
Results['w_lulcc' ] = {}
Results['misc']		= {} # for saving some comparitive results

def Threshold_and_Binary_Ar(data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], per = per ):
	"""
	returns the global percentile based thresholds and binary arrays of consecutive windows
	
	Parameters:
	-----------
	data : The anomalies whose threshold you want to calculate
	
	Universal:
	---------
	start_dates, idx_dates_win, per

	Returns:
	--------
	threshold_neg: 	the threshold for negative extremes; size = # windows
	threshold_pos: 	the threshold for positive extremes; size = # windows
	bin_ext_neg:	the binary array 1's are extremes based on the threshold_neg; shape = same as data	
	bin_ext_pos:	the binary array 1's are extremes based on the threshold_pos; shape = same as data

	"""
	thresholds_1= []	#list of the anomalies threshold for consecutive win_len yr corresponding to 'per'
	thresholds_2= []  	#list of the anomalies threshold for consecutive win_len yr corresponding to '100-per'
	bin_ext_neg	= np.ma.zeros((data.shape)) #3d array to capture the True binaray extmalies w.r.t. gpp loss events
	bin_ext_pos	= np.ma.zeros((data.shape)) #3d array to capture the True binaray extmalies w.r.t. gpp gain events

	for i in range(len(start_dates)):
		ano_loc 		= data[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:]
		threshold_loc_1	= np.percentile(ano_loc[ano_loc.mask == False],per) # calculation of threshold for the local anomalies	
		thresholds_1	. append(threshold_loc_1)
		if plt_global_ts == 'pgts' or plt_freq == 'pf' :
			if per <=50:
				bin_ext_neg[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:] = ano_loc < threshold_loc_1
			else:
				bin_ext_pos[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:] = ano_loc > threshold_loc_1
			
		if per_both		== 'b': #in case we want to analyse both neg and positive percentiles
			threshold_loc_2	= np.percentile(ano_loc[ano_loc.mask == False],(100-per))
			thresholds_2	. append(threshold_loc_2)
			if plt_global_ts == 'pgts' or plt_freq == 'pf' :
				if 100-per <=50:
					bin_ext_neg[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:] = ano_loc	< threshold_loc_2
				else:
					bin_ext_pos[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:] = ano_loc > threshold_loc_2
					#print threshold_loc_2/10**9, bin_ext_pos[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:].sum()

		if per < 50: 
			threshold_neg = np.ma.array(thresholds_1)
			threshold_pos = np.ma.array(thresholds_2)
		elif per > 50:
			threshold_neg = np.ma.array(thresholds_2)
			threshold_pos = np.ma.array(thresholds_1)

	return threshold_neg, threshold_pos, bin_ext_neg, bin_ext_pos
if (plt_threshold == "pth") or (plt_global_ts == "pgts") or (plt_freq == "pf") :
	Results['wo_lulcc']['th_neg'],Results['wo_lulcc']['th_pos'], Results['wo_lulcc']['bin_ext_neg'], Results['wo_lulcc']['bin_ext_pos'] = Threshold_and_Binary_Ar(data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], per = per )
	Results['w_lulcc' ]['th_neg'],Results['w_lulcc' ]['th_pos'], Results['w_lulcc' ]['bin_ext_neg'], Results['w_lulcc' ]['bin_ext_pos'] = Threshold_and_Binary_Ar(data = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...], per = per )

# Based on a positive threshold of wo_lulcc calculation the bin_ext_pos in w_lulcc scenaio : 'bin_ext_pos_w_lulcc_based_pos_wo_lulcc'
# Based on a mean threshold of all scenarios calculation the bin_ext_pos in w_lulcc scenaio : 'bin_ext_pos_w_lulcc_based_av_all'
av_th_tmp		= np.ma.zeros((4,len(start_dates)))
av_th_tmp[0] 	= abs(Results['wo_lulcc']['th_neg'])
av_th_tmp[1] 	= abs(Results['wo_lulcc']['th_pos'])
av_th_tmp[2] 	= abs(Results['w_lulcc' ]['th_neg'])
av_th_tmp[3] 	= abs(Results['w_lulcc' ]['th_pos'])

Results['misc']['av_th_all_sce'] = np.average (av_th_tmp, axis = 0)
del av_th_tmp

#Making a ts of threshold for ploting
Results['wo_lulcc']['ts_th_neg'] = np.ma.array([np.ma.array([Results['wo_lulcc']['th_neg'][i]]*win_len) for i in range(len(Results['wo_lulcc']['th_neg']))]).flatten()
Results['wo_lulcc']['ts_th_pos'] = np.ma.array([np.ma.array([Results['wo_lulcc']['th_pos'][i]]*win_len) for i in range(len(Results['wo_lulcc']['th_pos']))]).flatten()
Results['w_lulcc' ]['ts_th_neg'] = np.ma.array([np.ma.array([Results['w_lulcc' ]['th_neg'][i]]*win_len) for i in range(len(Results['w_lulcc' ]['th_neg']))]).flatten()
Results['w_lulcc' ]['ts_th_pos'] = np.ma.array([np.ma.array([Results['w_lulcc' ]['th_pos'][i]]*win_len) for i in range(len(Results['w_lulcc' ]['th_pos']))]).flatten()

Results['misc']['ts_av_th_all_sce'] = np.ma.array([np.ma.array([Results['misc']['av_th_all_sce'][i]]*win_len) for i in range(len(Results['misc']['av_th_all_sce']))]).flatten()

def Binary_Ar_based_on_Threshold(th, data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'neg' ):
	"""
	returns the binary arrays on consecutive windows based on an input threshold
	
	Parameters:
	-----------
	data : The anomalies whose threshold you want to calculate, the data.shape[0] == len(th)
	th	 : Given threshold; shape = # windows
	type_extreme : [neg, pos], what type of extremes do you want to calculate? 
	if you want to calculate the binary array of negtative extremes based on the input of postive threshold enter 'neg'

	
	Universal:
	---------
	start_dates, idx_dates_win

	Returns:
	--------
	bin_ext:	the binary array 1's are extremes based on the input threshold 'th'; shape = same as data
	"""
	bin_ext	= np.ma.zeros((data.shape)) #3d array to capture the True binaray extmalies based on the input threshold
	
	#to make sure that the supplied threshold is performing calculation correctly i.e. neg threshold should calc negative extremes and vice versa
	if type_extreme == 'neg':
		if np.ma.average(th) < 0: pass
		if np.ma.average(th) > 0: th = -th
	if type_extreme == 'pos':
		if np.ma.average(th) > 0: pass
		if np.ma.average(th) < 0: th = abs(th)
	

	for i in range(len(th)):
		ano_loc 		= data[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:]
		if type_extreme == 'neg':
			bin_ext[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:] = ano_loc < th[i]
		if type_extreme == 'pos':
			bin_ext[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:] = ano_loc > th[i]
	return bin_ext

Results['misc']['bin_ext'] = {}

Results['misc']['bin_ext']['neg_w_lulcc_based_pos_wo_lulcc' ] = Binary_Ar_based_on_Threshold( data = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...], type_extreme = 'neg' , th = Results['wo_lulcc']['th_pos'])
Results['misc']['bin_ext']['pos_w_lulcc_based_pos_wo_lulcc' ] = Binary_Ar_based_on_Threshold( data = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...], type_extreme = 'pos' , th = Results['wo_lulcc']['th_pos'])
Results['misc']['bin_ext']['neg_wo_lulcc_based_pos_wo_lulcc'] = Binary_Ar_based_on_Threshold( data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'neg' , th = Results['wo_lulcc']['th_pos'])
Results['misc']['bin_ext']['pos_wo_lulcc_based_pos_wo_lulcc'] = Binary_Ar_based_on_Threshold( data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'pos' , th = Results['wo_lulcc']['th_pos'])
	
Results['misc']['bin_ext']['neg_w_lulcc_based_av_all' ] 	= Binary_Ar_based_on_Threshold( data = nc_data ['w_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'neg' , th = Results['misc']['av_th_all_sce'])
Results['misc']['bin_ext']['pos_w_lulcc_based_av_all' ] 	= Binary_Ar_based_on_Threshold( data = nc_data ['w_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'pos' , th = Results['misc']['av_th_all_sce'])
Results['misc']['bin_ext']['neg_wo_lulcc_based_av_all'] 	= Binary_Ar_based_on_Threshold( data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'neg' , th = Results['misc']['av_th_all_sce'])
Results['misc']['bin_ext']['pos_wo_lulcc_based_av_all'] 	= Binary_Ar_based_on_Threshold( data = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], type_extreme = 'pos' , th = Results['misc']['av_th_all_sce'])


keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
keys_wrt_av_all		  = ['neg_w_lulcc_based_av_all','pos_w_lulcc_based_av_all','neg_wo_lulcc_based_av_all','pos_wo_lulcc_based_av_all']


if save_binary_matrix in ['y','yy','Y','yes']:
	"""
	If you want to save the binary matrix of extremes as nc files
	this was done so that this coulld be used as input the attribution analysis
	the following code is only for the cases in "keys_wrt_wo_lulcc_pos"
	it can be modified later if needed
	"""
	for case in keys_wrt_wo_lulcc_pos:
		print ("Saving Binary Matrix for the misc case: %s"%case)
		with nc4.Dataset(paths['out' ]['comp_lulcc'] + 'nc_files/bin_ext_'+case+'.nc', mode = 'w') as dset:
			dset        .createDimension( "time",size = time.size)
			dset        .createDimension( "lat" ,size = lat.size)
			dset        .createDimension( "lon" ,size = lon.size)
			t   =   dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time"), fill_value = 1e+36)
			x   =   dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon") , fill_value = 1e+36)
			y   =   dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat") , fill_value = 1e+36)
			z   =   dset.createVariable(varname = variable+'_bin_ext'  ,datatype = float, dimensions = ("time","lat","lon"),fill_value = 1e+36) #varible = gpp_bin_ext
			t.axis  =   "T"
			x.axis  =   "X"
			y.axis  =   "Y"
			t[...]  =   time  [...]
			x[...]  =   lon   [...]
			y[...]  =   lat   [...]
			z[...]  =   Results['misc']['bin_ext'][case]
			z.missing_value = 1e+36
			z.stardard_name = variable+"binary extremes matrix"
			z.units         = "0,1"
			x.units         =   lat.units
			x.missing_value =   1e+36
			x.setncattr         ("long_name",lat.long_name)
			y.units         =   lon.units
			y.missing_value =   1e+36
			y.setncattr         ("long_name",lon.long_name)
			t.units         =   time.units
			t.setncattr         ("calendar",time.calendar)
			t.setncattr         ("long_name",time.long_name)
			t.missing_value =   1e+36

# For TCE in misc case lags are used for print bin mat
lags_TCE = np.asarray([0,1,2,3,4], dtype = int)
def Binary_Mat_TCE_Win (bin_ar, win_start_year=2000,lags = lags_TCE):
	"""
	To save the binary matrix of the Time Continuous Extremes(TCEs) so that the location and duration of the extremes can be identified.

	returns:
bin_TCE_01s: are the binary values of extreme values in a TCE only at qualified locations with gaps ( actual as value 0) [hightlight extreme values]
	bin_TCE_1s : are the binary values of extreme values in a TCE only at qualified locations with gaps ( 0 replaced with value 1) [selecting full TCE with only 1s]
	bin_TCE_len : are the len of TCE extreme events, the length  of TCE is captured at the trigger locations
	shape : These matrix are of shape (5,300,192,288) i.e. lags(0-4 months), time(300 months or 25 years {2000-24}), lat(192) and lon(288).

	"""
	from functions import create_seq_mat

	for i,date in enumerate(start_dates):
		if date.year in [win_start_year]:
			start_yr_idx = i
	
	data = bin_ar[start_yr_idx*win_len: (start_yr_idx+1)*win_len]
	del bin_ar
	
	
	bin_TCE_1s 	= np.ma.zeros((len(lags), data.shape[0],data.shape[1],data.shape[2]))
	bin_TCE_01s = np.ma.zeros((len(lags), data.shape[0],data.shape[1],data.shape[2]))
	bin_TCE_len = np.ma.zeros((len(lags), data.shape[0],data.shape[1],data.shape[2]))
	
	for lag in lags:
		for lat_i in range(len(lat)):
			for lon_i in range(len(lon)):
				if land_frac[...].mask[lat_i,lon_i] != True:
					#print lag, lat_i, lon_i
					try:
						tmp = patch_with_gaps_and_eventsize (data[:,lat_i,lon_i], max_gap =2, min_cont_event_size=3, lag=lag)
						for idx, trig in enumerate (tmp[1]):
							bin_TCE_01s [lag, trig:trig+len(tmp[0][idx]), lat_i, lon_i] = tmp[0][idx] 
							bin_TCE_1s  [lag, trig:trig+len(tmp[0][idx]), lat_i, lon_i] = np.ones(tmp[0][idx].shape)
							bin_TCE_len [lag, trig, lat_i, lon_i] 						= np.sum(np.ones(tmp[0][idx].shape))
					except:
						bin_TCE_01s[lag, :, lat_i, lon_i]  = np.ma.masked_all(data.shape[0])
						bin_TCE_1s [lag, :, lat_i, lon_i]  = np.ma.masked_all(data.shape[0])
						bin_TCE_len[lag, :, lat_i, lon_i]  = np.ma.masked_all(data.shape[0])

				else:
					bin_TCE_01s[lag, :, lat_i, lon_i]  = np.ma.masked_all(data.shape[0])
					bin_TCE_1s [lag, :, lat_i, lon_i]  = np.ma.masked_all(data.shape[0])
					bin_TCE_len[lag, :, lat_i, lon_i]  = np.ma.masked_all(data.shape[0])

#					try:
#						sp_count_TCE[lat_i,lon_i]=len(patch_with_gaps_and_eventsize(bin_ar = data[:,lat_i,lon_i], max_gap =2,min_cont_event_size =3, lag=lag)[0])
#					except:
##						sp_count_TCE[lat_i,lon_i] = 0
	return bin_TCE_01s, bin_TCE_1s, bin_TCE_len

win_start_years = np.arange(1850,2300,25)

if save_TCE_binary in ['y','yy','Y','yes']:
	"""
	To save the binary matrix of the Time Continuous Extremes(TCEs) so that the location and duration of the extremes can be identified.

	If you want to save the binary matrix of extremes as nc files
	this was done so that this coulld be used as input the attribution analysis
	the following code is only for the cases in "keys_wrt_wo_lulcc_pos"
	it can be modified later if needed
	"""
	for start_yr in win_start_years:
		Results['misc']['bin_TCE_01s']  = {}
		Results['misc']['bin_TCE_1s' ]  = {}
		Results['misc']['bin_TCE_len']  = {}
		for key in keys_wrt_wo_lulcc_pos:
			Results['misc']['bin_TCE_01s'][key], Results['misc']['bin_TCE_1s'][key], Results['misc']['bin_TCE_len'][key] = Binary_Mat_TCE_Win (Results['misc']['bin_ext'][key], win_start_year = start_yr)

		for case in keys_wrt_wo_lulcc_pos:
			print ("Saving Binary Matrix for the misc 01s case: %s"%case)
			with nc4.Dataset(paths['out' ]['comp_lulcc'] + 'nc_files/bin_TCE_01s_'+case+'_%d.nc'%start_yr, mode = 'w') as dset:
				dset        .createDimension( "lag",size = lags_TCE.size)
				dset        .createDimension( "time",size = win_len)
				dset        .createDimension( "lat" ,size = lat.size)
				dset        .createDimension( "lon" ,size = lon.size)
				w   =   dset.createVariable(varname = "lag"  ,datatype = float, dimensions = ("lag") , fill_value = 1e+36)
				t   =   dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time"), fill_value = 1e+36)
				x   =   dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon") , fill_value = 1e+36)
				y   =   dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat") , fill_value = 1e+36)
				z   =   dset.createVariable(varname = variable+'_TCE_01s'  ,datatype = float, dimensions = ("lag","time","lat","lon"),fill_value = 1e+36) #varible = gpp_bin_ext
				w.axis  =   "T"
				t.axis  =   "T"
				x.axis  =   "X"
				y.axis  =   "Y"
				w[...]  =   lags_TCE
				t[...]  =   time  [...][6*win_len : 7*win_len]
				x[...]  =   lon   [...]
				y[...]  =   lat   [...]
				z[...]  =   Results['misc']['bin_TCE_01s'][case]
				z.missing_value = 1e+36
				z.stardard_name = variable+" binary TCE (01s) matrix for 25 years starting at the  year %d"%start_yr
				z.units         = "0,1"
				x.units         =   lat.units
				x.missing_value =   1e+36
				x.setncattr         ("long_name",lat.long_name)
				y.units         =   lon.units
				y.missing_value =   1e+36
				y.setncattr         ("long_name",lon.long_name)
				t.units         =   time.units
				t.setncattr         ("calendar",time.calendar)
				t.setncattr         ("long_name",time.long_name)
				t.missing_value =   1e+36
				w.units         =   "month"
				w.setncattr         ("long_name","lags in months")
				w.missing_value =   1e+36

			print ("Saving Binary Matrix for the misc TCE_1s case: %s"%case)
			with nc4.Dataset(paths['out' ]['comp_lulcc'] + 'nc_files/bin_TCE_1s_'+case+'_%d.nc'%start_yr, mode = 'w') as dset:
				dset        .createDimension( "lag",size = lags_TCE.size)
				dset        .createDimension( "time",size = win_len)
				dset        .createDimension( "lat" ,size = lat.size)
				dset        .createDimension( "lon" ,size = lon.size)
				w   =   dset.createVariable(varname = "lag"  ,datatype = float, dimensions = ("lag") , fill_value = 1e+36)
				t   =   dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time"), fill_value = 1e+36)
				x   =   dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon") , fill_value = 1e+36)
				y   =   dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat") , fill_value = 1e+36)
				z   =   dset.createVariable(varname = variable+'_TCE_1s'  ,datatype = float, dimensions = ("lag","time","lat","lon"),fill_value = 1e+36) #varible = gpp_bin_ext
				w.axis  =   "T"
				t.axis  =   "T"
				x.axis  =   "X"
				y.axis  =   "Y"
				w[...]  =   lags_TCE
				t[...]  =   time  [...][6*win_len:7*win_len]
				x[...]  =   lon   [...]
				y[...]  =   lat   [...]
				z[...]  =   Results['misc']['bin_TCE_1s'][case]
				z.missing_value = 1e+36
				z.stardard_name = variable+"binary TCE (1s) matrix for 25 years starting at the  year %d"%start_yr
				z.units         = "0,1"
				x.units         =   lat.units
				x.missing_value =   1e+36
				x.setncattr         ("long_name",lat.long_name)
				y.units         =   lon.units
				y.missing_value =   1e+36
				y.setncattr         ("long_name",lon.long_name)
				t.units         =   time.units
				t.setncattr         ("calendar",time.calendar)
				t.setncattr         ("long_name",time.long_name)
				t.missing_value =   1e+36
				w.units         =   "month"
				w.setncattr         ("long_name","lags in months")
				w.missing_value =   1e+36

			print ("Saving Binary Matrix for the misc TCE_len case: %s"%case)
			with nc4.Dataset(paths['out' ]['comp_lulcc'] + 'nc_files/bin_TCE_len_'+case+'_%d.nc'%start_yr, mode = 'w') as dset:
				dset        .createDimension( "lag",size = lags_TCE.size)
				dset        .createDimension( "time",size = win_len)
				dset        .createDimension( "lat" ,size = lat.size)
				dset        .createDimension( "lon" ,size = lon.size)
				w   =   dset.createVariable(varname = "lag"  ,datatype = float, dimensions = ("lag") , fill_value = 1e+36)
				t   =   dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time"), fill_value = 1e+36)
				x   =   dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon") , fill_value = 1e+36)
				y   =   dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat") , fill_value = 1e+36)
				z   =   dset.createVariable(varname = variable+'_TCE_len'  ,datatype = float, dimensions = ("lag","time","lat","lon"),fill_value = 1e+36) #varible = gpp_bin_ext
				w.axis  =   "T"
				t.axis  =   "T"
				x.axis  =   "X"
				y.axis  =   "Y"
				w[...]  =   lags_TCE
				t[...]  =   time  [...][6*win_len:7*win_len]
				x[...]  =   lon   [...]
				y[...]  =   lat   [...]
				z[...]  =   Results['misc']['bin_TCE_len'][case]
				z.missing_value = 1e+36
				z.stardard_name = variable+"TCE length matrix 25 years starting at the  year %d"%start_yr
				z.units         = "0,1"
				x.units         =   lat.units
				x.missing_value =   1e+36
				x.setncattr         ("long_name",lat.long_name)
				y.units         =   lon.units
				y.missing_value =   1e+36
				y.setncattr         ("long_name",lon.long_name)
				t.units         =   time.units
				t.setncattr         ("calendar",time.calendar)
				t.setncattr         ("long_name",time.long_name)
				t.missing_value =   1e+36
				w.units         =   "month"
				w.setncattr         ("long_name","lags in months")
				w.missing_value =   1e+36


def Spatial_Freq_Extremes (bin_ar):
	"""
	returns the spatial frequency of extreme events; shape (# wins, nlat, nlon)

	Parameters:
	-----------
	bin_ar : a binary array of shape (# wins, #)
	
	Method:
	------
	it loops over the bin_ar.shape[0] and sums all 1's over axis = 0; calculating the spatial freq per time window
	"""
	return np.ma.array([np.ma.sum(bin_ar[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:],axis=0) for i in range(idx_dates_win.shape[0])])

Results['wo_lulcc']['sp_freq_neg'] = Spatial_Freq_Extremes(Results['wo_lulcc']['bin_ext_neg'])
Results['wo_lulcc']['sp_freq_pos'] = Spatial_Freq_Extremes(Results['wo_lulcc']['bin_ext_pos'])
Results['w_lulcc' ]['sp_freq_neg'] = Spatial_Freq_Extremes(Results['w_lulcc' ]['bin_ext_neg'])
Results['w_lulcc' ]['sp_freq_pos'] = Spatial_Freq_Extremes(Results['w_lulcc' ]['bin_ext_pos'])

Results['misc']['sp_freq'] = {}
for key in Results['misc']['bin_ext'].keys():
	Results['misc']['sp_freq'][key] = Spatial_Freq_Extremes (Results['misc']['bin_ext'][key])

def Spatial_Freq_TCE_Extremes_Win (bin_ar, win_start_year=2000):
	"""
	returns the spatial frequency of TCE extreme events (ref functions.py); shape ( nlat, nlon)

	Parameters:
	-----------
	bin_ar : a binary array of shape (win_len,nlat,nlon)
	
	Method:
	------
	ref: patch_with_gaps_and_eventsize in 'functions.py'
	"""
	
	for i,date in enumerate(start_dates):
		if date.year in [win_start_year]:
			start_yr_idx = i
	
	data = bin_ar[start_yr_idx*win_len: (start_yr_idx+1)*win_len]
	del bin_ar
	sp_count_TCE = np.ma.zeros(data.mask[0].shape)
	sp_count_TCE.mask = data.mask[0]
	for lat_i in range(len(lat)):
		for lon_i in range(len(lon)):
			if sp_count_TCE.mask[lat_i,lon_i] != True:
				try:
					sp_count_TCE[lat_i,lon_i]=len(patch_with_gaps_and_eventsize(bin_ar = data[:,lat_i,lon_i], max_gap =2,min_cont_event_size =3, lag=None)[0])
				except:
					sp_count_TCE[lat_i,lon_i] = 0

	return np.asarray(sp_count_TCE, dtype = int)

Results['wo_lulcc']['sp_freq_neg_TCE'] = Spatial_Freq_TCE_Extremes_Win(Results['wo_lulcc']['bin_ext_neg'])
Results['wo_lulcc']['sp_freq_pos_TCE'] = Spatial_Freq_TCE_Extremes_Win(Results['wo_lulcc']['bin_ext_pos'])
Results['w_lulcc' ]['sp_freq_neg_TCE'] = Spatial_Freq_TCE_Extremes_Win(Results['w_lulcc' ]['bin_ext_neg'])
Results['w_lulcc' ]['sp_freq_pos_TCE'] = Spatial_Freq_TCE_Extremes_Win(Results['w_lulcc' ]['bin_ext_pos'])

Results['misc']['sp_freq_TCE'] = {}
for key in Results['misc']['bin_ext'].keys():
	Results['misc']['sp_freq_TCE'][key] = Spatial_Freq_TCE_Extremes_Win (Results['misc']['bin_ext'][key])

def Global_TS_of_Extremes(bin_ar, ano_gC):
	"""
	Returns the global TS of :
	1. total carbon loss/gain associated neg/pos extremes	
	2. total freq of extremes
	3. total area affected by extremes

	Parameters:
	-----------
	bin_ar : the binary array of extremes (pos/neg)
	ano_gC : the array which will use the mask or binary arrays to calc the carbon loss/gain

	Universal:
	----------
	2-d area array (nlat, nlon), dates_win (# wins, win_size)


	Returns:
	--------
	1d array of length # wins x win_size for all : ext_gC_ts, ext_freq_ts, ext_area_ts
	"""
	ext_ar		= bin_ar * ano_gC 	 # extremes array
	ext_area_ar	= bin_ar * ( area[...] * land_frac[...]) # area array of extremes
	ext_gC_ts	= []
	ext_freq_ts = []
	ext_area_ts = []
	for i in range(dates_win.flatten().size):
		ext_gC_ts	. append(np.ma.sum(ext_ar[i]))
		ext_freq_ts	. append(np.ma.sum(bin_ar[i]))
		ext_area_ts	. append(np.ma.sum(ext_area_ar[i]))
	
	return np.ma.array(ext_gC_ts), np.ma.array(ext_freq_ts),np.ma.array(ext_area_ts)
	
def Slope_Intercept_Pv_Trend_Increase (ts):
	"""
	Returns the slope, intercept, r value , p value and trend line points for time period 1850-2100 (as '_21') and 2101-2300 ('_23')

	Parameters:
	-----------
	One dimentional time series of len 5400 from 1850 through 2299

	Returns:
	--------
	single values for slope, intercept, r value , p value, increase percentage**
	1d array for same legnth as 'ts' for 'trend'

	** it return the percent increase of trend line relavtive to the year 1850 (mean trend line value),..
	"""
	# calculation of the magnitudes of global gpp loss and trend from 1850-2100
 	slope_21, intercept_21,rv_21,pv_21,std_e21  = stats.linregress(time[...][:idx_yr_2100],ts[:idx_yr_2100])
	trend_21       								= slope_21*time[...][:idx_yr_2100]+intercept_21
	increase_21    								= (trend_21[-1]-trend_21[0])*100/trend_21[0]

	# calculation of the magnitudes of global gpp loss and trend from 2101-2300
 	slope_23, intercept_23,rv_23,pv_23,std_e23  = stats.linregress(time[...][idx_yr_2100:idx_yr_2299],ts[idx_yr_2100:idx_yr_2299])
	trend_23       								= slope_23*time[...][idx_yr_2100:idx_yr_2299]+intercept_23
	increase_23    								= (trend_23[-1]-trend_23[0])*100/trend_23[0]
	increase_23_r1850							= (trend_23[-1]-trend_21[0])*100/trend_21[0]

	return slope_21,intercept_21,pv_21,trend_21,increase_21,slope_23,intercept_23,pv_23,trend_23,increase_23,increase_23_r1850

Results['wo_lulcc']['ts_global_gC'] = {}
Results['w_lulcc' ]['ts_global_gC'] = {}
Results['wo_lulcc']['ts_global_freq'] = {}
Results['w_lulcc' ]['ts_global_freq'] = {}
Results['wo_lulcc']['ts_global_area'] = {}
Results['w_lulcc' ]['ts_global_area'] = {}

Results['wo_lulcc']['ts_global_gC']['neg_ext']={}
Results['wo_lulcc']['ts_global_gC']['pos_ext']={}
Results['wo_lulcc']['ts_global_freq']['neg_ext']={}
Results['wo_lulcc']['ts_global_freq']['pos_ext']={}
Results['wo_lulcc']['ts_global_area']['neg_ext']={}
Results['wo_lulcc']['ts_global_area']['pos_ext']={}
Results['w_lulcc']['ts_global_gC']['neg_ext']={}
Results['w_lulcc']['ts_global_gC']['pos_ext']={}
Results['w_lulcc']['ts_global_freq']['neg_ext']={}
Results['w_lulcc']['ts_global_freq']['pos_ext']={}
Results['w_lulcc']['ts_global_area']['neg_ext']={}
Results['w_lulcc']['ts_global_area']['pos_ext']={}

Results['wo_lulcc']['ts_global_gC']['neg_ext']['ts'], Results['wo_lulcc']['ts_global_freq']['neg_ext']['ts'],Results['wo_lulcc']['ts_global_area']['neg_ext']['ts'] = Global_TS_of_Extremes ( bin_ar = Results['wo_lulcc']['bin_ext_neg'], ano_gC = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...])
Results['wo_lulcc']['ts_global_gC']['pos_ext']['ts'], Results['wo_lulcc']['ts_global_freq']['pos_ext']['ts'],Results['wo_lulcc']['ts_global_area']['pos_ext']['ts'] = Global_TS_of_Extremes ( bin_ar = Results['wo_lulcc']['bin_ext_pos'], ano_gC = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...])
Results['w_lulcc' ]['ts_global_gC']['neg_ext']['ts'], Results['w_lulcc' ]['ts_global_freq']['neg_ext']['ts'],Results['w_lulcc' ]['ts_global_area']['neg_ext']['ts'] = Global_TS_of_Extremes ( bin_ar = Results['w_lulcc' ]['bin_ext_neg'], ano_gC = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...])
Results['w_lulcc' ]['ts_global_gC']['pos_ext']['ts'], Results['w_lulcc' ]['ts_global_freq']['pos_ext']['ts'],Results['w_lulcc' ]['ts_global_area']['pos_ext']['ts'] = Global_TS_of_Extremes ( bin_ar = Results['w_lulcc' ]['bin_ext_pos'], ano_gC = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...])

Results['wo_lulcc']['ts_global_gC']['neg_ext']['s21'],_,Results['wo_lulcc']['ts_global_gC']['neg_ext']['pv21'],Results['wo_lulcc']['ts_global_gC']['neg_ext']['trend_21'],Results['wo_lulcc']['ts_global_gC']['neg_ext']['inc_21'],Results['wo_lulcc']['ts_global_gC']['neg_ext']['s23'],_,Results['wo_lulcc']['ts_global_gC']['neg_ext']['pv_23'],Results['wo_lulcc']['ts_global_gC']['neg_ext']['trend_23'],Results['wo_lulcc']['ts_global_gC']['neg_ext']['inc_23'], Results['wo_lulcc']['ts_global_gC']['neg_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['wo_lulcc']['ts_global_gC']['neg_ext']['ts']) 
Results['wo_lulcc']['ts_global_gC']['pos_ext']['s21'],_,Results['wo_lulcc']['ts_global_gC']['pos_ext']['pv21'],Results['wo_lulcc']['ts_global_gC']['pos_ext']['trend_21'], Results['wo_lulcc']['ts_global_gC']['pos_ext']['inc_21'] , Results['wo_lulcc']['ts_global_gC']['pos_ext']['s23'],_,Results['wo_lulcc']['ts_global_gC']['pos_ext']['pv_23'],Results['wo_lulcc']['ts_global_gC']['pos_ext']['trend_23'], Results['wo_lulcc']['ts_global_gC']['pos_ext']['inc_23'],Results['wo_lulcc']['ts_global_gC']['pos_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['wo_lulcc']['ts_global_gC']['pos_ext']['ts']) 
Results['wo_lulcc']['ts_global_freq']['neg_ext']['s21'],_,Results['wo_lulcc']['ts_global_freq']['neg_ext']['pv21'],Results['wo_lulcc']['ts_global_freq']['neg_ext']['trend_21'], Results['wo_lulcc']['ts_global_freq']['neg_ext']['inc_21'], Results['wo_lulcc']['ts_global_freq']['neg_ext']['s23'],_,Results['wo_lulcc']['ts_global_freq']['neg_ext']['pv_23'],Results['wo_lulcc']['ts_global_freq']['neg_ext']['trend_23'], Results['wo_lulcc']['ts_global_freq']['neg_ext']['inc_23'],Results['wo_lulcc']['ts_global_freq']['neg_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['wo_lulcc']['ts_global_freq']['neg_ext']['ts']) 
Results['wo_lulcc']['ts_global_freq']['pos_ext']['s21'],_,Results['wo_lulcc']['ts_global_freq']['pos_ext']['pv21'],Results['wo_lulcc']['ts_global_freq']['pos_ext']['trend_21'], Results['wo_lulcc']['ts_global_freq']['pos_ext']['inc_21'],Results['wo_lulcc']['ts_global_freq']['pos_ext']['s23'],_,Results['wo_lulcc']['ts_global_freq']['pos_ext']['pv_23'],Results['wo_lulcc']['ts_global_freq']['pos_ext']['trend_23'],Results['wo_lulcc']['ts_global_freq']['pos_ext']['inc_23'], Results['wo_lulcc']['ts_global_freq']['pos_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['wo_lulcc']['ts_global_freq']['pos_ext']['ts']) 
Results['wo_lulcc']['ts_global_area']['neg_ext']['s21'],_,Results['wo_lulcc']['ts_global_area']['neg_ext']['pv21'],Results['wo_lulcc']['ts_global_area']['neg_ext']['trend_21'],Results['wo_lulcc']['ts_global_area']['neg_ext']['inc_21'],Results['wo_lulcc']['ts_global_area']['neg_ext']['s23'],_,Results['wo_lulcc']['ts_global_area']['neg_ext']['pv_23'],Results['wo_lulcc']['ts_global_area']['neg_ext']['trend_23'], Results['wo_lulcc']['ts_global_area']['neg_ext']['inc_23'], Results['wo_lulcc']['ts_global_area']['neg_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['wo_lulcc']['ts_global_area']['neg_ext']['ts']) 
Results['wo_lulcc']['ts_global_area']['pos_ext']['s21'],_,Results['wo_lulcc']['ts_global_area']['pos_ext']['pv21'],Results['wo_lulcc']['ts_global_area']['pos_ext']['trend_21'], Results['wo_lulcc']['ts_global_area']['pos_ext']['inc_21'],Results['wo_lulcc']['ts_global_area']['pos_ext']['s23'],_,Results['wo_lulcc']['ts_global_area']['pos_ext']['pv_23'],Results['wo_lulcc']['ts_global_area']['pos_ext']['trend_23'], Results['wo_lulcc']['ts_global_area']['pos_ext']['inc_23'],Results['wo_lulcc']['ts_global_area']['pos_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['wo_lulcc']['ts_global_area']['pos_ext']['ts']) 

Results['w_lulcc']['ts_global_gC']['neg_ext']['s21'],_,Results['w_lulcc']['ts_global_gC']['neg_ext']['pv21'],Results['w_lulcc']['ts_global_gC']['neg_ext']['trend_21'], Results['w_lulcc']['ts_global_gC']['neg_ext']['inc_21'], Results['w_lulcc']['ts_global_gC']['neg_ext']['s23'],_,Results['w_lulcc']['ts_global_gC']['neg_ext']['pv_23'],Results['w_lulcc']['ts_global_gC']['neg_ext']['trend_23'], Results['w_lulcc']['ts_global_gC']['neg_ext']['inc_23'], Results['w_lulcc']['ts_global_gC']['neg_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['w_lulcc']['ts_global_gC']['neg_ext']['ts']) 
Results['w_lulcc']['ts_global_gC']['pos_ext']['s21'],_,Results['w_lulcc']['ts_global_gC']['pos_ext']['pv21'],Results['w_lulcc']['ts_global_gC']['pos_ext']['trend_21'], Results['w_lulcc']['ts_global_gC']['pos_ext']['inc_21'], Results['w_lulcc']['ts_global_gC']['pos_ext']['s23'],_,Results['w_lulcc']['ts_global_gC']['pos_ext']['pv_23'],Results['w_lulcc']['ts_global_gC']['pos_ext']['trend_23'], Results['w_lulcc']['ts_global_gC']['pos_ext']['inc_23'], Results['w_lulcc']['ts_global_gC']['pos_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['w_lulcc']['ts_global_gC']['pos_ext']['ts']) 
Results['w_lulcc']['ts_global_freq']['neg_ext']['s21'],_,Results['w_lulcc']['ts_global_freq']['neg_ext']['pv21'],Results['w_lulcc']['ts_global_freq']['neg_ext']['trend_21'], Results['w_lulcc']['ts_global_freq']['neg_ext']['inc_21'],Results['w_lulcc']['ts_global_freq']['neg_ext']['s23'],_,Results['w_lulcc']['ts_global_freq']['neg_ext']['pv_23'],Results['w_lulcc']['ts_global_freq']['neg_ext']['trend_23'], Results['w_lulcc']['ts_global_freq']['neg_ext']['inc_23'],  Results['w_lulcc']['ts_global_freq']['neg_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['w_lulcc']['ts_global_freq']['neg_ext']['ts']) 
Results['w_lulcc']['ts_global_freq']['pos_ext']['s21'],_,Results['w_lulcc']['ts_global_freq']['pos_ext']['pv21'],Results['w_lulcc']['ts_global_freq']['pos_ext']['trend_21'], Results['w_lulcc']['ts_global_freq']['pos_ext']['inc_21'],Results['w_lulcc']['ts_global_freq']['pos_ext']['s23'],_,Results['w_lulcc']['ts_global_freq']['pos_ext']['pv_23'],Results['w_lulcc']['ts_global_freq']['pos_ext']['trend_23'], Results['w_lulcc']['ts_global_freq']['pos_ext']['inc_23'],Results['w_lulcc']['ts_global_freq']['pos_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['w_lulcc']['ts_global_freq']['pos_ext']['ts']) 
Results['w_lulcc']['ts_global_area']['neg_ext']['s21'],_,Results['w_lulcc']['ts_global_area']['neg_ext']['pv21'],Results['w_lulcc']['ts_global_area']['neg_ext']['trend_21'], Results['w_lulcc']['ts_global_area']['neg_ext']['inc_21'],Results['w_lulcc']['ts_global_area']['neg_ext']['s23'],_,Results['w_lulcc']['ts_global_area']['neg_ext']['pv_23'],Results['w_lulcc']['ts_global_area']['neg_ext']['trend_23'], Results['w_lulcc']['ts_global_area']['neg_ext']['inc_23'], Results['w_lulcc']['ts_global_area']['neg_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['w_lulcc']['ts_global_area']['neg_ext']['ts']) 
Results['w_lulcc']['ts_global_area']['pos_ext']['s21'],_,Results['w_lulcc']['ts_global_area']['pos_ext']['pv21'],Results['w_lulcc']['ts_global_area']['pos_ext']['trend_21'],Results['w_lulcc']['ts_global_area']['pos_ext']['inc_21'],Results['w_lulcc']['ts_global_area']['pos_ext']['s23'],_,Results['w_lulcc']['ts_global_area']['pos_ext']['pv_23'],Results['w_lulcc']['ts_global_area']['pos_ext']['trend_23'], Results['w_lulcc']['ts_global_area']['pos_ext']['inc_23'], Results['w_lulcc']['ts_global_area']['pos_ext']['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['w_lulcc']['ts_global_area']['pos_ext']['ts']) 

#Misc
Results['misc']['ts_global_gC'] 	= {}
Results['misc']['ts_global_freq'] 	= {}
Results['misc']['ts_global_area'] 	= {}

keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
keys_wrt_av_all		  = ['neg_w_lulcc_based_av_all','pos_w_lulcc_based_av_all','neg_wo_lulcc_based_av_all','pos_wo_lulcc_based_av_all']
for key in keys_wrt_wo_lulcc_pos:
	Results['misc']['ts_global_gC'][key]   	= {}
	Results['misc']['ts_global_freq'][key] 	= {}
	Results['misc']['ts_global_area'][key] 	= {}
	Results['misc']['ts_global_gC'][key]['ts'],Results['misc']['ts_global_freq'][key]['ts'], Results['misc']['ts_global_area'][key]['ts'] = Global_TS_of_Extremes ( bin_ar = Results['misc']['bin_ext'][key ], ano_gC = nc_data[key.split('_')[1]+'_'+key.split('_')[2]]['ano_gC' ].variables[variable][...])

for key in keys_wrt_wo_lulcc_pos:
	Results['misc']['ts_global_gC'][key]['s21'],_,Results['misc']['ts_global_gC'][key]['pv_21'],Results['misc']['ts_global_gC'][key]['trend_21'], Results['misc']['ts_global_gC'][key]['inc_21'],Results['misc']['ts_global_gC'][key]['s23'],_,Results['misc']['ts_global_gC'][key]['pv_23'],Results['misc']['ts_global_gC'][key]['trend_23'],Results['misc']['ts_global_gC'][key]['inc_23'], Results['misc']['ts_global_gC'][key]['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['misc']['ts_global_gC'][key]['ts'])
	Results['misc']['ts_global_freq'][key]['s21'],_,Results['misc']['ts_global_freq'][key]['pv_21'],Results['misc']['ts_global_freq'][key]['trend_21'], Results['misc']['ts_global_freq'][key]['inc_21'],Results['misc']['ts_global_freq'][key]['s23'],_,Results['misc']['ts_global_freq'][key]['pv_23'],Results['misc']['ts_global_freq'][key]['trend_23'], Results['misc']['ts_global_freq'][key]['inc_23'], Results['misc']['ts_global_freq'][key]['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['misc']['ts_global_freq'][key]['ts'])
	Results['misc']['ts_global_area'][key]['s21'],_,Results['misc']['ts_global_area'][key]['pv_21'],Results['misc']['ts_global_area'][key]['trend_21'], Results['misc']['ts_global_area'][key]['inc_21'] , Results['misc']['ts_global_area'][key]['s23'],_,Results['misc']['ts_global_area'][key]['pv_23'],Results['misc']['ts_global_area'][key]['trend_23'] , Results['misc']['ts_global_area'][key]['inc_23'], Results['misc']['ts_global_area'][key]['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['misc']['ts_global_area'][key]['ts'])

for key in keys_wrt_av_all:
	Results['misc']['ts_global_gC'][key] 	= {}
	Results['misc']['ts_global_freq'][key] 	= {}
	Results['misc']['ts_global_area'][key] 	= {}
	Results['misc']['ts_global_gC'][key]['ts'],Results['misc']['ts_global_freq'][key]['ts'], Results['misc']['ts_global_area'][key]['ts'] = Global_TS_of_Extremes ( bin_ar = Results['misc']['bin_ext'][key ], ano_gC = nc_data[key.split('_')[1]+'_'+key.split('_')[2]]['ano_gC' ].variables[variable][...])

for key in keys_wrt_av_all:
	Results['misc']['ts_global_gC'][key]['s21'],_,Results['misc']['ts_global_gC'][key]['pv_21'],Results['misc']['ts_global_gC'][key]['trend_21'], Results['misc']['ts_global_gC'][key]['inc_21'],Results['misc']['ts_global_gC'][key]['s23'],_,Results['misc']['ts_global_gC'][key]['pv_23'],Results['misc']['ts_global_gC'][key]['trend_23'], Results['misc']['ts_global_gC'][key]['inc_23'], Results['misc']['ts_global_gC'][key]['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['misc']['ts_global_gC'][key]['ts'])
	Results['misc']['ts_global_freq'][key]['s21'],_,Results['misc']['ts_global_freq'][key]['pv_21'],Results['misc']['ts_global_freq'][key]['trend_21'], Results['misc']['ts_global_freq'][key]['inc_21'],Results['misc']['ts_global_freq'][key]['s23'],_,Results['misc']['ts_global_freq'][key]['pv_23'],Results['misc']['ts_global_freq'][key]['trend_23'],Results['misc']['ts_global_freq'][key]['inc_23'], Results['misc']['ts_global_freq'][key]['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['misc']['ts_global_freq'][key]['ts'])
	Results['misc']['ts_global_area'][key]['s21'],_,Results['misc']['ts_global_area'][key]['pv_21'],Results['misc']['ts_global_area'][key]['trend_21'], Results['misc']['ts_global_area'][key]['inc_21'],Results['misc']['ts_global_area'][key]['s23'],_,Results['misc']['ts_global_area'][key]['pv_23'],Results['misc']['ts_global_area'][key]['trend_23'],Results['misc']['ts_global_area'][key]['inc_23'],Results['misc']['ts_global_area'][key]['inc_23_r18'] = Slope_Intercept_Pv_Trend_Increase(Results['misc']['ts_global_area'][key]['ts'])

def Sum_and_Diff_of_Fluxes_perWin(ano_gC, bin_ar = None, data_type = 'ext', diff_ref_yr = 1975):
	"""
	returns a 2-d array sum of fluxes and difference of the sum of fluxes with reference to the ref yr

	Parameters:
	----------
	bin_ar: the binary array of extremes (pos/neg)
	ano_gC : the array which will use the mask or binary arrays to calc the carbon loss/gain
	diff_ref_yr : the starting year of the reference time window for differencing
	data_type : do you want to calculate the sum and difference of extremes or original fluxes? ...
				'ext' is for extremes and will mask based on the 'bin_ar' in calculation ... 
				otherwise it will not multiply by bin_ar and the original flux difference will be calculated.

	Universal:
	----------
	start_dates :  the start_dates of every 25 year window, size = # wins

	Returns:
	--------
	sum_flux : shape (# wins, nlat,nlon), sum of fluxes per window
	diff_flux : shape (# wins, nlat,nlon), difference of sum of fluxes per window and reference window
	"""

	if data_type != 'ext': bin_ar = np.ma.ones(ano_gC.shape)
	sum_ext 	= []
	for i in range(len(start_dates)):
		ext_gC = bin_ar[idx_dates_win[i][0] : idx_dates_win [i][-1]+1,:,:] * ano_gC[idx_dates_win[i][0] : idx_dates_win [i][-1]+1,:,:]
		sum_ext . append (np.ma.sum(ext_gC, axis = 0))
	sum_ext 	= np.ma.asarray(sum_ext)
	
	#to calculate the index of the reference year starting window:
	for i,date in enumerate(start_dates):
		if date.year in [diff_ref_yr]:
			diff_yr_idx = i
	diff_ext	= []
	for i in range(len(start_dates)):
		diff	= sum_ext[i] - sum_ext[diff_yr_idx]
		diff_ext . append (diff)
	diff_ext	= np.ma.asarray(diff_ext) 
	return sum_ext , diff_ext

Results['wo_lulcc']['sum_neg_ext'], Results['wo_lulcc']['diff_neg_ext'] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar =  Results['wo_lulcc']['bin_ext_neg'], ano_gC = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], data_type = 'ext', diff_ref_yr = 1975 )
Results['wo_lulcc']['sum_pos_ext'], Results['wo_lulcc']['diff_pos_ext'] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar =  Results['wo_lulcc']['bin_ext_pos'], ano_gC = nc_data['wo_lulcc']['ano_gC' ].variables[variable][...], data_type = 'ext', diff_ref_yr = 1975 )
Results['w_lulcc' ]['sum_neg_ext'], Results['w_lulcc' ]['diff_neg_ext'] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar =  Results['w_lulcc' ]['bin_ext_neg'], ano_gC = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...], data_type = 'ext', diff_ref_yr = 1975 )
Results['w_lulcc' ]['sum_pos_ext'], Results['w_lulcc' ]['diff_pos_ext'] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar =  Results['w_lulcc' ]['bin_ext_pos'], ano_gC = nc_data['w_lulcc' ]['ano_gC' ].variables[variable][...], data_type = 'ext', diff_ref_yr = 1975 )

Results['wo_lulcc']['sum_%s'%variable], Results['wo_lulcc']['diff_%s'%variable] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar = None, ano_gC = nc_data['wo_lulcc']['var_gC' ].variables[variable][...], data_type = 'ori', diff_ref_yr = 1975 )
Results['w_lulcc' ]['sum_%s'%variable], Results['w_lulcc' ]['diff_%s'%variable] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar = None, ano_gC = nc_data['w_lulcc' ]['var_gC' ].variables[variable][...], data_type = 'ori', diff_ref_yr = 1975 )

Results ['misc']['diff_extremes'] = {}
Results ['misc']['sum_extremes'] = {}
for key in keys_wrt_wo_lulcc_pos:
	Results['misc']['sum_extremes'][key], Results['misc']['diff_extremes'][key] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar = Results['misc']['bin_ext'][key ], ano_gC = nc_data['w_lulcc']['ano_gC' ].variables[variable][...], data_type = 'ext', diff_ref_yr = 1975)
for key in keys_wrt_av_all:
	Results['misc']['sum_extremes'][key], Results['misc']['diff_extremes'][key] = Sum_and_Diff_of_Fluxes_perWin ( bin_ar = Results['misc']['bin_ext'][key ], ano_gC = nc_data['w_lulcc']['ano_gC' ].variables[variable][...], data_type = 'ext', diff_ref_yr = 1975)

# ================================================================================================================================================================================================================================= 
# ================================================================================================================================================================================================================================= 
                                       ##    #         ##     #######
                                       # #   #       ##  ##     ##
                                       ##    #       #    #     ##
                                       #     #       ##  ##     ##
									   #	 #####     ##	    ##
# ================================================================================================================================================================================================================================= 
# ================================================================================================================================================================================================================================= 
tmp_idx = np.arange(0,5401,600) #for x ticks
tmp_idx[-1]=tmp_idx[-1]-1
dates_ticks = []
for i in tmp_idx:
	a = dates_win.flatten()[i]
	dates_ticks.append(a)

# PLOTING THE THRESHOLD FOR QUALIFICATION OF EXTREME EVENTS:
# ==========================================================
if plt_threshold == 'pth':
	fig1,ax2  = plt.subplots(tight_layout = True, figsize = (9,5), dpi = 400)
	#----------------------------1
	#Hard coding the y limits
	if variable == 'gpp':
		ymin	= int(400)
		ymax	= int(851)
	#----------------------------1
	plt.title (r"The Threshold values of GPP anomalies (Percentile: %.2f)" %per)
	plt.style.use("classic")
	ax2.plot(dates_win.flatten(), abs(Results['wo_lulcc']['ts_th_neg'])/10**9 ,'r--', label = "Th$-$ without LULCC"	, alpha = .7, lw=1.5)
	ax2.plot(dates_win.flatten(), abs(Results['w_lulcc' ]['ts_th_neg'])/10**9 ,'r'  , label = "Th$-$ with LULCC" 	, alpha = .7, lw=1.5)
	ax2.set_ylabel("Negative Extremes (G%s)"%nc_data['wo_lulcc']['ano_gC' ].variables[variable].units, {'color': 'r'},fontsize =14)
	ax2.set_xlabel("Time", fontsize = 14)
	ax2.set_ylim([ymin,ymax])
	ax2.set_yticks(np.arange(ymin, ymax,50))
	ax2.set_yticklabels(-np.arange(ymin, ymax,50))
	# ax2.set_yticks(np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100),50))
	# ax2.set_yticklabels(-np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100),50))
	ax2.tick_params(axis='y', colors='red')
	ax2.set_xticks(dates_ticks)
	ax2.grid(which='both', linestyle=':', linewidth='0.3', color='gray')
	ax1=ax2.twinx()
	ax1.plot(dates_win.flatten(), Results['wo_lulcc']['ts_th_pos']/10**9,'g--', label = "Th+ without LULCC", alpha = .7, lw=1.5)
	ax1.plot(dates_win.flatten(), Results['w_lulcc' ]['ts_th_pos']/10**9,'g'  , label = "Th+ with LULCC" , alpha = .7, lw=1.5)
	#ax1.plot(dates_win.flatten(), Results['misc']['ts_av_th_all_sce']/10**9,'k'  , label = "av_all_sce" , alpha = .7)
	ax1.set_ylabel(" Positive Extremes (G%s)"%nc_data['wo_lulcc']['ano_gC' ].variables[variable].units,{'color': 'g'}, fontsize =14)
	ax1.set_ylim([ymin,ymax]) 
	ax1.set_yticks(np.arange(ymin, ymax,50))
	# ax1.set_yticks(np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100),50))
	ax1.tick_params(axis='y', colors='green')
	#ax1.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
	lines, labels = ax1.get_legend_handles_labels()
	lines2, labels2 = ax2.get_legend_handles_labels()
	ax2.legend(lines + lines2, labels + labels2, loc=0,fontsize =12)	
	fig1.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_all_scenario_%s_per_%s.pdf'%(variable,int(per)))
	fig1.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_all_scenario_%s_per_%s.png'%(variable,int(per)))
	#fig1.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_all_scenario_and_av_%s_per_%s.pdf'%(variable,int(per)))
	#fig1.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_all_scenario_and_av_%s_per_%s.png'%(variable,int(per)))
	plt.close(fig1)
	
	fig101,ax2  = plt.subplots(tight_layout = True, figsize = (9,5), dpi = 400)
	#----------------------------1
	#Hard coding the y limits
	if variable == 'gpp':
		ymin	= int(400)
		ymax	= int(701)
	#----------------------------1
	plt.title (r"The Threshold values of GPP anomalies (Percentile: %.2f)" %per)
	plt.style.use("classic")
	ax2.plot(dates_win.flatten(), abs(Results['wo_lulcc']['ts_th_pos'])/10**9 ,'k'	, alpha = .7, lw=1.5)
	ax2.set_ylabel("Negative Extremes (G%s)"%nc_data['wo_lulcc']['ano_gC' ].variables[variable].units, {'color': 'r'},fontsize =14)
	ax2.set_xlabel("Time", fontsize = 14)
	ax2.set_ylim([ymin,ymax])
	ax2.set_yticks(np.arange(ymin, ymax,50))
	ax2.set_yticklabels(-np.arange(ymin, ymax,50))
	# ax2.set_yticks(np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100),50))
	# ax2.set_yticklabels(-np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100),50))
	ax2.tick_params(axis='y', colors='red')
	ax2.set_xticks(dates_ticks)
	ax2.grid(which='both', linestyle=':', linewidth='0.3', color='gray')
	ax1=ax2.twinx()
	ax1.plot(dates_win.flatten(), Results['wo_lulcc']['ts_th_pos']/10**9,'g--',  alpha = .7, lw=1.5)
	#ax1.plot(dates_win.flatten(), Results['misc']['ts_av_th_all_sce']/10**9,'k'  , label = "av_all_sce" , alpha = .7)
	ax1.set_ylabel(" Positive Extremes (G%s)"%nc_data['wo_lulcc']['ano_gC' ].variables[variable].units,{'color': 'g'}, fontsize =14)
	ax1.set_ylim([ymin,ymax]) 
	ax1.set_yticks(np.arange(ymin, ymax,50))
	# ax1.set_yticks(np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100),50))
	ax1.tick_params(axis='y', colors='green')
	#ax1.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
	#lines, labels = ax1.get_legend_handles_labels()
	#lines2, labels2 = ax2.get_legend_handles_labels()
	#ax2.legend(lines + lines2, labels + labels2, loc=0,fontsize =12)	
	fig101.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_common_scenario_%s_per_%s.pdf'%(variable,int(per)))
	fig101.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_common_scenario_%s_per_%s.png'%(variable,int(per)))
	#fig1.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_all_scenario_and_av_%s_per_%s.pdf'%(variable,int(per)))
	#fig1.savefig(paths['out' ]['comp_lulcc']+'ts_threshold_all_scenario_and_av_%s_per_%s.png'%(variable,int(per)))
	plt.close(fig101)


# PLOTING THE GLOBAL TIMESERIES OF THE EXTREME EVENTS (INDEPENDENT ANALYSIS)
# ==========================================================================
import matplotlib.patheffects as path_effects
def plot_global_gC_timeseries(sce='wo_lulcc', ts_type = 'ts_global_gC', analysis_type = 'independent'):
	""" Plotting the global timeseries in the plot which has negative and positive y-axis
	"""
	if ts_type == 'ts_global_gC':
		Units			= nc_data[sce]['ano_gC' ].variables[variable].units # Units of gpp
		ts_div_factor	= 10**15
		y_label_text	= "Intensity of Extremes (PgC/mon)" # 'P' + Units
		slope_div_factor= 10**6
		slope_units		= 'M'+ Units
	if ts_type == 'ts_global_area':
		Units			= area.units # Units of area
		ts_div_factor	= 10**6
		y_label_text	= str(ts_div_factor) + Units
		slope_div_factor= 1
		slope_units		= Units

	if analysis_type ==  'independent':
		negext = 'neg_ext'
		posext = 'pos_ext'
		if (variable == 'gpp' and ts_type == 'ts_global_gC')  : ymin = -.9 ; ymax	= .9
		conf 		= sce
		save_path	= paths['out' ]['comp_lulcc']+'%s_%s_sce_%s_per_%d_changing_slopes'%(ts_type,variable,conf,per)
	if analysis_type == "rel_av_all":
		negext = 'neg_%s_based_av_all'%sce
		posext = 'pos_%s_based_av_all'%sce
   		conf   = 'misc'
		if (variable == 'gpp' and ts_type == 'ts_global_gC')  : ymin = -.9 ; ymax =  .9
		save_path	= paths['out' ]['comp_lulcc']+'misc/%s_%s_sce_%s_%s_per_%d_changing_slopes'%(ts_type,variable,sce,analysis_type,per)
	if analysis_type == "rel_pos_wo_lulcc":
		negext = 'neg_%s_based_pos_wo_lulcc'%sce
		posext = 'pos_%s_based_pos_wo_lulcc'%sce
		conf   = 'misc'
		if (variable == 'gpp' and ts_type == 'ts_global_gC'): ymin    = -1; ymax =  1
		save_path	= paths['out' ]['comp_lulcc']+'misc/%s_%s_sce_%s_%s_per_%d_changing_slopes'%(ts_type,variable,sce,analysis_type,per)
	
	fig2	= plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
	plt.style.use("classic")
	plt.title	("TS Global GPP Extremes for CESM1-BGC when percentile is 1.0")
	plt.ylim(ymin,ymax)
	plt.plot (dates_win.flatten()							, Results[conf][ts_type][negext]['ts'      ] / ts_div_factor, 
				'r'  , label = "Negative Extremes"     , linewidth = 0.5, alpha=0.7)
	plt.plot (dates_win.flatten()[:idx_yr_2100]				, Results[conf][ts_type][negext]['trend_21'] / ts_div_factor, 
				'k--', label = "Neg Trend 1850-2100", linewidth = 0.5, alpha=0.9)
	plt.plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]	, Results[conf][ts_type][negext]['trend_23'] / ts_div_factor, 
				'k--', label = "Neg Trend 2100-2300", linewidth = 0.5, alpha=0.9)
	plt.plot (dates_win.flatten()							, Results[conf][ts_type][posext]['ts'      ] / ts_div_factor, 
				'g'  , label = "Positive Extremes"     , linewidth = 0.5, alpha=0.7)
	plt.plot (dates_win.flatten()[:idx_yr_2100]				, Results[conf][ts_type][posext]['trend_21'] / ts_div_factor, 
				'k--', label = "Pos Trend 1850-2100", linewidth = 0.5, alpha=0.9)
	plt.plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]	, Results[conf][ts_type][posext]['trend_23'] / ts_div_factor, 
				'k--', label = "Pos Trend 2100-2300", linewidth = 0.5, alpha=0.9)
	
	plt.xlabel('Time', fontsize = 14)
	plt.ylabel(y_label_text, fontsize = 14)
	plt.xticks(dates_ticks,fontsize = 12)
	plt.yticks(fontsize = 12)
	plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
	
	text_1 = plt.text(dates_win.flatten()[600]		, ymin+0.3,
					"Slope = %d %s"%(int(Results[conf][ts_type][negext]['s21']/slope_div_factor), slope_units), size =14, color = 'r')
	text_2 = plt.text(dates_win.flatten()[600+idx_yr_2100]	, ymin+0.3,
					"Slope = %d %s"%(int(Results[conf][ts_type][negext]['s23']/slope_div_factor), slope_units), size =14, color = 'r')
	text_3 = plt.text(dates_win.flatten()[600]				, ymax-0.3,
					"Slope = %d %s"%(int(Results[conf][ts_type][posext]['s21']/slope_div_factor), slope_units), size =14, color = 'g')
	text_4 = plt.text(dates_win.flatten()[600+idx_yr_2100]	, ymax-0.3,
					"Slope = %d %s"%(int(Results[conf][ts_type][posext]['s23']/slope_div_factor), slope_units), size =14, color = 'g')
	# Setting the outline of the text
	text_1.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='black'),
                       path_effects.Normal()])
	text_2.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='black'),
                       path_effects.Normal()])
	text_3.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='black'),
                       path_effects.Normal()])
	text_4.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='black'),
                       path_effects.Normal()])
	
	fig2.savefig(save_path+'.pdf')
	fig2.savefig(save_path+'.png')
	plt.close(fig2)
if plt_global_ts == 'pgts':
	plot_global_gC_timeseries(sce = 'w_lulcc' , ts_type = 'ts_global_gC' )  
	plot_global_gC_timeseries(sce = 'wo_lulcc', ts_type = 'ts_global_gC' )
	plot_global_gC_timeseries(sce = 'wo_lulcc', ts_type = 'ts_global_gC', analysis_type = "rel_av_all")
	plot_global_gC_timeseries(sce = 'w_lulcc' , ts_type = 'ts_global_gC', analysis_type = "rel_av_all")
	plot_global_gC_timeseries(sce = 'wo_lulcc', ts_type = 'ts_global_gC', analysis_type = "rel_pos_wo_lulcc")
	plot_global_gC_timeseries(sce = 'w_lulcc' , ts_type = 'ts_global_gC', analysis_type = "rel_pos_wo_lulcc")

def plot_global_freq_area_timeseries(sce='both', ts_type = 'ts_global_freq', analysis_type = 'independent', trend_rel = 'r18'):
	""" Plotting the freq and area global timeseries, plotting two plots.
	"""
	if ts_type == 'ts_global_freq':
		Units			= ''
		ts_div_factor	= 1
		y_label_text	= 'Counts/month'
		slope_div_factor= 1
		slope_units		= '%'
	if ts_type == 'ts_global_area':
		Units			= area.units # Units of area
		ts_div_factor	= 10**6
		y_label_text	= 'Area ('+r'$10^6$ '+ r'$%s$'%Units+')'
		slope_div_factor= 1
		slope_units		= '%'

	if analysis_type ==  'independent':
		conf_0  = 'wo_lulcc'; conf_1 = 'w_lulcc'
		negext_0 = 'neg_ext'; negext_1 = 'neg_ext'
		posext_0 = 'pos_ext'; posext_1 = 'pos_ext'
		if (variable == 'gpp' and ts_type == 'ts_global_freq'): ymin = 0   ; ymax	=  650
		if (variable == 'gpp' and ts_type == 'ts_global_area'): ymin = 0   ; ymax	=  11
		conf= sce
		save_path	= paths['out' ]['comp_lulcc']+'%s_%s_sce_%s_per_%d_changing_slopes'%(ts_type,variable,sce,per)
	if analysis_type == "rel_av_all":
		conf_0   = 'wo_lulcc' ; conf_1 = 'w_lulcc'
		negext_0 = 'neg_%s_based_av_all'%conf_0 ; negext_1 = 'neg_%s_based_av_all'%conf_1
		posext_0 = 'pos_%s_based_av_all'%conf_0	; posext_1 = 'pos_%s_based_av_all'%conf_1
   		conf_0   = 'misc' ; conf_1   = 'misc'; conf='misc'
		if (variable == 'gpp' and ts_type == 'ts_global_freq'): ymin =  0  ; ymax = 750  
		if (variable == 'gpp' and ts_type == 'ts_global_area'): ymin =  0  ; ymax = 11  
		save_path	= paths['out' ]['comp_lulcc']+'misc/%s_%s_sce_%s_%s_per_%d_changing_slopes'%(ts_type,variable,sce,analysis_type,per)
	if analysis_type == "rel_pos_wo_lulcc":
		conf_0   = 'wo_lulcc' ; conf_1 = 'w_lulcc'
		negext_0 = 'neg_%s_based_pos_wo_lulcc'%conf_0 ; negext_1 = 'neg_%s_based_pos_wo_lulcc'%conf_1
		posext_0 = 'pos_%s_based_pos_wo_lulcc'%conf_0 ; posext_1 = 'pos_%s_based_pos_wo_lulcc'%conf_1
		conf_0   = 'misc' ; conf_1   = 'misc'; conf ='misc'
		if (variable == 'gpp' and ts_type == 'ts_global_freq'): ymin =  0  ; ymax = 800 
		if (variable == 'gpp' and ts_type == 'ts_global_area'): ymin =  0  ; ymax = 11
		save_path	= paths['out' ]['comp_lulcc']+'misc/%s_%s_sce_%s_%s_per_%d_changing_slopes'%(ts_type,variable,sce,analysis_type,per)
	
	fig3,ax	= plt.subplots (nrows=2,ncols=2, gridspec_kw = {'wspace':0.02, 'hspace':0.02},tight_layout = True, figsize = (9,9), dpi = 400)
	ax      = ax.ravel()
	
	ax[0].plot (dates_win.flatten()                           , Results[conf_0][ts_type][posext_0]['ts'      ] / ts_div_factor, 'g'  , label = "pos_ext"     , linewidth = 0.5, alpha=0.7)
	ax[0].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_0][ts_type][posext_0]['trend_21'] / ts_div_factor, 'k--', label = "pos_trend_21", linewidth = 0.5, alpha=0.9)
	ax[0].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_0][ts_type][posext_0]['trend_23'] / ts_div_factor, 'k--', label = "pos_trend_23", linewidth = 0.5, alpha=0.9)

	ax[1].plot (dates_win.flatten()                           , Results[conf_1][ts_type][posext_1]['ts'      ] / ts_div_factor, 'g'  , label = "pos_ext"     , linewidth = 0.5, alpha=0.7)
	ax[1].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_1][ts_type][posext_1]['trend_21'] / ts_div_factor, 'k--', label = "pos_trend_21", linewidth = 0.5, alpha=0.9)
	ax[1].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_1][ts_type][posext_1]['trend_23'] / ts_div_factor, 'k--', label = "pos_trend_23", linewidth = 0.5, alpha=0.9)

	ax[2].plot (dates_win.flatten()                           , Results[conf_0][ts_type][negext_0]['ts'      ] / ts_div_factor, 'r'  , label = "neg_ext"     , linewidth = 0.5, alpha=0.7)
	ax[2].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_0][ts_type][negext_0]['trend_21'] / ts_div_factor, 'k--', label = "neg_trend_21", linewidth = 0.5, alpha=0.9) 
	ax[2].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_0][ts_type][negext_0]['trend_23'] / ts_div_factor, 'k--', label = "neg_trend_23", linewidth = 0.5, alpha=0.9)

	ax[3].plot (dates_win.flatten()                           , Results[conf_1][ts_type][negext_1]['ts'      ] / ts_div_factor, 'r'  , label = "neg_ext"     , linewidth = 0.5, alpha=0.7)
	ax[3].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_1][ts_type][negext_1]['trend_21'] / ts_div_factor, 'k--', label = "neg_trend_21", linewidth = 0.5, alpha=0.9) 
	ax[3].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_1][ts_type][negext_1]['trend_23'] / ts_div_factor, 'k--', label = "neg_trend_23", linewidth = 0.5, alpha=0.9)
	
	for plt_idx in range(4):
		ax[plt_idx].set_xticks(dates_ticks)
		ax[plt_idx].tick_params(axis="x",direction="in")
		ax[plt_idx].set_ylim([ymin,ymax])
		ax[plt_idx].tick_params(axis="y",direction="in")
	ax[2].set_xlabel('Time (wo_lulcc)',fontsize =14)
	ax[3].set_xlabel('Time (w_lulcc)' ,fontsize =14)
	ax[0].set_ylabel('positive extremes', color='g',fontsize =12)
	ax[2].set_ylabel('negative extremes', color='r',fontsize =12)
	fig3.text(0.03, 0.5, y_label_text ,va='center', ha='center', rotation='vertical', fontsize=14)
	ax[1].set_yticklabels([])	
	ax[3].set_yticklabels([])	
	ax[0].set_xticklabels([])	
	ax[1].set_xticklabels([])
	for tick in ax[2].get_xticklabels():
		tick.set_rotation(45)
	for tick in ax[3].get_xticklabels():
		tick.set_rotation(45)

 	 #ax[0].legend(loc = 'upper center', ncol=2, bbox_to_anchor=(.44,1.15),frameon =False,fontsize=9,handletextpad = 0.1) 
	#plt.yticks(fontsize = 14)
	#plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
	ax[0].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][posext_0]['inc_21']/slope_div_factor), slope_units), size =12, color = 'g')
	ax[1].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][posext_1]['inc_21']/slope_div_factor), slope_units), size =12, color = 'g')
	ax[2].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][negext_0]['inc_21']/slope_div_factor), slope_units), size =12, color = 'r')
	ax[3].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][negext_1]['inc_21']/slope_div_factor), slope_units), size =12, color = 'r')
	ax[0].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][posext_0]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'g')
	ax[1].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][posext_1]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'g')
	ax[2].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][negext_0]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'r')
	ax[3].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][negext_1]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'r')

	if trend_rel== 'r18' : 
		fig3.savefig(save_path+'_%s.pdf'%trend_rel)
	else:
		fig3.savefig(save_path+'.pdf')
	plt.close(fig3)
if plt_global_ts == 'pgts':
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_freq' )
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_freq', analysis_type = "rel_av_all", trend_rel= 'n')
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_freq', analysis_type = "rel_av_all", trend_rel = 'r18')
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_freq', analysis_type = "rel_pos_wo_lulcc", trend_rel= 'n')
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_freq', analysis_type = "rel_pos_wo_lulcc", trend_rel = 'r18')
	
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_area', analysis_type = "rel_av_all", trend_rel= 'n')
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_area', analysis_type = "rel_av_all", trend_rel = 'r18')
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_area', analysis_type = "rel_pos_wo_lulcc", trend_rel= 'n')
	plot_global_freq_area_timeseries(sce = 'both', ts_type = 'ts_global_area', analysis_type = "rel_pos_wo_lulcc", trend_rel = 'r18')

def plot_global_freq_area_timeseries_win_step(sce='both', ts_type = 'ts_global_freq', analysis_type = 'independent', trend_rel = 'r18'):
	""" Plotting the freq and area global timeseries, plotting two plots.
	"""
	def func(a):
		b=[]
		for i in range(len(start_dates)):
			b.append( [a[i*win_len:(i+1)*win_len].sum()] *win_len)
		return np.asarray(b).flatten()

	if ts_type == 'ts_global_freq':
		Units			= ''
		ts_div_factor	= 1
		y_label_text	= 'Counts/month'
		slope_div_factor= 1
		slope_units		= '%'
	if ts_type == 'ts_global_area':
		Units			= area.units # Units of area
		ts_div_factor	= 10**6
		y_label_text	=  'Area (' r'$10^6$ '+ r'$%s$'%Units+ ')'
		slope_div_factor= 1
		slope_units		= '%'

	if analysis_type ==  'independent':
		conf_0  = 'wo_lulcc'; conf_1 = 'w_lulcc'
		negext_0 = 'neg_ext'; negext_1 = 'neg_ext'
		posext_0 = 'pos_ext'; posext_1 = 'pos_ext'
		#if (variable == 'gpp' and ts_type == 'ts_global_freq'): ymin = 0   ; ymax	=  650
		#if (variable == 'gpp' and ts_type == 'ts_global_area'): ymin = 0   ; ymax	=  11
		conf= sce
		save_path	= paths['out' ]['comp_lulcc']+'%s_%s_sce_%s_per_%d_win_step'%(ts_type,variable,sce,per)
	if analysis_type == "rel_av_all":
		conf_0   = 'wo_lulcc' ; conf_1 = 'w_lulcc'
		negext_0 = 'neg_%s_based_av_all'%conf_0 ; negext_1 = 'neg_%s_based_av_all'%conf_1
		posext_0 = 'pos_%s_based_av_all'%conf_0	; posext_1 = 'pos_%s_based_av_all'%conf_1
   		conf_0   = 'misc' ; conf_1   = 'misc'; conf='misc'
		if (variable == 'gpp' and ts_type == 'ts_global_freq'): ymin =  46500  ; ymax = 70500  
		if (variable == 'gpp' and ts_type == 'ts_global_area'): ymin =  640  ; ymax = 1050  
		save_path	= paths['out' ]['comp_lulcc']+'misc/%s_%s_sce_%s_%s_per_%d_win_step'%(ts_type,variable,sce,analysis_type,per)
	if analysis_type == "rel_pos_wo_lulcc":
		conf_0   = 'wo_lulcc' ; conf_1 = 'w_lulcc'
		negext_0 = 'neg_%s_based_pos_wo_lulcc'%conf_0 ; negext_1 = 'neg_%s_based_pos_wo_lulcc'%conf_1
		posext_0 = 'pos_%s_based_pos_wo_lulcc'%conf_0 ; posext_1 = 'pos_%s_based_pos_wo_lulcc'%conf_1
		conf_0   = 'misc' ; conf_1   = 'misc'; conf ='misc'
		if (variable == 'gpp' and ts_type == 'ts_global_freq'): ymin =  58500; ymax = 97500 
		if (variable == 'gpp' and ts_type == 'ts_global_area'): ymin =  825  ; ymax = 1350
		save_path	= paths['out' ]['comp_lulcc']+'misc/%s_%s_sce_%s_%s_per_%d_win_step'%(ts_type,variable,sce,analysis_type,per)
	
	fig3_1,ax	= plt.subplots (nrows=2,ncols=2, gridspec_kw = {'wspace':0.02, 'hspace':0.02},tight_layout = True, figsize = (9,9), dpi = 400)
	ax      = ax.ravel()
	
	ax[0].plot (dates_win.flatten()                           , func(Results[conf_0][ts_type][posext_0]['ts'      ]) / ts_div_factor, 'g'  , label = "pos_ext"     , linewidth = 0.5, alpha=0.7)
	#ax[0].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_0][ts_type][posext_0]['trend_21'] / ts_div_factor, 'k--', label = "pos_trend_21", linewidth = 0.5, alpha=0.9)
	#ax[0].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_0][ts_type][posext_0]['trend_23'] / ts_div_factor, 'k--', label = "pos_trend_23", linewidth = 0.5, alpha=0.9)

	ax[1].plot (dates_win.flatten()                           , func(Results[conf_1][ts_type][posext_1]['ts'      ]) / ts_div_factor, 'g'  , label = "pos_ext"     , linewidth = 0.5, alpha=0.7)
	#ax[1].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_1][ts_type][posext_1]['trend_21'] / ts_div_factor, 'k--', label = "pos_trend_21", linewidth = 0.5, alpha=0.9)
	#ax[1].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_1][ts_type][posext_1]['trend_23'] / ts_div_factor, 'k--', label = "pos_trend_23", linewidth = 0.5, alpha=0.9)

	ax[2].plot (dates_win.flatten()                           , func(Results[conf_0][ts_type][negext_0]['ts'      ]) / ts_div_factor, 'r'  , label = "neg_ext"     , linewidth = 0.5, alpha=0.7)
	#ax[2].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_0][ts_type][negext_0]['trend_21'] / ts_div_factor, 'k--', label = "neg_trend_21", linewidth = 0.5, alpha=0.9) 
	#ax[2].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_0][ts_type][negext_0]['trend_23'] / ts_div_factor, 'k--', label = "neg_trend_23", linewidth = 0.5, alpha=0.9)

	ax[3].plot (dates_win.flatten()                           , func(Results[conf_1][ts_type][negext_1]['ts'      ]) / ts_div_factor, 'r'  , label = "neg_ext"     , linewidth = 0.5, alpha=0.7)
	#ax[3].plot (dates_win.flatten()[:idx_yr_2100]             , Results[conf_1][ts_type][negext_1]['trend_21'] / ts_div_factor, 'k--', label = "neg_trend_21", linewidth = 0.5, alpha=0.9) 
	#ax[3].plot (dates_win.flatten()[idx_yr_2100:idx_yr_2300]  , Results[conf_1][ts_type][negext_1]['trend_23'] / ts_div_factor, 'k--', label = "neg_trend_23", linewidth = 0.5, alpha=0.9)
	
	for plt_idx in range(4):
		ax[plt_idx].set_xticks(dates_ticks)
		ax[plt_idx].tick_params(axis="x",direction="in")
		ax[plt_idx].set_ylim([ymin,ymax])
		ax[plt_idx].tick_params(axis="y",direction="in")
	ax[2].set_xlabel('Time (wo_lulcc)',fontsize =14)
	ax[3].set_xlabel('Time (w_lulcc)' ,fontsize =14)
	ax[0].set_ylabel('positive extremes', color='g',fontsize =12)
	ax[2].set_ylabel('negative extremes', color='r',fontsize =12)
	fig3_1.text(0.03, 0.5, y_label_text ,va='center', ha='center', rotation='vertical', fontsize=14)
	ax[1].set_yticklabels([])	
	ax[3].set_yticklabels([])	
	ax[0].set_xticklabels([])	
	ax[1].set_xticklabels([])
	for tick in ax[2].get_xticklabels():
		tick.set_rotation(45)
	for tick in ax[3].get_xticklabels():
		tick.set_rotation(45)

 	 #ax[0].legend(loc = 'upper center', ncol=2, bbox_to_anchor=(.44,1.15),frameon =False,fontsize=9,handletextpad = 0.1) 
	#plt.yticks(fontsize = 14)
	#plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
	#ax[0].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][posext_0]['inc_21']/slope_div_factor), slope_units), size =12, color = 'g')
	#ax[1].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][posext_1]['inc_21']/slope_div_factor), slope_units), size =12, color = 'g')
	#ax[2].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][negext_0]['inc_21']/slope_div_factor), slope_units), size =12, color = 'r')
	#ax[3].text(dates_win.flatten()[500]             , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][negext_1]['inc_21']/slope_div_factor), slope_units), size =12, color = 'r')
	#ax[0].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][posext_0]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'g')
	#ax[1].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][posext_1]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'g')
	#ax[2].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_0][ts_type][negext_0]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'r')
	#ax[3].text(dates_win.flatten()[500+idx_yr_2100] , ymax *.9,"growth= %d %s"%(int(Results[conf_1][ts_type][negext_1]['inc_23_r18' if trend_rel== 'r18' else 'inc_23']/slope_div_factor), slope_units), size =12, color = 'r')

	if trend_rel== 'r18' : 
		fig3_1.savefig(save_path+'_%s.pdf'%trend_rel)
	else:
		fig3_1.savefig(save_path+'.pdf')
	print (save_path.split('/')[-1])
	plt.close(fig3_1)
if plt_global_ts == 'pgts':
	#plot_global_freq_area_timeseries_win_step(sce = 'both', ts_type = 'ts_global_freq' )
	plot_global_freq_area_timeseries_win_step(sce = 'both', ts_type = 'ts_global_freq', analysis_type = "rel_av_all", trend_rel = 'r18')
	plot_global_freq_area_timeseries_win_step(sce = 'both', ts_type = 'ts_global_freq', analysis_type = "rel_pos_wo_lulcc", trend_rel = 'r18')
	
	plot_global_freq_area_timeseries_win_step(sce = 'both', ts_type = 'ts_global_area', analysis_type = "rel_av_all", trend_rel = 'r18')
	plot_global_freq_area_timeseries_win_step(sce = 'both', ts_type = 'ts_global_area', analysis_type = "rel_pos_wo_lulcc", trend_rel = 'r18')

# Registering a color map
import colorsys as cs
val         = 0.8
Rd          = cs.rgb_to_hsv(1,0,0)
Rd          = cs.hsv_to_rgb(Rd[0],Rd[1],val)
Gn          = cs.rgb_to_hsv(0,1,0)
Gn          = cs.hsv_to_rgb(Gn[0],Gn[0],val)
RdGn        = {'red'  : ((0.0,  0.0,    Rd[0]),
		                 (0.5,  1.0,    1.0  ),
						 (1.0,  Gn[0],  0.0  )),
			   'green': ((0.0,  0.0,    Rd[1]),
					     (0.5,  1.0,    1.0  ),
						 (1.0,  Gn[1],  0.0  )),
			   'blue' : ((0.0,  0.0,    Rd[2]),
					     (0.5,  1.0,    1.0  ),
						 (1.0,  Gn[2],  0.0  ))}
plt.register_cmap(name  = 'RdGn',data = RdGn)

def Plot_Diff_Plot(diff_array, datatype='neg_ext', sce = 'wo_lulcc', analysis_type='independent'):
	""" Plotting the change or difference maps for the given data and passing the following filters for plotting
	
	Parameters:
	-----------
	diff_array: the difference array of shape  (nwin,nlat,nlon)
	datatype : 'ext' or 'var_ori'
	analysis type: 'independent or 'rel_av_all' or 'rel_pos_wo_wulcc'

	"""
	text = []
	for i in start_dates:
		a = i.year
		b = a+24
		c = str(a)+'-'+str(b)[2:]
		text.append(c)
	#to calculate the index of the reference year starting window:
	diff_ref_yr = 1975
	for i,date in enumerate(start_dates):
		if date.year in [diff_ref_yr]:
			diff_yr_idx = i

	if (variable=='gpp' and datatype=='var_ori' and analysis_type=='independent'): 
		ymin = -500
		ymax = 500
		Units			= nc_data[sce]['var_gC' ].variables[variable].units # Units of gpp
		ts_div_factor	= 10**12
		y_label_text	= 'T' + Units
		save_path   = paths['out' ]['comp_lulcc']+'diff_plots/change_%s_sce_%s_'%(variable,sce)

	#if (variable=='gpp' and datatype=='neg_ext' and analysis_type=='independent'): 
	if (variable=='gpp' and datatype=='neg_ext'):
		ymin = -14
		ymax = 14
		Units			= nc_data[sce]['ano_gC' ].variables[variable].units # Units of gpp
		ts_div_factor	= 10**12
		y_label_text	= 'T' + Units

		if analysis_type == 'independent'		: save_path   = paths['out' ]['comp_lulcc']+'diff_plots/change_%s_negext_sce_%s_'%(variable,sce)
		if analysis_type == "rel_pos_wo_lulcc"	: save_path   = paths['out' ]['comp_lulcc']+'misc/diff_plots/change_%s_negext_sce_%s_rel_pos_wo_lulcc_'%(variable,sce)	
		if analysis_type == "rel_av_all"		: save_path   = paths['out' ]['comp_lulcc']+'misc/diff_plots/change_%s_negext_sce_%s_rel_av_all_'%(variable,sce)	
	if (variable=='gpp' and datatype=='pos_ext'): 
		ymin = -14
		ymax = 14
		Units			= nc_data[sce]['ano_gC' ].variables[variable].units # Units of gpp
		ts_div_factor	= 10**12
		y_label_text	= 'T' + Units
		if analysis_type == 'independent'		: save_path   = paths['out' ]['comp_lulcc']+'diff_plots/change_%s_posext_sce_%s_'%(variable,sce)
		if analysis_type == "rel_pos_wo_lulcc"	: save_path   = paths['out' ]['comp_lulcc']+'misc/diff_plots/change_%s_posext_sce_%s_rel_pos_wo_lulcc_'%(variable,sce)
		if analysis_type == "rel_av_all"		: save_path   = paths['out' ]['comp_lulcc']+'misc/diff_plots/change_%s_posext_sce_%s_rel_av_all_'%(variable,sce)	
	for i, data in enumerate (diff_array):
		print "Plotting %s %s %s for the win %d"%(datatype,sce,analysis_type,i)
		fig4,ax	= plt.subplots(figsize = (7,2.8),tight_layout=True,dpi=500)
		bmap    = Basemap(  projection  =   'eck4',
							lon_0       =   0.,
							resolution  =   'c')
		LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij')
		ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_invalid(data/ts_div_factor),latlon=True,cmap= 'RdGn',vmax= ymax, vmin= ymin)
		cbar    = plt.colorbar(ax)
		cbar    .ax.set_ylabel(y_label_text)
		bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
		bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
		bmap    .drawcoastlines(linewidth = .2)
		plt.title(" %s minus %s" %(text[i], text[diff_yr_idx]))
		fig4.savefig(save_path + '%s.pdf'%format(i,'02'))
		plt.close(fig4)


# Calling the functions:
Plot_Diff_Plot(diff_array = Results['wo_lulcc']['diff_%s'%variable], datatype='var_ori',  sce = 'wo_lulcc', analysis_type='independent')
Plot_Diff_Plot(diff_array = Results['w_lulcc' ]['diff_%s'%variable], datatype='var_ori',  sce = 'w_lulcc' , analysis_type='independent')

Plot_Diff_Plot(diff_array = Results['wo_lulcc']['diff_neg_ext']    , datatype='neg_ext',  sce = 'wo_lulcc', analysis_type='independent')
Plot_Diff_Plot(diff_array = Results['w_lulcc' ]['diff_neg_ext']    , datatype='neg_ext',  sce = 'w_lulcc' , analysis_type='independent')

Plot_Diff_Plot(diff_array = Results['wo_lulcc']['diff_pos_ext']    , datatype='pos_ext',  sce = 'wo_lulcc', analysis_type='independent')
Plot_Diff_Plot(diff_array = Results['w_lulcc' ]['diff_pos_ext']    , datatype='pos_ext',  sce = 'w_lulcc' , analysis_type='independent')

Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['neg_wo_lulcc_based_pos_wo_lulcc']    , datatype='neg_ext',  sce = 'wo_lulcc', analysis_type='rel_pos_wo_lulcc')
Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['neg_w_lulcc_based_pos_wo_lulcc' ]    , datatype='neg_ext',  sce = 'w_lulcc' , analysis_type='rel_pos_wo_lulcc')

Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['pos_wo_lulcc_based_pos_wo_lulcc']    , datatype='pos_ext',  sce = 'wo_lulcc', analysis_type='rel_pos_wo_lulcc')
Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['pos_w_lulcc_based_pos_wo_lulcc' ]    , datatype='pos_ext',  sce = 'w_lulcc' , analysis_type='rel_pos_wo_lulcc')

Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['neg_wo_lulcc_based_av_all']    , datatype='neg_ext',  sce = 'wo_lulcc', analysis_type='rel_av_all')
Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['neg_w_lulcc_based_av_all' ]    , datatype='neg_ext',  sce = 'w_lulcc' , analysis_type='rel_av_all')

Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['pos_wo_lulcc_based_av_all']    , datatype='pos_ext',  sce = 'wo_lulcc', analysis_type='rel_av_all')
Plot_Diff_Plot(diff_array = Results['misc']['diff_extremes']['pos_w_lulcc_based_av_all' ]    , datatype='pos_ext',  sce = 'w_lulcc' , analysis_type='rel_av_all')


def Plot_Spatial_Freq_Win(binary_array, win_start_year = 2000, filename_out = "spatial_freq/spatial_freq_win_2000-24_wo_lulcc"):
	#to calculate the index of the reference year starting window:
	for i,date in enumerate(start_dates):
		if date.year in [win_start_year]:
			start_yr_idx = i
	save_path   = paths['out' ]['comp_lulcc'] + filename_out 

	data = binary_array[start_yr_idx]

	fig5,ax	= plt.subplots(figsize = (7,2.8),tight_layout=True,dpi=500)
	bmap    = Basemap(  projection  =   'eck4',
						lon_0       =   0.,
						resolution  =   'c')
	LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij')
	ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_invalid(data),latlon=True,cmap= cm.afmhot_r)#,vmax= ymax, vmin= ymin)
	cbar    = plt.colorbar(ax)
	cbar    .ax.set_ylabel('Counts')
	bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
	bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
	bmap    .drawcoastlines(linewidth = .2)
	plt.title("2000-24 and total : %d"%int(data.sum()))
	fig5.savefig(save_path + '.pdf')
	plt.close(fig5)

Plot_Spatial_Freq_Win ( binary_array = Results['wo_lulcc']['sp_freq_neg'], filename_out = "spatial_freq/spatial_freq_win_2000-24_wo_lulcc_neg_ext")
Plot_Spatial_Freq_Win ( binary_array = Results['wo_lulcc']['sp_freq_pos'], filename_out = "spatial_freq/spatial_freq_win_2000-24_wo_lulcc_pos_ext")
Plot_Spatial_Freq_Win ( binary_array = Results['w_lulcc' ]['sp_freq_neg'], filename_out = "spatial_freq/spatial_freq_win_2000-24_w_lulcc_neg_ext" )
Plot_Spatial_Freq_Win ( binary_array = Results['w_lulcc' ]['sp_freq_pos'], filename_out = "spatial_freq/spatial_freq_win_2000-24_w_lulcc_pos_ext" )
for key in Results['misc']['bin_ext'].keys():
	Plot_Spatial_Freq_Win ( binary_array = Results['misc']['sp_freq'][key], filename_out = "misc/spatial_freq/spatial_freq_win_2000-24_%s"%key )


def Plot_Spatial_Freq_TCE_Win(sp_count_TCE, win_start_year = 2000, filename_out = "spatial_freq/spatial_freq_TCE_win_2000-24_wo_lulcc"):
	#to calculate the index of the reference year starting window:
	for i,date in enumerate(start_dates):
		if date.year in [win_start_year]:
			start_yr_idx = i
	save_path   = paths['out' ]['comp_lulcc'] + filename_out 

	fig6,ax	= plt.subplots(figsize = (7,2.8),tight_layout=True,dpi=500)
	bmap    = Basemap(  projection  =   'eck4',
						lon_0       =   0.,
						resolution  =   'c')
	LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij')
	ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_invalid(sp_count_TCE),latlon=True,cmap= cm.afmhot_r)#,vmax= ymax, vmin= ymin)
	cbar    = plt.colorbar(ax)
	cbar    .ax.set_ylabel('Counts')
	bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
	bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
	bmap    .drawcoastlines(linewidth = .2)
	plt.title("2000-24 and total : %d"%int(sp_count_TCE.sum()))
	fig6.savefig(save_path + '.pdf')
	plt.close(fig6)

Plot_Spatial_Freq_TCE_Win ( sp_count_TCE = Results['wo_lulcc']['sp_freq_neg_TCE'], filename_out = "spatial_freq/spatial_freq_TCE_win_2000-24_wo_lulcc_neg_ext")
Plot_Spatial_Freq_TCE_Win ( sp_count_TCE = Results['wo_lulcc']['sp_freq_pos_TCE'], filename_out = "spatial_freq/spatial_freq_TCE_win_2000-24_wo_lulcc_pos_ext")
Plot_Spatial_Freq_TCE_Win ( sp_count_TCE = Results['w_lulcc' ]['sp_freq_neg_TCE'], filename_out = "spatial_freq/spatial_freq_TCE_win_2000-24_w_lulcc_neg_ext" )
Plot_Spatial_Freq_TCE_Win ( sp_count_TCE = Results['w_lulcc' ]['sp_freq_pos_TCE'], filename_out = "spatial_freq/spatial_freq_TCE_win_2000-24_w_lulcc_pos_ext" )
for key in Results['misc']['bin_ext'].keys():
	Plot_Spatial_Freq_TCE_Win ( sp_count_TCE = Results['misc']['sp_freq_TCE'][key], filename_out = "misc/spatial_freq/spatial_freq_TCE_win_2000-24_%s"%key )


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	
### def Calc_anomalies(pft_idx_in):
### 	""" This function is designed to calculate the both positive and negative anomalies thresholds, global time series, count of frequency of events relative to the (1850 - 1990) average
### 		The threshold is calculated using consecutive 'window'years 
### 		if pft_idx_in is in [0,16] ano_gC * ano_pft : only anomalies pertaining to the pft type would be considered 
### 		else if pft_idx_in == 99 then, all_pft case will 
### 		"""
### 	idx_dates_win= []	#indicies of time in 30yr windows
### 	dates_win	= []	#sel dates from time variables in win_len windows
### 	thresholds_1= []	#list of the anomalies threshold for consecutive win_len yr corresponding to 'per'
### 	thresholds_2= []  	#list of the anomalies threshold for consecutive win_len yr corresponding to '100-per'
### 	bin_ano_neg	= np.zeros((ano_gC.shape)) #3d array to capture the True binaray anomalies w.r.t. gpp loss events
### 	bin_ano_pos	= np.zeros((ano_gC.shape)) #3d array to capture the True binaray anomalies w.r.t. gpp gain events
### 
### 	if pft_idx_in in pft_idx[1:-2]:
### 		pft_n 		= pft_rel_pct[pft_idx_in].copy()
### 		ano_gC_pft	= pft_n*ano_gC
### 	elif pft_idx_in == 99:
### 	 	ano_gC_pft	= ano_gC.copy()
### 	else:
### 		ano_gC_pft  = ano_gC.copy()
### 
### 	for i in range(len(start_dates)):
### 		idx_loc, dates_loc 	= index_and_dates_slicing(dates_ar,start_dates[i],end_dates[i]) # see functions.py
### 		idx_dates_win		. append	(idx_loc)
### 		dates_win			. append	(dates_loc)
### 		ano_loc 			= ano_gC_pft[idx_loc[0]:idx_loc[-1]+1,:,:]
### 	
### 		threshold_loc_1	= np.percentile(ano_loc[ano_loc.mask == False],per) # calculation of threshold for the local anomalies	
### 		thresholds_1	. append(threshold_loc_1)
### 		if plt_global_ts == 'pgts' or plt_freq == 'pf' :
### 			if per <=50:
### 				bin_ano_neg[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc < threshold_loc_1
### 			else:
### 				bin_ano_pos[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc > threshold_loc_1
### 			
### 		if per_both		== 'b': #in case we want to analyse both neg and positive percentiles
### 			threshold_loc_2	= np.percentile(ano_loc[ano_loc.mask == False],(100-per))
### 			thresholds_2	. append(threshold_loc_2)
### 			if plt_global_ts == 'pgts' or plt_freq == 'pf' :
### 				if 100-per <=50:
### 					bin_ano_neg[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc	< threshold_loc_2
### 				else:
### 					bin_ano_pos[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc > threshold_loc_2
### 
### 	dates_win           =   np.array(dates_win)     
### 	idx_dates_win       =   np.array(idx_dates_win)
### 	
### 	#Making the spatial map of frequncy of events:
### 	if plt_freq == 'pf' and per_both == 'b':
### 		freq_pos	= np.array([np.ma.equal(np.sum(bin_ano_pos[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:],axis=0),0) for i in range(idx_dates_win.shape[0])])
### 		freq_neg	= np.array([np.ma.equal(np.sum(bin_ano_neg[idx_dates_win[i][0]:idx_dates_win[i][-1]+1,:,:],axis=0),0) for i in range(idx_dates_win.shape[0])])
### 
### 	#calculating the global time series if the 4th arg/plt_global_ts = 'pgts'
### 	if plt_global_ts == 'pgts' and per_both != 'b':
### 		if per <50:
### 			gpp_loss = bin_ano_neg * ano_gC_pft
### 			gpp_gain = None
### 		else:
### 			gpp_gain = bin_ano_pos * ano_gC_pft
### 			gpp_loss = None
### 
### 	if plt_global_ts == 'pgts' and per_both == 'b':
### 		gpp_loss = bin_ano_neg * ano_gC_pft
### 		gpp_gain = bin_ano_pos * ano_gC_pft
### 		area_loss= bin_ano_neg * area[...]		#area matrix corresponding to the negative anomalies
### 		area_gain= bin_ano_pos * area[...]		#area matrix corresponding to the positive anomalies	
### 
### 	elif plt_global_ts == 'pgts' and per_both != 'b':
### 		if per <= 50:
### 			gpp_loss = bin_ano_neg * ano_gC_pft
### 		else:
### 			gpp_gain = bin_ano_pos * ano_gC_pft
### 
### 	# calculating the data for plotting the threshold vs dates if the 3rd argument is 'pth'
### 	if plt_threshold    == 'pth':
### 		threshold_ts1   = np.array([np.array([thresholds_1[i]]*win_len)for i in range(len(thresholds_1))])
### 		threshold_ts2   =  None
### 		if per_both		== 'b':
### 			threshold_ts2 = np.array([np.array([thresholds_2[i]]*win_len)for i in range(len(thresholds_2))])
### 	#calculation for the global time series if the 4th arg/plt_global_ts = 'pgts'
### 	global_gpp_loss	= None
### 	global_gpp_gain = None
### 	slope_loss		= None
### 	slope_gain		= None
### 	trend_loss		= None
### 	trend_gain		= None
### 	if plt_global_ts == 'pgts' and per_both =='b':
### 		global_gpp_loss	= []
### 		global_gpp_gain	= []
### 		freq_loss		= []
### 		freq_gain		= []
### 		for i in range(dates_win.flatten().size):
### 			global_gpp_loss .append(np.sum(gpp_loss[i]))	#global gpp loss TS, extremes only, adding over all the spatial points for every timestamp
### 			global_gpp_gain .append(np.sum(gpp_gain[i]))	#global gpp gain TS, extremes only
### 			freq_loss    	.append(np.sum(bin_ano_neg[i]))	#count of gpp loss events
### 			freq_gain		.append(np.sum(bin_ano_pos[i]))	##count of gpp gain events
### 		""" calculation of plotting the real magnitudes of global gpp loss and trend"""
### 	 	slope_loss, intercept_loss,_,_,_  = stats.linregress(time[...][:idx_yr_2299],global_gpp_loss)
### 		trend_loss       			=   slope_loss*time[...][:idx_yr_2299]+intercept_loss
### 		increase_loss    			=   (trend_loss[-1]-trend_loss[0])*100/trend_loss[0]
### 	
### 		""" calculation of plotting the real magnitudes of global gpp gain and trend"""
### 		slope_gain, intercept_gain,_,_,_  = stats.linregress(time[...][:idx_yr_2299],global_gpp_gain)
### 		trend_gain	      			=   slope_gain*time[...][:idx_yr_2299]+intercept_gain
### 		increase_gain    			=   (trend_gain[-1]-trend_gain[0])*100/trend_gain[0]
### 	
### 		""" calculation of the magnitudes of global gpp loss and trend from 1850-2100 """
### 	 	slope_loss_21, intercept_loss_21,_,_,_  = stats.linregress(time[...][:idx_yr_2100],global_gpp_loss[:idx_yr_2100])
### 		trend_loss_21       			=   slope_loss_21*time[...][:idx_yr_2100]+intercept_loss_21
### 		increase_loss_21    			=   (trend_loss_21[-1]-trend_loss_21[0])*100/trend_loss_21[0]
### 	
### 		""" calculation of the magnitudes of global gpp gain and trend from 1850-2100"""
### 		slope_gain_21, intercept_gain_21,_,_,_  = stats.linregress(time[...][:idx_yr_2100],global_gpp_gain[:idx_yr_2100])
### 		trend_gain_21	      			=   slope_gain_21*time[...][:idx_yr_2100]+intercept_gain_21
### 		increase_gain_21    			=   (trend_gain_21[-1]-trend_gain_21[0])*100/trend_gain_21[0]
### 
### 		""" calculation of the magnitudes of global gpp loss and trend from 2101-2300 """
### 	 	slope_loss_23, intercept_loss_23,_,_,_  = stats.linregress(time[...][idx_yr_2100:idx_yr_2299],global_gpp_loss[idx_yr_2100:idx_yr_2299])
### 		trend_loss_23       			=   slope_loss_23*time[...][idx_yr_2100:idx_yr_2299]+intercept_loss_23
### 		increase_loss_23    			=   (trend_loss_23[-1]-trend_loss_23[0])*100/trend_loss_23[0]
### 	
### 		""" calculation of the magnitudes of global gpp gain and trend from 2101-2300"""
### 		slope_gain_23, intercept_gain_23,_,_,_  = stats.linregress(time[...][idx_yr_2100:idx_yr_2299],global_gpp_gain[idx_yr_2100:idx_yr_2299])
### 		trend_gain_23	      			=   slope_gain_23*time[...][idx_yr_2100:idx_yr_2299]+intercept_gain_23
### 		increase_gain_23    			=   (trend_gain_23[-1]-trend_gain_23[0])*100/trend_gain_23[0]
### 
### 	elif plt_global_ts == 'pgts'and per_both !='b':
### 		global_gpp_loss = []
### 		freq_loss       = []
### 		global_gpp_gain = []
### 		freq_gain       = []
### 		if per <= 50:
### 			for i in range(dates_win.flatten().size):
### 				global_gpp_loss .append(np.sum(gpp_loss[i]))    #global gpp loss TS, extremes only
### 				freq_loss       .append(np.sum(bin_ano_neg[i])) #count of gpp loss events
### 			""" calculation of plotting the real magnitudes of global gpp loss and trend"""
### 			slope_loss, intercept_loss,_,_,_  = stats.linregress(time[...][:dates_win.flatten().size],global_gpp_loss)
### 			trend_loss       			=   slope_loss*time[...][:dates_win.flatten().size]+intercept_loss
### 			increase_loss    			=   (trend_loss[-1]-trend_loss[0])*100/trend_loss[0]
### 			trend_gain, slope_gain		=	None, None
### 		if per > 50:
### 		  	for i in range(dates_win.flatten().size):
### 				global_gpp_gain .append(np.sum(gpp_gain[i]))    #global gpp gain TS, extremes only 
### 				freq_gain       .append(np.sum(bin_ano_pos[i])) ##count of gpp gain events
### 				""" calculation of plotting the real magnitudes of global gpp gain and trend"""
### 			slope_gain, intercept_gain,_,_,_  = stats.linregress(time[...][:dates_win.flatten().size],global_gpp_gain)
### 			trend_gain	      			=   slope_gain*time[...][:dates_win.flatten().size]+intercept_gain
### 			increase_gain    			=   (trend_gain[-1]-trend_gain[0])*100/trend_gain[0]
### 			trend_loss, slope_loss		=	None, None
### 
### 	#Calculation of the area associated with the extremes
### 	if plt_global_ts == 'pgts' and per_both =='b':
### 		global_area_loss = []
### 		global_area_gain = []
### 		for i in range(dates_win.flatten().size): 
### 			global_area_loss .append(np.sum(area_loss[i]))    #global area loss TS, extremes only, adding over all the spatial points for every timestamp
### 			global_area_gain .append(np.sum(area_gain[i]))    #global area gain TS, extremes only
### 
### 	#Calculate the extremes if the thresholds of uptill year 1999 is considered as the bench mark
### 	freq_loss_r1999 = None
### 	freq_gain_r1999 = None
### 	trend_gain_r1999= None
### 	trend_loss_r1999= None
### 	slope_loss_r1999= None
### 	slope_gain_r1999= None
### 
### 	if plt_count_r1999 == 'pcr1999':
### 		threshold_r1999_1	=	np.mean(thresholds_1[0:5])
### 		if per	<= 50:
### 			bin_ano_neg_r1999	=	ano_gC_pft < threshold_r1999_1
### 		else:
### 			bin_ano_pos_r1999    =   ano_gC_pft > threshold_r1999_1
### 		
### 		if per_both == 'b':
### 			threshold_r1999_2	=	np.mean(thresholds_2[0:5])
### 			if 100-per	<= 50:
### 				bin_ano_neg_r1999	=	ano_gC_pft < threshold_r1999_2
### 			else:
### 				bin_ano_pos_r1999    =   ano_gC_pft > threshold_r1999_2
### 
### 		#calculating the count of extreme events timeseries relative to the 1850-1999 average
### 		if per_both == 'b':
### 			freq_loss_r1999	=	[]
### 			freq_gain_r1999 =   []
### 			for i in range(time.size):
### 				freq_loss_r1999	.append(float(np.sum(bin_ano_neg_r1999[i])))
### 				freq_gain_r1999 .append(float(np.sum(bin_ano_pos_r1999[i])))
### 			""" calculation for plotting the count of gpp loss and trend"""
### 			slope_loss_r1999, intercept_loss_r1999,_,_,_  = stats.linregress(time[...],freq_loss_r1999) 		#time[...] is in days so slope is counts/days
### 			trend_loss_r1999       			=   slope_loss_r1999*time[...]+intercept_loss_r1999
### 			increase_loss_r1999    			=   (trend_loss_r1999[-1]-trend_loss_r1999[0])*100/trend_loss_r1999[0]
### 			""" calculation for plotting the count of gpp gain and trend"""
### 			slope_gain_r1999, intercept_gain_r1999,_,_,_  = stats.linregress(time[...],freq_gain_r1999)
### 			trend_gain_r1999       			=   slope_gain_r1999*time[...]+intercept_gain_r1999
### 			increase_gain_r1999    			=   (trend_gain_r1999[-1]-trend_gain_r1999[0])*100/trend_gain_r1999[0]
### 
### 			""" calculation for plotting the count of gpp loss and trend for 1850-2100"""
### 			slope_loss_r1999_21, intercept_loss_r1999_21,_,_,_  = stats.linregress(time[...][:idx_yr_2100],freq_loss_r1999[:idx_yr_2100]) 		#time[...] is in days so slope is counts/days
### 			trend_loss_r1999_21       			=   slope_loss_r1999_21*time[...][:idx_yr_2100]+intercept_loss_r1999_21
### 			increase_loss_r1999_21    			=   (trend_loss_r1999_21[-1]-trend_loss_r1999_21[0])*100/trend_loss_r1999_21[0]
### 			""" calculation for plotting the count of gpp gain and trend for 1850-2100"""
### 			slope_gain_r1999_21, intercept_gain_r1999_21,_,_,_  = stats.linregress(time[...][:idx_yr_2100],freq_gain_r1999[:idx_yr_2100])
### 			trend_gain_r1999_21       			=   slope_gain_r1999_21*time[...][:idx_yr_2100]+intercept_gain_r1999_21
### 			increase_gain_r1999_21    			=   (trend_gain_r1999_21[-1]-trend_gain_r1999_21[0])*100/trend_gain_r1999_21[0]
### 
### 			""" calculation for plotting the count of gpp loss and trend for 2100-2300"""
### 			slope_loss_r1999_23, intercept_loss_r1999_23,_,_,_  = stats.linregress(time[...][idx_yr_2100:idx_yr_2299],freq_loss_r1999[idx_yr_2100:idx_yr_2299]) 		#time[...] is in days so slope is counts/days
### 			trend_loss_r1999_23       			=   slope_loss_r1999_23*time[...][idx_yr_2100:idx_yr_2299]+intercept_loss_r1999_23
### 			increase_loss_r1999_23    			=   (trend_loss_r1999_23[-1]-trend_loss_r1999_23[0])*100/trend_loss_r1999_23[0]
### 			""" calculation for plotting the count of gpp gain and trend for 2100-2300"""
### 			slope_gain_r1999_23, intercept_gain_r1999_23,_,_,_  = stats.linregress(time[...][idx_yr_2100:idx_yr_2299],freq_gain_r1999[idx_yr_2100:idx_yr_2299])
### 			trend_gain_r1999_23       			=   slope_gain_r1999_23*time[...][idx_yr_2100:idx_yr_2299]+intercept_gain_r1999_23
### 			increase_gain_r1999_23    			=   (trend_gain_r1999_23[-1]-trend_gain_r1999_23[0])*100/trend_gain_r1999_23[0]
### 
### 
### 
### 		else:
### 			if per <50:	
### 				freq_loss_r1999 = [np.sum(bin_ano_neg_r1999[i]) for i in range (time.size)]
### 			else:
### 				freq_gain_r1999	= [np.sum(bin_ano_pos_r1999[i]) for i in range (time.size)]
### 
### 	#calculating the per pixel timeseries
### 	pp_ts	=	None
### 	if -90 <= lat_in <=90 	and 	0 <=lon_in <=360:
### 		lat_idx =   geo_idx(float(lat_in),lat[...])
### 		lon_idx =   geo_idx(float(lon_in),lon[...])
### 		pp_ts  =   ano_gC_pft[:-(time.size - dates_win.size)][:,lat_idx,lon_idx] #per pixel timeseries
### 	
### 	slope_loss_tsw = None
### 	slope_gain_tsw = None
### 	nrows_tsw	   = None
### 	ncols_tsw	   = None
### 
### 	#plotting the timeseries for window years
### 	global_loss_tsw = None; global_gain_tsw = None; trend_loss_tsw  = None; trend_gain_tsw  = None
### 	if p_win	== 	'pwin' and per_both== 'b':
### 		global_loss_tsw	= []
### 		global_gain_tsw = []
### 		trend_loss_tsw	= []
### 		trend_gain_tsw	= []
### 		slope_loss_tsw	= []
### 		slope_gain_tsw	= []
### 		subs			= dates_win.shape[0]
### 		nrows_tsw		= int(np.sqrt(subs))
### 		ncols_tsw		= int(subs/nrows_tsw) if subs%nrows_tsw == 0 else (int(subs/nrows_tsw) +1)
### 		for i in range(subs):
### 			global_loss_loc	= global_gpp_loss[idx_dates_win[i][0]:idx_dates_win[i][-1]+1]
### 			slope_loss_loc,intercept_loss_loc,_,_,_ = stats.linregress(time[...][idx_dates_win[i][0]:idx_dates_win[i][-1]+1],global_loss_loc) 
### 			trend_loss_loc	= slope_loss_loc*time[...][idx_dates_win[i][0]:idx_dates_win[i][-1]+1] + intercept_loss_loc
### 			increase_loss	= (trend_loss_loc[-1]-trend_loss_loc[0])*100/trend_loss_loc[0]	
### 	
### 			global_gain_loc	= global_gpp_gain[idx_dates_win[i][0]:idx_dates_win[i][-1]+1]
### 			slope_gain_loc,intercept_gain_loc,_,_,_ = stats.linregress(time[...][idx_dates_win[i][0]:idx_dates_win[i][-1]+1],global_gain_loc) 
### 			trend_gain_loc	= slope_gain_loc*time[...][idx_dates_win[i][0]:idx_dates_win[i][-1]+1] + intercept_gain_loc
### 			increase_gain	= (trend_gain_loc[-1]-trend_gain_loc[0])*100/trend_gain_loc[0]	
### 
### 			global_loss_tsw	.append(global_loss_loc)
### 			global_gain_tsw .append(global_gain_loc)
### 			trend_loss_tsw 	.append(trend_loss_loc)
### 			trend_gain_tsw 	.append(trend_gain_loc)
### 			slope_loss_tsw	.append(slope_loss_loc)
### 			slope_gain_tsw	.append(slope_gain_loc)
### 
### 
### 	return dates_win,idx_dates_win,threshold_ts1,threshold_ts2,np.array(global_gpp_loss),np.array(global_gpp_gain),slope_loss,slope_gain,np.array(trend_loss),np.array(trend_gain),slope_loss_21,slope_gain_21,np.array(trend_loss_21),np.array(trend_gain_21),slope_loss_23,slope_gain_23,np.array(trend_loss_23),np.array(trend_gain_23),freq_loss_r1999,freq_gain_r1999,slope_loss_r1999,slope_gain_r1999,trend_loss_r1999,trend_gain_r1999,slope_loss_r1999_21,slope_gain_r1999_21,trend_loss_r1999_21,trend_gain_r1999_21,slope_loss_r1999_23,slope_gain_r1999_23,trend_loss_r1999_23,trend_gain_r1999_23,pp_ts,global_loss_tsw,global_gain_tsw,trend_loss_tsw,trend_gain_tsw, slope_loss_tsw, slope_gain_tsw, nrows_tsw,ncols_tsw
### 
### Results_calc_pft_idx=	Calc_anomalies(pft_idx_in)
### Results_calc_pft_all=	Calc_anomalies(99)
### Results_var_name	=	['dates_win','idx_dates_win','threshold_ts1','threshold_ts2','global_gpp_loss','global_gpp_gain','slope_loss','slope_gain','trend_loss','trend_gain','slope_loss_21','slope_gain_21','trend_loss_21','trend_gain_21','slope_loss_23','slope_gain_23','trend_loss_23','trend_gain_23','freq_loss_r1999','freq_gain_r1999','slope_loss_r1999','slope_gain_r1999','trend_loss_r1999','trend_gain_r1999','slope_loss_r1999_21','slope_gain_r1999_21','trend_loss_r1999_21','trend_gain_r1999_21','slope_loss_r1999_23','slope_gain_r1999_23','trend_loss_r1999_23','trend_gain_r1999_23', 'pp_ts','global_loss_tsw','global_gain_tsw','trend_loss_tsw','trend_gain_tsw', 'slope_loss_tsw' , 'slope_gain_tsw', 'nrows_tsw','ncols_tsw']
### #calling all results to capture:
### Results						=	{}
### Results	['all_pft']			=	{}  # the analysis on total anomalies (i.e. all pfts are contributing to it
### Results	['sel_pft']			=	{}	# the analysis on a selected pft type anomalies, nothing to to with the all_pft_anomalies
### Results	['comp_pft']		=	{}	# the analysis when the threshold of the all_pft anomalies is used and sel_pft is analysed accordingly
### Results	['cumsum_neg_pft']	=	{} 	# the cummulative sum of neg extremes global timeseries when keys(pft_type) are in decending order of mean
### Results	['cumsum_pos_pft']  =   {}	# simmilar to Results ['cumsum_neg_pft'], but for positive anomalies
### Results ['compwin_cumsum_neg']=	{}	# the cummulative sum of neg extremes global timeseries for all time windows  when keys(pft_type) are in decending order of mean 
### Results ['compwin_cumsum_pos']=	{}	# the cummulative sum of pos extremes global timeseries for all time windows  when keys(pft_type) are in decending order of mean
### 
### for i in range(len(Results_var_name)):
### 	Results ['all_pft'][Results_var_name[i]] = Results_calc_pft_all[i]
### 	Results ['sel_pft'][Results_var_name[i]] = Results_calc_pft_idx[i]
### def Impact(ar):
### 	ar=np.array(ar)
### 	impact	= (ar[-1]-ar[0])*100/ar[0]
### 	return impact
### 
### #comparision using the all_pft threshold
### def Calc_pft_comp(pft_idx_in):
### 	""" if the thresholds of the 'all_pft' analysis is used to define the anomalies of the 'sel_pft',
### 		how many anomalies will be calculated """
### 	ano_gC_pft	=	None
### 	if pft_comp == 'comp':
### 		bin_ano_neg = np.zeros((ano_gC.shape)) #3d array to capture the True binaray anomalies w.r.t. gpp loss events
### 		bin_ano_pos = np.zeros((ano_gC.shape)) #3d array to capture the True binaray anomalies w.r.t. gpp gain events
### 		if pft_idx_in in pft_idx[1:-2]:	
### 			pft_n       = pft_rel_pct[pft_idx_in].copy()
### 			ano_gC_pft  = pft_n*ano_gC
### 		for i in range(Results['all_pft']['idx_dates_win'].shape[0]):
### 			idx_loc			= Results['all_pft']['idx_dates_win'][i]
### 			ano_loc     	= ano_gC_pft[idx_loc[0]:idx_loc[-1]+1,:,:]
### 			threshold_loc_1	= np.mean(Results['all_pft']['threshold_ts1'][i])
### 			if plt_global_ts == 'pgts' or plt_freq == 'pf' :
### 				if per <=50:
### 					bin_ano_neg[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc < threshold_loc_1
### 				else:
### 					bin_ano_pos[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc > threshold_loc_1
### 				
### 			if per_both		== 'b': #in case we want to analyse both neg and positive percentiles
### 				threshold_loc_2	= np.mean(Results['all_pft']['threshold_ts2'][i])
### 				if plt_global_ts == 'pgts' or plt_freq == 'pf' :
### 					if 100-per <=50:
### 						bin_ano_neg[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc	< threshold_loc_2
### 					else:
### 						bin_ano_pos[idx_loc[0]:idx_loc[-1]+1,:,:] = ano_loc > threshold_loc_2
### 	
### 	if plt_global_ts == 'pgts' and per_both == 'b':
### 		gpp_loss = bin_ano_neg * ano_gC_pft
### 		gpp_gain = bin_ano_pos * ano_gC_pft
### 	elif plt_global_ts == 'pgts' and per_both != 'b':
### 		if per <= 50:
### 			gpp_loss = bin_ano_neg * ano_gC_pft
### 		else:
### 			gpp_gain = bin_ano_pos * ano_gC_pft
### 	
### 	#calculation for the global time series if the 4th arg/plt_global_ts = 'pgts'
### 	if plt_global_ts == 'pgts' and per_both =='b':
### 		global_gpp_loss	= []
### 		global_gpp_gain	= []
### 		freq_loss		= []
### 		freq_gain		= []
### 		for i in range(Results['all_pft']['dates_win'].flatten().size):
### 			global_gpp_loss .append(np.sum(gpp_loss[i]))	#global gpp loss TS, extremes only, adding over all the spatial points for every timestamp
### 			global_gpp_gain .append(np.sum(gpp_gain[i]))	#global gpp gain TS, extremes only
### 			freq_loss    	.append(np.sum(bin_ano_neg[i]))	#count of gpp loss events
### 			freq_gain		.append(np.sum(bin_ano_pos[i]))	##count of gpp gain events
### 		""" calculation of plotting the real magnitudes of global gpp loss and trend"""
### 	 	slope_loss, intercept_loss,_,_,_  = stats.linregress(time[...][:Results['all_pft']['dates_win'].flatten().size],global_gpp_loss)
### 		trend_loss       			=   slope_loss*time[...][:Results['all_pft']['dates_win'].flatten().size]+intercept_loss
### 		increase_loss    			=   (trend_loss[-1]-trend_loss[0])*100/trend_loss[0]
### 	
### 		""" calculation of plotting the real magnitudes of global gpp gain and trend"""
### 		slope_gain, intercept_gain,_,_,_  = stats.linregress(time[...][:Results['all_pft']['dates_win'].flatten().size],global_gpp_gain)
### 		trend_gain	      			=   slope_gain*time[...][:Results['all_pft']['dates_win'].flatten().size]+intercept_gain
### 		increase_gain    			=   (trend_gain[-1]-trend_gain[0])*100/trend_gain[0]
### 
### 	elif plt_global_ts == 'pgts'and per_both !='b':
###  		for i in range(dates_win.flatten().size):
### 			if per <= 50:
### 				global_gpp_loss .append(np.sum(gpp_loss[i]))    #global gpp loss TS, extremes only
### 				freq_loss       .append(np.sum(bin_ano_neg[i])) #count of gpp loss events
### 				""" calculation of plotting the real magnitudes of global gpp loss and trend"""
### 				slope_loss, intercept_loss,_,_,_  = stats.linregress(time[...][:Results['all_pft']['dates_win'].flatten().size],global_gpp_loss)
### 				trend_loss       			=   slope_loss*time[...][:Results['all_pft']['dates_win'].flatten().size]+intercept_loss
### 				increase_loss    			=   (trend_loss[-1]-trend_loss[0])*100/trend_loss[0]
### 				trend_gain, slope_gain		=	None, None
### 			else:
### 				global_gpp_gain .append(np.sum(gpp_gain[i]))    #global gpp gain TS, extremes only 
### 				freq_gain       .append(np.sum(bin_ano_pos[i])) ##count of gpp gain events
### 				""" calculation of plotting the real magnitudes of global gpp gain and trend"""
### 				slope_gain, intercept_gain,_,_,_  = stats.linregress(time[...][:Results['all_pft']['dates_win'].flatten().size],global_gpp_gain)
### 				trend_gain	      			=   slope_gain*time[...][:Results['all_pft']['dates_win'].flatten().size]+intercept_gain
### 				increase_gain    			=   (trend_gain[-1]-trend_gain[0])*100/trend_gain[0]
### 				trend_loss, slope_loss		=	None, None
### 
### 	if pft_comp	== 'comp':
### 		#plotting the timeseries for window years                                                     
### 		global_loss_tsw = None; global_gain_tsw = None; trend_loss_tsw  = None; trend_gain_tsw  = None
### 		global_loss_tsw = []
### 		global_gain_tsw = []
### 		trend_loss_tsw  = []
### 		trend_gain_tsw  = []
### 		slope_loss_tsw  = []
### 		slope_gain_tsw  = []
### 		for i in range(Results['all_pft']['idx_dates_win'].shape[0]):
### 			global_loss_loc = global_gpp_loss[Results['all_pft']['idx_dates_win'][i][0]:Results['all_pft']['idx_dates_win'][i][-1]+1]
### 			slope_loss_loc,intercept_loss_loc,_,_,_ = stats.linregress(time[...][Results['all_pft']['idx_dates_win'][i][0]:Results['all_pft']['idx_dates_win'][i][-1]+1],global_loss_loc)
### 			trend_loss_loc  = slope_loss_loc*time[...][Results['all_pft']['idx_dates_win'][i][0]:Results['all_pft']['idx_dates_win'][i][-1]+1] + intercept_loss_loc   
### 			increase_loss   = (trend_loss_loc[-1]-trend_loss_loc[0])*100/trend_loss_loc[0]
### 		
### 			global_gain_loc = global_gpp_gain[Results['all_pft']['idx_dates_win'][i][0]:Results['all_pft']['idx_dates_win'][i][-1]+1]
### 			slope_gain_loc,intercept_gain_loc,_,_,_ = stats.linregress(time[...][Results['all_pft']['idx_dates_win'][i][0]:Results['all_pft']['idx_dates_win'][i][-1]+1],global_loss_loc)
### 			trend_gain_loc  = slope_gain_loc*time[...][Results['all_pft']['idx_dates_win'][i][0]:Results['all_pft']['idx_dates_win'][i][-1]+1] + intercept_gain_loc
### 			increase_gain   = (trend_gain_loc[-1]-trend_gain_loc[0])*100/trend_gain_loc[0]
### 
### 			global_loss_tsw .append(global_loss_loc)
### 			global_gain_tsw .append(global_gain_loc)
### 			trend_loss_tsw  .append(trend_loss_loc) 
### 			trend_gain_tsw  .append(trend_gain_loc) 
### 			slope_loss_tsw  .append(slope_loss_loc) 
### 			slope_gain_tsw  .append(slope_gain_loc) 
### 
### 	return global_gpp_loss,global_gpp_gain,slope_loss,slope_gain,trend_loss,trend_gain, global_loss_tsw,global_gain_tsw,trend_loss_tsw,trend_gain_tsw, slope_loss_tsw, slope_gain_tsw
### 
### 
### if pft_comp == 'comp':
### 	Results_var_name_pft_comp   = ['global_gpp_loss', 'global_gpp_gain','slope_loss','slope_gain','trend_loss','trend_gain', 'global_loss_tsw','global_gain_tsw','trend_loss_tsw','trend_gain_tsw','slope_loss_tsw','slope_gain_tsw'] 
### 	if pft_idx_in in pft_idx[1:-2]:
### 		Results_calc_pft_comp       = Calc_pft_comp(pft_idx_in)
### 		for i in range(len(Results_var_name_pft_comp)):
### 			Results['comp_pft'][Results_var_name_pft_comp[i]]	= Results_calc_pft_comp[i]
### 	if pft_idx_in == 99:
### 	 	sort_pft = {}
### 	 	for p in range(len(pft_names[1:-2])):			#ignoring the Crop2 because its all nans
### 			Results_calc_pft_comp   = Calc_pft_comp(p+1)
### 			sort_pft[pft_names[p+1]]	= {}
### 			for i,var_name in enumerate (Results_var_name_pft_comp):
### 				sort_pft[pft_names[p+1]][var_name] = Results_calc_pft_comp[i] 		
### 		#Sorting the keys of positive and negative extremes by mean
### 		keys_pos_ext	= sort_pft.keys()
### 		keys_neg_ext	= sort_pft.keys()
### 		keys_neg_ext	. sort(key = lambda x: np.array(sort_pft[x]['global_gpp_loss']).mean())
### 		keys_pos_ext    . sort(key = lambda x: np.array(sort_pft[x]['global_gpp_gain']).mean(),reverse = True)
### 		"""
### 		Making a cummulative sum dictionary for the negative and postive extremes,
### 		The extremes are sorted in decending order of their mean losses,
### 		The cummulative sum is done in the same order as extremes, """
### 	
### 		#Negative extremes:
### 		ts 		= np.zeros(Results['all_pft']['dates_win'].flatten().shape[0])
### 		for key in keys_neg_ext:
### 			ts 	= ts + np.array(sort_pft[key]['global_gpp_loss'])
### 			Results['cumsum_neg_pft'][key]	= ts
### 		#Positive extremes
### 		ts      = np.zeros(Results['all_pft']['dates_win'].flatten().shape[0])
### 		for key in keys_pos_ext:
### 			ts  = ts + np.array(sort_pft[key]['global_gpp_gain'])
### 			Results['cumsum_pos_pft'][key]  = ts
### 		
### 		#if you want to see the consicutive window plots of pft contribution using pft_comp_win = 'compwin'
### 		if pft_comp_win == 'compwin':
### 			sort_win_pft = {}
### 			for win in range(Results['all_pft']['dates_win'].shape[0]):
### 				sort_win_pft[str(win)]={}
### 				for p in range(len(pft_names[1:-2])): #ignoring the Crop2 because its all nans
### 					Results_calc_pft_comp   = Calc_pft_comp(p+1)
### 					sort_win_pft[str(win)][pft_names[p+1]]  = {}
### 					for i,var_name in enumerate (Results_var_name_pft_comp):
### 						sort_win_pft[str(win)][pft_names[p+1]][var_name] = Results_calc_pft_comp[i]         
### 				#Sorting the keys of positive and negative extremes by mean
### 				keys_pos_ext    = sort_win_pft[str(win)].keys()
### 				keys_neg_ext    = sort_win_pft[str(win)].keys()
### 				keys_neg_ext    . sort(key = lambda x: np.array(sort_win_pft[str(win)][x]['global_loss_tsw']).mean())
### 				keys_pos_ext    . sort(key = lambda x: np.array(sort_win_pft[str(win)][x]['global_gain_tsw']).mean(),reverse = True)
### 					#Negative Extremes
### 				ts	= np.zeros(Results['all_pft']['dates_win'][win].shape[0])
### 				for key in keys_neg_ext:
### 					ts  = ts + np.array(sort_win_pft[str(win)][key]['global_loss_tsw'])
### 					Results ['compwin_cumsum_neg']  = ts
### 					
### 					#Positive Extremes
### 				ts	= np.zeros(Results['all_pft']['dates_win'][win].shape[0])
### 				for key in keys_pos_ext:
### 					ts  = ts + np.array(sort_win_pft[str(win)][key]['global_gain_tsw'])
### 					Results ['compwin_cumsum_pos'] 	= ts
### 
### def Tables():
### 	if pft_comp == 'comp':
### 		#Capturing the slope of positive and negative extremes for multiple time windows
### 		idx				= np.array([Results['all_pft']['dates_win'][i][0].year for i in range(Results['all_pft']['dates_win'].shape[0])])
### 		#data_df_tsw		= {'Slope_Neg_Ext(gC/yr)' : Results['all_pft']['slope_loss_tsw'], 'Slope_Pos_Ext(gC/yr)' : Results['all_pft']['slope_gain_tsw'], 'Slope_Neg_Ext(kgC/dec)' : np.array(Results['all_pft']['slope_loss_tsw'])/100, 'Slope_Pos_Ext(kgC/dec)' : np.array(Results['all_pft']['slope_gain_tsw'])/100}
### 		data_df_tsw		= {'Slope_Neg_Ext(kgC/decade)' : np.array(Results['all_pft']['slope_loss_tsw'])/100, 'Slope_Pos_Ext(kgC/decade)' : np.array(Results['all_pft']['slope_gain_tsw'])/100}
### 		df_tsw			= pd.DataFrame(data_df_tsw)
### 		df_tsw.index	= idx
### 		acc_neg,int_acc_neg,_,_,_ 	= stats.linregress(idx,df_tsw['Slope_Neg_Ext(kgC/decade)'])	#negative means increasing acceleration for negative extremes g/decade/decade
### 		acc_pos,int_acc_pos,_,_,_  = stats.linregress(idx,df_tsw['Slope_Pos_Ext(kgC/decade)'])
### 
### 		fig_t1			= plt.figure()
### 		df_tsw			. plot.bar()
### 		plt.title		("Acc of negative and positive extremes are  : %.2f  and %.2f gC/sq.decade rest" % (acc_neg*1000, acc_pos*1000))
### 		plt.ylabel		("kgC/decade")
### 		plt.xlabel		("decades")
### 		
### 		if pft_idx_in in pft_idx[1:-2]:
### 			data_df_tsw_pft	= {'Slope_Neg_Ext(kgC/decade)_pft' : np.array(Results['sel_pft']['slope_loss_tsw'])/100, 'Slope_Pos_Ext(kgC/decade)_pft' : np.array(Results['sel_pft']['slope_gain_tsw'])/100}
### 			df_tsw			= pd.DataFrame(data_df_tsw_pft)
### 			df_tsw.index	= idx
### 			acc_neg_pft,int_acc_neg_pft,_,_,_ 	= stats.linregress(idx,df_tsw['Slope_Neg_Ext(kgC/decade)_pft'])	#negative means increasing acceleration for negative extremes g/dec/decade
### 			acc_pos_pft,int_acc_pos_pft,_,_,_  = stats.linregress(idx,df_tsw['Slope_Pos_Ext(kgC/decade)_pft'])
### 
### 			fig_t2			= plt.figure()
### 			df_tsw			. plot.bar()
### 			plt.title		("Acc of negative and positive extremes of a selected pft are  : %.2f  and %.2f gC/sq.decade rest" % (acc_neg_pft*1000, acc_pos_pft*1000))
### 			plt.ylabel		("kgC/decade")
### 			plt.xlabel		("decades")
### 
### 		plt.show()
### 
### 
### def Ploting():
### 	print ("Plotting Thresholds......")
### 	#ploting the threshold
### 	if plt_threshold    == 'pth':
### 		fig1            = plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
### #		plt.title ("The changes in the threshold values of anomalies when per is %.2f" %per)
### 		if per_both     == 'b':
### 			plt.plot    (Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts1'].flatten()/10**9,'r',Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts2'].flatten()/10**9,'g',Results['sel_pft']['dates_win'].flatten(),[0]*len(Results['sel_pft']['dates_win'].flatten()),'k')
### 			ymax = np.max(np.absolute(Results['sel_pft']['threshold_ts1'].flatten()/10**9))  #maximum absolute value for ploting the y axis
### 			ymin = np.min(np.absolute(Results['sel_pft']['threshold_ts1'].flatten()/10**9))  #maximum absolute value for ploting the y axis
### 			if ymax < np.max(np.absolute(Results['sel_pft']['threshold_ts2'].flatten()/10**9)):
### 				ymax = np.max(np.absolute(Results['sel_pft']['threshold_ts2'].flatten()/10**9))
### 			if ymin > np.min(np.absolute(Results['sel_pft']['threshold_ts2'].flatten()/10**9)):
### 				ymin = np.min(np.absolute(Results['sel_pft']['threshold_ts2'].flatten()/10**9))
### 			ymax_ll = round(-ymax/100.,0)*100   # lower limit on y axis
### 			ymax_ul = abs(ymax_ll)+100  # lower limit on y axis
### 			
### 			if pft_comp == 'comp' : plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['threshold_ts1'].flatten()/10**9,'royalblue',Results['all_pft']['dates_win'].flatten(),Results['all_pft']['threshold_ts2'].flatten()/10**9,'royalblue')
### 		else:
### 			plt.plot    (Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts1'].flatten()/10**9,'r',Results['sel_pft']['dates_win'].flatten(),[0]*len(Results['sel_pft']['dates_win'].flatten()),'k')
### 			if pft_comp == 'comp' : plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts1'].flatten()/10**9,'r',Results['sel_pft']['dates_win'].flatten(),[0]*len(Results['sel_pft']['dates_win'].flatten()),'k')
###    		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('G'+ano.units, fontsize = 16)
### 
### 		tmp_idx = np.arange(0,5401,600) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		#-------------------------2
### 		#Hard coding the y limits
### 		ymax_ll	= 400
###    		ymax_ul	= 850
### 		plt.ylim(400,850)
### 		#-------------------------2
### 		
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = Results['sel_pft']['dates_win'].flatten()[i]
### 			dates_ticks.append(a)
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		if per_both == 'b':
### 			plt.yticks(np.arange(ymax_ll,ymax_ul,100),fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		fig1.savefig(in_path+'results/timeseries_threshold_per_%d.png' %(per))
###    		fig1.savefig(in_path+'results/timeseries_threshold_per_%d.pdf' %(per))
### 		
### #ploting the threshold with sum line
### 		fig101            = plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
### 		plt.plot    (Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts1'].flatten()/10**9,'r',Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts2'].flatten()/10**9,'g',Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts1'].flatten()/10**9+Results['sel_pft']['threshold_ts2'].flatten()/10**9,'k')
### 		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('G'+ano.units, fontsize = 16)
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		plt.yticks(np.arange(ymax_ll,ymax_ul,100),fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		fig101.savefig(in_path+'results/timeseries_threshold_sum_per_%d.png' %(per))
### 		fig101.savefig(in_path+'results/timeseries_threshold_sum_per_%d.pdf' %(per))
### 		
### 
### #Ploting the Threshold with secondary axis
### 		fig111,ax2=plt.subplots()
### 		#----------------------------1
### 		#Hard coding the y limits
### 		ymin	= 400
###    		ymax	= 850
### 		#----------------------------1
### 		ax2.plot(Results['sel_pft']['dates_win'].flatten(),np.abs(Results['sel_pft']['threshold_ts1'].flatten())/10**9,'r')
### 		ax2.set_ylabel("Negative Extremes (G%s)"%ano.units, {'color': 'r'},fontsize =14)
### 		ax2.set_xlabel("Time", fontsize = 14)
### 		ax2.set_ylim([ymin,ymax])
### 		ax2.set_yticks(np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100+50),100))
### 		ax2.set_xticks(dates_ticks)
### 		ax2.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		ax1=ax2.twinx()
### 		ax1.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['threshold_ts2'].flatten()/10**9,'g')
### 		ax1.set_ylabel(" Positive Extremes (G%s)"%ano.units,{'color': 'g'}, fontsize =14)
### 		ax1.set_ylim([ymin,ymax])   
### 		ax1.set_yticks(np.arange(int(np.floor(ymin/100)*100),int(np.ceil(ymax/100)*100+50),100))
### 		ax1.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		fig111.savefig(in_path+'results/timeseries_secaxis_threshold_per_%d.pdf' %(per))
### 		fig111.savefig(in_path+'results/timeseries_secaxis_threshold_per_%d.png' %(per))
### 		
### 		for idx in np.arange(0,5400,300):
### 			print "neg: ",idx, ":" , (Results['sel_pft']['threshold_ts1'].flatten()/10**9)[idx]
### 			print "pos: ",idx, ":" , (Results['sel_pft']['threshold_ts2'].flatten()/10**9)[idx]
### 		#	plt.show()
### 	#ploting the global timeseries 
### 	if plt_global_ts == 'pgts':
### 		fig2	= plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
### 		if per_both == 'b':
### #			plt.title("Global Time Series of pos and neg extremes when percentile = %.1f: slope (%.1f & %.1f) kgC/yr  and impact (%.1f & %.1f)%% resp"%(per,Results['sel_pft']['slope_gain']*.365,Results['sel_pft']['slope_loss']*.365,Impact(Results['sel_pft']['trend_gain']),Impact(Results['sel_pft']['trend_loss'])))
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_loss']/10**15, 	'r',label = 'sel pft neg extremes',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_loss']/10**15,		'k--',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_gain']/10**15, 	'g' ,label = 'sel pft pos extremes',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_gain']/10**15,		'k--',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = 0.5)
### 			if pft_comp == 'comp':
### 				plt.title("Global Time Series of pos and neg extremes of the pft : %s  when percentile = %.1f: slope (%.1f & %.1f) kgC/yr  and impact (%.1f & %.1f)%% resp \n All_pft slope (%.1f & %.1f)kgC/yr  and impact (%.1f & %.1f)%% "%(pft_names[pft_idx_in],per,Results['sel_pft']['slope_gain']*.365,Results['sel_pft']['slope_loss']*.365,Impact(Results['sel_pft']['trend_gain']),Impact(Results['sel_pft']['trend_loss']),Results['all_pft']['slope_gain']*.365,Results['all_pft']['slope_loss']*.365,Impact(Results['all_pft']['trend_gain']),Impact(Results['all_pft']['trend_loss'])),fontsize = 16 )
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['global_gpp_loss']/10**15, 	'orangered',label = 'all pft neg extremes',linewidth = 0.5)
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['trend_loss']/10**15,		'b--',linewidth = 0.5)
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['global_gpp_gain']/10**15, 	'lightgreen', label = 'all pft pos extremes',linewidth = 0.5)
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['trend_gain']/10**15,		'b--',linewidth = 0.5)
### 				plt.legend()
### 
### 				
### 		else:
### 			if per <50:
### #				plt.title("Global Time Series of neg extremes when percentile = %.1f slope : %.1f kgC/yr and impact :%.1f %%" %(per,Results['sel_pft']['slope_loss']*.365, Impact(Results['sel_pft']['trend_loss'])))
### 				plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_loss']/10**15, 	'r',
### 						 Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_loss']/10**15,		'k--',
### 						 Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k--',linewidth = .5)
### 				plt.xlabel('Time', fontsize = 12)
### 				plt.ylabel(r'$10^{9} gC$', fontsize = 12)
### 
### 			else:
### 				plt.title("Global Time Series of pos extremes when percentile = %.1f slope : %.1f kgC/yr and impact :%.1f %%"%(per,Results['sel_pft']['slope_gain']*.365,Impact(Results['sel_pft']['trend_gain'])))
### 				plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_gain']/10**15, 	'g',                              
### 				 		 Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_gain']/10**15,		'k--',                                  
### 						 Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = 0.5) 
### 		tmp_idx = np.arange(0,5401,600) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = Results['sel_pft']['dates_win'].flatten()[i]
### 			dates_ticks.append(a)
### 		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('P'+ano.units, fontsize = 16)
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		plt.yticks(fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		
### 		plt.text(Results['sel_pft']['dates_win'].flatten()[600],.8,"Slope = %d MgC/month"%(Results['sel_pft']['slope_gain']*.365/10**6),size = 14, color = 'g')
### 		plt.text(Results['sel_pft']['dates_win'].flatten()[600],-.9,"Slope = %d MgC/month"%(Results['sel_pft']['slope_loss']*.365/10**6),size = 14, color = 'r')
### 		fig2.savefig(in_path+'results/global_%s_timeseries_per_%d.png' %(variable,per))
### 		fig2.savefig(in_path+'results/global_%s_timeseries_per_%d.pdf' %(variable,per))
### 
### 	#ploting the global timeseries with slopes from 1850-2100 and 2100-2299
### 	print ("Plotting the global timeseries with slopes from 1850-2100 and 2100-2299")
### 	if plt_global_ts == 'pgts':
### 		fig201	= plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
### 		if per_both == 'b':
### #			plt.title("Global Time Series of pos and neg extremes when percentile = %.1f: slope (%.1f & %.1f) kgC/yr  and impact (%.1f & %.1f)%% resp"%(per,Results['sel_pft']['slope_gain']*.365,Results['sel_pft']['slope_loss']*.365,Impact(Results['sel_pft']['trend_gain']),Impact(Results['sel_pft']['trend_loss'])))
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_loss']/10**15, 	'r',label = 'sel pft neg extremes',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten()[:idx_yr_2100],Results['sel_pft']['trend_loss_21']/10**15,		'k--',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten()[idx_yr_2100:idx_yr_2300],Results['sel_pft']['trend_loss_23']/10**15,		'k--',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_gain']/10**15, 	'g' ,label = 'sel pft pos extremes',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten()[:idx_yr_2100],Results['sel_pft']['trend_gain_21']/10**15,		'k--',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten()[idx_yr_2100:idx_yr_2300],Results['sel_pft']['trend_gain_23']/10**15,		'k--',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = 0.5)
### 			if pft_comp == 'comp':
### 				plt.title("Global Time Series of pos and neg extremes of the pft : %s  when percentile = %.1f: slope (%.1f & %.1f) kgC/yr  and impact (%.1f & %.1f)%% resp \n All_pft slope (%.1f & %.1f)kgC/yr  and impact (%.1f & %.1f)%% "%(pft_names[pft_idx_in],per,Results['sel_pft']['slope_gain']*.365,Results['sel_pft']['slope_loss']*.365,Impact(Results['sel_pft']['trend_gain']),Impact(Results['sel_pft']['trend_loss']),Results['all_pft']['slope_gain']*.365,Results['all_pft']['slope_loss']*.365,Impact(Results['all_pft']['trend_gain']),Impact(Results['all_pft']['trend_loss'])),fontsize = 16 )
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['global_gpp_loss']/10**15, 	'orangered',label = 'all pft neg extremes',linewidth = 0.5)
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['trend_loss']/10**15,		'b--',linewidth = 0.5)
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['global_gpp_gain']/10**15, 	'lightgreen', label = 'all pft pos extremes',linewidth = 0.5)
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['trend_gain']/10**15,		'b--',linewidth = 0.5)
### 				plt.legend()
### 
### 				
### 		else:
### 			if per <50:
### #				plt.title("Global Time Series of neg extremes when percentile = %.1f slope : %.1f kgC/yr and impact :%.1f %%" %(per,Results['sel_pft']['slope_loss']*.365, Impact(Results['sel_pft']['trend_loss'])))
### 				plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_loss']/10**15, 	'r',
### 						 Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_loss']/10**15,		'k--',
### 						 Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k--',linewidth = .5)
### 				plt.xlabel('Time', fontsize = 12)
### 				plt.ylabel(r'$10^{9} %s$'%ano.units, fontsize = 12)
### 
### 			else:
### 				plt.title("Global Time Series of pos extremes when percentile = %.1f slope : %.1f kgC/yr and impact :%.1f %%"%(per,Results['sel_pft']['slope_gain']*.365,Impact(Results['sel_pft']['trend_gain'])))
### 				plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_gain']/10**15, 	'g',                              
### 				 		 Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_gain']/10**15,		'k--',                                  
### 						 Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = 0.5) 
### 		tmp_idx = np.arange(0,5401,600) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = Results['sel_pft']['dates_win'].flatten()[i]
### 			dates_ticks.append(a)
### 		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('P'+ano.units, fontsize = 16)
### 
### 		#-------------------------3
### 		#Hard coding the y limits
### 		if variable == 'gpp':
### 			ymin	= -.9
###    			ymax	= .9
### 			plt.ylim(ymin,ymax)
### 		#-------------------------3
### 
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		plt.yticks(fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		plt.text(Results['sel_pft']['dates_win'].flatten()[600], ymax - 0.1,"Slope = %d MgC/month"%(Results['sel_pft']['slope_gain_21']*.365/10**6),size = 14, color = 'g')
### 		plt.text(Results['sel_pft']['dates_win'].flatten()[600+idx_yr_2100], ymax-0.1,"Slope = %d MgC/month"%(Results['sel_pft']['slope_gain_23']*.365/10**6),size = 14, color = 'g')
### 		plt.text(Results['sel_pft']['dates_win'].flatten()[600], ymin+0.1,"Slope = %d MgC/month"%(Results['sel_pft']['slope_loss_21']*.365/10**6),size = 14, color = 'r')
### 		plt.text(Results['sel_pft']['dates_win'].flatten()[600+idx_yr_2100], ymin+0.1,"Slope = %d MgC/month"%(Results['sel_pft']['slope_loss_23']*.365/10**6),size = 14, color = 'r')
### 		fig201.savefig(in_path+'results/global_%s_timeseries_per_%d_changing_slopes.png' %(variable,per))
### 		fig201.savefig(in_path+'results/global_%s_timeseries_per_%d_changing_slopes.pdf' %(variable,per))
### 		#print in_path+'results/global_%s_timeseries_per_%d_changing_slopes.pdf' %(variable,per)
### 	#ploting the global cumulative sum timeseries with values at 2100 and 2300 
### 	print ("Plotting the global cumsum timeseries with values at 1850,2100 and 2299")
### 	if plt_global_ts == 'pgts':
### 		fig211	= plt.figure(tight_layout = True, figsize = (9,5), dpi = 400)
### 		if per_both == 'b':
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),np.cumsum(Results['sel_pft']['global_gpp_loss']/10**15), 	'r',label = 'sel pft neg extremes',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),np.cumsum(Results['sel_pft']['global_gpp_gain']/10**15), 	'g' ,label = 'sel pft pos extremes',linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = 0.5)
### 				
### 		tmp_idx = np.arange(0,5401,600) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = Results['sel_pft']['dates_win'].flatten()[i]
### 			dates_ticks.append(a)
### 		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('P'+ano.units, fontsize = 16)
### 
### 		#-------------------------3
### 		#Hard coding the y limits
### 		if variable == 'gpp':
### 			ymin	= -1300
###    			ymax	=  1100
### 			plt.ylim(ymin,ymax)
### 		#-------------------------3
### 
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		plt.yticks(fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		fig211.savefig(in_path+'results/global_%s_cumsum_timeseries_per_%d_changing_slopes.png' %(variable,per))
### 		fig211.savefig(in_path+'results/global_%s_cumsum_timeseries_per_%d_changing_slopes.pdf' %(variable,per))
### 
### 	#ploting the global yearly and 5yr rolling mean sum timeseries 
### 	print ("Plotting the global yearly Mean and 5yr Roalling mean timeseries with values at 1850,2100 and 2299")
### 	if plt_global_ts == 'pgts':
### 		fig221	= plt.figure(tight_layout = True, figsize = (10,5), dpi = 400)
### 		if per_both == 'b':
### 			yr_global_tot_negext	= np.array([np.sum(Results['sel_pft']['global_gpp_loss'][i*12:(i*12)+12]) for i in range(len(range(1850,2300)))])
### 			yr_global_tot_posext	= np.array([np.sum(Results['sel_pft']['global_gpp_gain'][i*12:(i*12)+12]) for i in range(len(range(1850,2300)))])
### 			global_tot_negext_df	= pd.Series(yr_global_tot_negext)
### 			global_tot_posext_df	= pd.Series(yr_global_tot_posext)
### 			rm_5yr_global_tot_negext= global_tot_negext_df.rolling(window =5, center = True).mean()
### 			rm_5yr_global_tot_posext= global_tot_posext_df.rolling(window =5, center = True).mean()
### 		
### 		yrs                 = np.array(range(1850,2300))
### 		plt.plot(yrs, yr_global_tot_negext/10**15,	'r', label = 'neg extremes', linewidth = 0.5)
### 		plt.plot(yrs, yr_global_tot_posext/10**15,	'g', label = 'pos extremes', linewidth = 0.5)
### 		plt.plot(yrs,[0]*yr_global_tot_negext.flatten().size,'k',linewidth = 0.5)
### 
### 		tmp_idx = np.arange(0,450,50) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = yrs[i]
### 			dates_ticks.append(a)
### 		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('P'+ano.units, fontsize = 16)
### 
### 		#-------------------------3
### 		#Hard coding the y limits
### 		ymin	= -5
###    		ymax	=  5
### 		plt.ylim(ymin,ymax)
### 		#-------------------------3
### 
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		plt.yticks(fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 		fig221.savefig(in_path+'results/global_%s_yr_timeseries_per_%d_changing_slopes.png' %(variable,per))
### 		fig221.savefig(in_path+'results/global_%s_yr_timeseries_per_%d_changing_slopes.pdf' %(variable,per))
### 		
### 		#Rolling Mean
### 		fig222	= plt.figure(tight_layout = True, figsize = (10,5), dpi = 400)
### 		plt.plot(yrs, rm_5yr_global_tot_negext/10**15,	'r', label = 'neg extremes', linewidth = 0.5)
### 		plt.plot(yrs, rm_5yr_global_tot_posext/10**15,	'g', label = 'pos extremes', linewidth = 0.5)
### 		plt.plot(yrs,[0]*yr_global_tot_negext.flatten().size,'k',linewidth = 0.5)
### 		plt.xlabel('Time', fontsize = 16)
### 		plt.ylabel('P'+ano.units, fontsize = 16)
### 
### 		#-------------------------4
### 		#Hard coding the y limits
### 		ymin	= -5
###    		ymax	=  5
### 		plt.ylim(ymin,ymax)
### 		#-------------------------4
### 		plt.xticks(dates_ticks,fontsize = 14) 
### 		plt.yticks(fontsize = 14)
### 		plt.grid(which='major', linestyle=':', linewidth='0.3', color='gray')
### 
### 		fig221.savefig(in_path+'results/global_%s_5yr_rm_timeseries_per_%d_changing_slopes.png' %(variable,per))
### 		fig221.savefig(in_path+'results/global_%s_5yr_rm_timeseries_per_%d_changing_slopes.pdf' %(variable,per))
### 
### 		#plt.show()
### 		#ploting the counts of gpp gain, loss and ratio when the threshold is defined as average of first five windows of 30 years i.e. 1850-1999
### 	print ("plotting the frequency of extreme events relative to threshold defined from 1850-1999")
### 	if per_both == 'b':
### 		fig3, ax    = plt.subplots(nrows=3,ncols= 2 if pft_comp == 'comp' else 1,tight_layout = True, figsize = (9,7), dpi = 400)
### 		ax          = ax.ravel()
### 		tmp_idx = np.arange(0,5401,600) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = Results['sel_pft']['dates_win'].flatten()[i]
### 			dates_ticks.append(a)
### 		#plt.suptitle("Count of extreme events relative to the threshold 1850-1999 when per is %.2f (left: sel_pft, right:all_pfts) " %(per))
### 		ax[0].plot(dates_ar, Results['sel_pft']['freq_gain_r1999'],'g--',linewidth =.5)
### 		#ax[0].set_title('Count of pos gpp extremes with slope of %.2f events per decade' %(Results['sel_pft']['slope_gain_r1999']*365*10))
### 		ax[0].plot(dates_ar, Results['sel_pft']['trend_gain_r1999'], 'k--',linewidth =.5)
### 		ax[0].set_ylabel('Positive \nExtremes', fontsize = 16)
### 		ax[0].set_xticks(dates_ticks)
### 		ax[0].tick_params(labelsize=12)
### 		ax[0].set_xticklabels([])
### 		ax[0].text(Results['sel_pft']['dates_win'].flatten()[300],2000,"Slope = %.2f events/decade"%(Results['sel_pft']['slope_gain_r1999']*3650),size = 14, color = 'g')
### 		ax[2 if pft_comp == 'comp' else 1].plot(dates_ar, Results['sel_pft']['freq_loss_r1999'],'r--',linewidth =.5)
### 		#ax[2 if pft_comp == 'comp' else 1].set_title('Count of neg gpp extrems with slope of %.2f events per decade' %(Results['sel_pft']['slope_loss_r1999']*365*10))   
### 		ax[2 if pft_comp == 'comp' else 1].plot(dates_ar, Results['sel_pft']['trend_loss_r1999'], 'k--',linewidth =.5) 
### 		ax[2 if pft_comp == 'comp' else 1].set_ylabel('Negative \nExtremes', fontsize = 16)  
### 		ax[2 if pft_comp == 'comp' else 1].set_xticks(dates_ticks)
### 		ax[2 if pft_comp == 'comp' else 1].set_xticklabels([])
### 		ax[2 if pft_comp == 'comp' else 1].tick_params(labelsize=12)
### 		ax[2 if pft_comp == 'comp' else 1].text(Results['sel_pft']['dates_win'].flatten()[300],2000,"Slope = %.2f events/decade"%(Results['sel_pft']['slope_loss_r1999']*3650),size = 14, color = 'r')
### 		ax[4 if pft_comp == 'comp' else 2].plot(dates_ar, np.array(Results['sel_pft']['freq_loss_r1999'])/np.array(Results['sel_pft']['freq_gain_r1999']),'skyblue',linewidth = 0.5)
### 		#ax[4 if pft_comp == 'comp' else 1].set_title('Ratio of counts of neg to pos extremes')
### 		slope_ratio,intercept_ratio,_,_,_ = stats.linregress(time[...],(np.array(Results['all_pft']['freq_loss_r1999'])/np.array(Results['all_pft']['freq_gain_r1999'])))
### 		s,c,_,_,_=stats.linregress(time[...], np.array(Results['sel_pft']['freq_loss_r1999'])/np.array(Results['sel_pft']['freq_gain_r1999']))
### 		#if 's' and 'c' are 'nan' use s=999 and c =999
### 		if s*0 != 0: s = 999
### 		if c*0 != 0: c = 999
### 		trend_ratio     =   slope_ratio*time[...]+intercept_ratio
### 		ax[4 if pft_comp == 'comp' else 2].plot(dates_ar,trend_ratio,'k--', linewidth =.6)
### 		ax[4 if pft_comp == 'comp' else 2].plot(dates_ar,[1]*len(dates_ar),'c--', linewidth =.2)
### 		ax[4 if pft_comp == 'comp' else 2].set_xticks(dates_ticks)#,fontsize = 14)
### 		ax[4 if pft_comp == 'comp' else 2].set_ylabel('Ratio of\nNegative to \nPositive \nExtremes', fontsize = 16) 
### 		ax[4 if pft_comp == 'comp' else 2].set_xlabel('Time', fontsize = 16) 
### 		ax[4 if pft_comp == 'comp' else 2].tick_params(labelsize=12)
### 		ax[4 if pft_comp == 'comp' else 2].text(Results['sel_pft']['dates_win'].flatten()[300],2,"Slope = %.1f x $10^{3}$ /decade"%(int(s*3650*10**3)),size = 14, color = 'gray')
### 		if pft_comp	== 'comp':
### 			ax[1].plot(dates_ar, Results['all_pft']['freq_gain_r1999'],'g--',linewidth =.5)
### 			ax[1].set_title('Count of pos gpp extremes with slope of %.2f events per decade'%(Results['all_pft']['slope_gain_r1999']*365*10))
### 			ax[1].plot(dates_ar, Results['all_pft']['trend_gain_r1999'], 'k--',linewidth =.5)
### 			ax[3].plot(dates_ar, Results['all_pft']['freq_loss_r1999'],'r--',linewidth =.5)
### 			ax[3].set_title('Count of neg gpp extremes with slope of %.2f events per decade' %(Results['all_pft']['slope_loss_r1999']*365*10))
### 			ax[3].plot(dates_ar, Results['all_pft']['trend_loss_r1999'], 'k--',linewidth =.5)  
### 			ax[5].plot(dates_ar, np.array(Results['all_pft']['freq_loss_r1999'])/np.array(Results['all_pft']['freq_gain_r1999']),'b',linewidth =.5)
### 			ax[5].set_title('Ratio of counts of neg to pos extremes')
### 			ax[5].plot(dates_ar,[1]*len(dates_ar),'k', linewidth =.8)
### 			
### 		fig3.savefig(in_path+'results/count_relative_to_threshold_1850-1999_per_%d.pdf' %(per))
### 		fig3.savefig(in_path+'results/count_relative_to_threshold_1850-1999_per_%d.png' %(per))
### 		#plt.show()
### 	else:
### 		if per <50:
### 			plt.title("Count of neg extreme events relative to the threshold 1850-1999 when per is %.2f" %(per))
### 			plt.plot(dates_ar, Results['sel_pft']['freq_loss_r1999'],'r--',linewidth =.5)
### 		else:
### 			plt.title("Count of pos extreme events relative to the threshold 1850-1999 when per is %.2f" %(per))
### 			plt.plot(dates_ar, Results['sel_pft']['freq_gain_r1999'],'r--',linewidth =.5)
### 
### 		fig3.savefig(in_path+'results/count_relative_to_threshold_1850-1999_per_%d.png' %(per))
### 		fig3.savefig(in_path+'results/count_relative_to_threshold_1850-1999_per_%d.pdf' %(per))
### 		#plt.show()
### 	
### #plotting counts with changing slope
### 	if per_both == 'b':
### 		fig301, ax    = plt.subplots(nrows=2,ncols= 1,tight_layout = True, figsize = (9,7), dpi = 400)
### 		ax          = ax.ravel()
### 		tmp_idx = np.arange(0,5401,600) #for x ticks
### 		tmp_idx[-1]=tmp_idx[-1]-1
### 		dates_ticks = []
### 		for i in tmp_idx:
### 			a = Results['sel_pft']['dates_win'].flatten()[i]
### 			dates_ticks.append(a)
### 		#plt.suptitle("Count of extreme events relative to the threshold 1850-1999 when per is %.2f (left: sel_pft, right:all_pfts) " %(per))
### 		ax[0].plot(dates_ar, Results['sel_pft']['freq_gain_r1999'],'g--',linewidth =.5)
### 		#ax[0].set_title('Count of pos gpp extremes with slope of %.2f events per decade' %(Results['sel_pft']['slope_gain_r1999']*365*10))
### 		ax[0].plot(dates_ar[:idx_yr_2100], Results['sel_pft']['trend_gain_r1999_21'], 'k--',linewidth =.5)
### 		ax[0].plot(dates_ar[idx_yr_2100:idx_yr_2299], Results['sel_pft']['trend_gain_r1999_23'], 'k--',linewidth =.5)
### 		ax[0].set_ylabel('Positive \nExtremes', fontsize = 16)
### 		ax[0].set_xticks(dates_ticks)
### 		ax[0].tick_params(labelsize=12)
### 		ax[0].set_xticklabels([])
### 		ax[0].text(dates_ar[400], 1100 if variable == 'gpp' else np.max(Results['sel_pft']['freq_gain_r1999'])*.9,"Slope = %.2f events/decade"%(Results['sel_pft']['slope_gain_r1999_21']*3650),size = 12, color = 'g')
### 		ax[0].text(dates_ar[400+idx_yr_2100], 1100 if variable=='gpp' else np.max(Results['sel_pft']['freq_gain_r1999'])*.9,"Slope = %.2f events/decade"%(Results['sel_pft']['slope_gain_r1999_23']*3650),size = 12, color = 'g')
### 	
### 		ax[1].plot(dates_ar, Results['sel_pft']['freq_loss_r1999'],'r--',linewidth =.5)
### 		#ax[0].set_title('Count of pos gpp extremes with slope of %.2f events per decade' %(Results['sel_pft']['slope_gain_r1999']*365*10))
### 		ax[1].plot(dates_ar[:idx_yr_2100], Results['sel_pft']['trend_loss_r1999_21'], 'k--',linewidth =.5)
### 		ax[1].plot(dates_ar[idx_yr_2100:idx_yr_2299], Results['sel_pft']['trend_loss_r1999_23'], 'k--',linewidth =.5)
### 		ax[1].set_ylabel('Negative \nExtremes', fontsize = 16)
### 		ax[1].set_xlabel('Time', fontsize = 16)
### 		ax[1].set_xticks(dates_ticks)
### 		ax[1].tick_params(labelsize=12)
### 		#ax[1].set_xticklabels([])
### 		ax[1].text(dates_ar[400], 1100 if variable =='gpp' else np.max(Results['sel_pft']['freq_loss_r1999'])*.9,"Slope = %.2f events/decade"%(Results['sel_pft']['slope_loss_r1999_21']*3650),size = 12, color = 'r')
### 		ax[1].text(dates_ar[400+idx_yr_2100], 1100 if variable =='gpp' else np.max(Results['sel_pft']['freq_loss_r1999'])*.9,"Slope = %.2f events/decade"%(Results['sel_pft']['slope_loss_r1999_23']*3650),size = 12, color = 'r')
### 
### 		#-------------------------4
### 		#Hard coding the y limits
### 		if variable == 'gpp':
### 			ymin	= 0
###    			ymax	= 1250
### 			ax[0].set_ylim([ymin,ymax])
### 			ax[1].set_ylim([ymin,ymax])
### 		#-------------------------4
### 
### 
### 		fig301.savefig(in_path+'results/count_relative_to_threshold_1850-1999_changingslopes_per_%d.pdf' %(per))
### 		fig301.savefig(in_path+'results/count_relative_to_threshold_1850-1999_changingslopes_per_%d.png' %(per))
### 
### #	Plotting the per pixel timeseries	
### 	if -90 <= lat_in <=90 	and 	0 <=lon_in <=360:
### 		fig4 = plt.figure()
### 		if per_both == 'b':
### 			plt.plot(dates_ar[:-(time.size - Results['sel_pft']['dates_win'].size)],Results['sel_pft']['pp_ts'],'c',
### 					 dates_ar[:-(time.size - Results['sel_pft']['dates_win'].size)],Results['sel_pft']['threshold_ts1'].flatten(),'g',
### 					 dates_ar[:-(time.size - Results['sel_pft']['dates_win'].size)],Results['sel_pft']['threshold_ts2'].flatten(),'k',
### 					 linewidth =1)
### 		else:
### 			plt.plot(dates_ar[:-(time.size - Results['sel_pft']['dates_win'].size)],Results['sel_pft']['pp_ts'],'c',
### 					 dates_ar[:-(time.size - Results['sel_pft']['dates_win'].size)],Results['sel_pft']['threshold_ts1'].flatten(),'k',
### 					 linewidth= .8)	
### 		plt.show()
### 
### #plotting the timeseries for 50 yer window
### 	if p_win	== 	'pwin' and per_both== 'b':
### 		fig5,ax	=	plt.subplots (Results['sel_pft']['nrows_tsw'],Results['sel_pft']['ncols_tsw'])
### 		ax		=	ax.ravel()
### 		plt.suptitle ("Timeseries of %d year slices" %window)
### 		for i in range(Results['sel_pft']['dates_win'].shape[0]):
### 			ax[i].plot(Results['sel_pft']['dates_win'][i],Results['sel_pft']['global_loss_tsw'][i],'r',
### 					   Results['sel_pft']['dates_win'][i],Results['sel_pft']['trend_loss_tsw'][i],'k--',
### 					   Results['sel_pft']['dates_win'][i],Results['sel_pft']['global_gain_tsw'][i],'g',
### 					   Results['sel_pft']['dates_win'][i],Results['sel_pft']['trend_gain_tsw'][i],'k--',
### 					   Results['sel_pft']['dates_win'][i],[0]*Results['sel_pft']['dates_win'][i].size,'k',linewidth = 0.5)
### 			if pft_comp == 'comp':
### 				ax[i].plot(	Results['all_pft']['dates_win'][i],Results['all_pft']['global_loss_tsw'][i],'skyblue',
### 					   		Results['all_pft']['dates_win'][i],Results['all_pft']['trend_loss_tsw'][i],'b--',
### 					   		Results['all_pft']['dates_win'][i],Results['all_pft']['global_gain_tsw'][i],'skyblue',
### 					   		Results['all_pft']['dates_win'][i],Results['all_pft']['trend_gain_tsw'][i],'b--',
### 					   		Results['all_pft']['dates_win'][i],[0]*Results['all_pft']['dates_win'][i].size,'k',linewidth = 0.5)
### 			
### 			ax[i].set_title("per : %.1f and loss and gain trends: %.2f and %.2f")# %(per,increase_loss,increase_gain),fontsize = 10)
### 		plt.show()
### 	#ploting the global timeseries when the threshold of the all pft is used to find components 
### 	if plt_global_ts == 'pgts' and pft_comp == 'comp' and pft_idx_in in pft_idx[1:-2]:
### 		fig6		=  plt.figure()
### 		if per_both	== 'b':
### 			plt.title("Global Time Series of pos and neg extremes(threshold : All pft when percentile = %.1f: slope (%.1f & %.1f) kgC/yr  and impact (%.1f & %.1f)%% resp"%(per,Results['comp_pft']['slope_gain']*.365,Results['comp_pft']['slope_loss']*.365,Impact(Results['comp_pft']['trend_gain']),Impact(Results['comp_pft']['trend_loss'])))
### 			plt.plot(Results['all_pft']['dates_win'].flatten(),Results['comp_pft']['global_gpp_loss'], 	'r',
### 					 Results['all_pft']['dates_win'].flatten(),Results['comp_pft']['trend_loss'],		'k--',
### 					 Results['all_pft']['dates_win'].flatten(),Results['comp_pft']['global_gpp_gain'], 	'g',
### 					 Results['all_pft']['dates_win'].flatten(),Results['comp_pft']['trend_gain'],		'k--',
### 					 Results['all_pft']['dates_win'].flatten(),[0]*Results['all_pft']['dates_win'].flatten().size,'k',linewidth = 0.5)
### 			if pft_comp != 'comp':
### 				plt.title("Global Time Serie of pos and neg extremes when percentile = %.1f: slope (%.1f & %.1f) kgC/yr  and impact (%.1f & %.1f)%% resp \n ALL_pft slope (%.1f & %.1f)kgC/yr  and impact (%.1f & %.1f)%% "%(per,Results['comp_pft']['slope_gain']*.365,Results['comp_pft']['slope_loss']*.365,Impact(Results['comp_pft']['trend_gain']),Impact(Results['comp_pft']['trend_loss']),Results['all_pft']['slope_gain']*.365,Results['all_pft']['slope_loss']*.365,Impact(Results['all_pft']['trend_gain']),Impact(Results['all_pft']['trend_loss'])))
### 				plt.plot(Results['all_pft']['dates_win'].flatten(),Results['all_pft']['global_gpp_loss'], 	'skyblue',
### 					 Results['all_pft']['dates_win'].flatten(),Results['all_pft']['trend_loss'],		'b--',
### 					 Results['all_pft']['dates_win'].flatten(),Results['all_pft']['global_gpp_gain'], 	'skyblue',
### 					 Results['all_pft']['dates_win'].flatten(),Results['all_pft']['trend_gain'],		'b--',linewidth = 0.5)
### 
### 				
### 		else:
### 			if per <50:
### 				plt.title("Global Time Series of neg extremes when percentile = %.1f slope : %.1f kgC/yr and impact :%.1f %%" %(per,Results['sel_pft']['slope_loss'], Impact(Results['sel_pft']['trend_loss'])))
### 				plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_loss'], 	'r',
### 						 Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_loss'],		'k--',
### 						 Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = .5)
### 			else:
### 				plt.title("Global Time Series of pos extremes when percentile = %.1f slope : %.1f kgC/yr and impact :%.1f %%"%(per,Results['sel_pft']['slope_gain'],Impact(Results['sel_pft']['trend_gain'])))
### 				plt.plot(Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['global_gpp_gain'], 	'g',                              
### 				 		 Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_gain'],		'k--',                                  
### 						 Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k',linewidth = 0.5) 
### #		fig2.savefig(in_path+'results/global_timeseries_per_%.1f.png' %(per))
### 		plt.show()
### 	
### 	if plt_global_ts == 'pgts' and pft_comp == 'comp' and pft_idx_in == 99:
### 		fig7 = plt.plot()
### 		plt.title ("Cummulative global time series of negative extremes for all pfts")
### 		for key in keys_neg_ext:
### 			plt.plot(Results['all_pft']['dates_win'].flatten(),Results['cumsum_neg_pft'][key],linewidth = 0.6,label = key)
### 			if Results['cumsum_neg_pft'][key].mean() <= 0.9 * np.array(Results['all_pft']['global_gpp_loss']).mean():
### 				break
### 		plt.plot(Results['all_pft']['dates_win'].flatten(), np.array(Results['all_pft']['global_gpp_loss']),'k--', linewidth = 0.6, label = 'all_pfts')
### 		plt.legend()
### 		plt.show()
### 	
### 	if plt_global_ts == 'pgts' and pft_comp == 'comp' and pft_idx_in == 99:
### 		fig8 = plt.plot()
### 		plt.title ("Cummulative global time series of negative extremes for all pfts")
### 		for key in keys_pos_ext:
### 			plt.plot(Results['all_pft']['dates_win'].flatten(),Results['cumsum_pos_pft'][key],linewidth = 0.6,label = key)
### 			if Results['cumsum_pos_pft'][key].mean() >= 0.9 * np.array(Results['all_pft']['global_gpp_gain']).mean():
### 				break
### 		plt.plot(Results['all_pft']['dates_win'].flatten(), np.array(Results['all_pft']['global_gpp_gain']),'k--', linewidth = 0.6,label = 'all_pfts')
### 		plt.legend()
### 		plt.show()
### 	
### 	if plt_global_ts == 'pgts' and pft_comp == 'comp' and per_both == 'b':
### 		fig10 = plt.figure( tight_layout = True, figsize = (9,5), dpi = 300)
### 		for i in range(Results['sel_pft']['dates_win'].shape[0]):
### 			plt.plot(Results['sel_pft']['dates_win'][i],np.array(Results['sel_pft']['global_loss_tsw'][i])/10**15,'r--',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'][i],np.array(Results['sel_pft']['trend_loss_tsw'][i])/10**15,'k--',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.text(Results['sel_pft']['dates_win'][i][0], -.6, "%d" %(Results['sel_pft']['slope_loss_tsw'][i]*.365/10**6),size = 8, color = 'k')
### 			plt.plot(Results['sel_pft']['dates_win'][i],[0]*Results['sel_pft']['dates_win'][i].size,'k',linewidth = 0.5)
### 
### 			plt.plot(Results['sel_pft']['dates_win'][i],np.array(Results['sel_pft']['global_gain_tsw'][i])/10**15,'g--',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.plot(Results['sel_pft']['dates_win'][i],np.array(Results['sel_pft']['trend_gain_tsw'][i])/10**15,'k--',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.text(Results['sel_pft']['dates_win'][i][0], .5, "%d" %(Results['sel_pft']['slope_gain_tsw'][i]*.365/10**6),size = 8, color = 'k')
### 
### 			plt.plot(Results['all_pft']['dates_win'][i],np.array(Results['all_pft']['global_loss_tsw'][i])/10**15,'orangered',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.plot(Results['all_pft']['dates_win'][i],np.array(Results['all_pft']['trend_loss_tsw'][i])/10**15,'k--',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.text(Results['all_pft']['dates_win'][i][0], -.8, "%d" %(Results['all_pft']['slope_loss_tsw'][i]*.365/10**6),size = 8, color = 'k')
### 			plt.plot(Results['all_pft']['dates_win'][i],[0]*Results['all_pft']['dates_win'][i].size,'k',linewidth = 0.5)
### 
### 			plt.plot(Results['all_pft']['dates_win'][i],np.array(Results['all_pft']['global_gain_tsw'][i])/10**15,'lightgreen',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.plot(Results['all_pft']['dates_win'][i],np.array(Results['all_pft']['trend_gain_tsw'][i])/10**15,'k--',label = pft_names[pft_idx_in],linewidth = 0.5)
### 			plt.text(Results['all_pft']['dates_win'][i][0], .7, "%d" %(Results['all_pft']['slope_gain_tsw'][i]*.365/10**6),size = 8, color = 'k')
### 		fig10.savefig('global_time_series_with_changing_slopes')	
### 
### 	"""
### 	#plotting the timeseries for 50 yer window
### 	if pft_comp_win == 'compwin':
### 		fig9,ax	=	plt.subplots (Results['all_pft']['nrows_tsw'],Results['all_pft']['ncols_tsw'])
### 		ax		=	ax.ravel()
### 		plt.suptitle ("Timeseries of %d year slices" %window)
### 		for i in range(Results['sel_pft']['dates_win'].shape[0]):
### 			ax[i].plot(	Results['all_pft']['dates_win'][i],Results['all_pft']['global_loss_tsw'][i],'skyblue',
### 				   		Results['all_pft']['dates_win'][i],Results['all_pft']['trend_loss_tsw'][i],'b--',
### 				   		Results['all_pft']['dates_win'][i],Results['all_pft']['global_gain_tsw'][i],'skyblue',
### 				   		Results['all_pft']['dates_win'][i],Results['all_pft']['trend_gain_tsw'][i],'b--',
### 				   		Results['all_pft']['dates_win'][i],[0]*Results['all_pft']['dates_win'][i].size,'k',linewidth = 0.5)
### 			
### 			ax[i].set_title("per : %.1f and loss and gain trends: %.2f and %.2f")# %(per,increase_loss,increase_gain),fontsize = 10)
### 	"""
### Ploting()
### Tables()
### """
### Ploting the shaded windows: 
### tmp_idx = np.arange(0,5401,300)
### tmp_idx[-1]=tmp_idx[-1]-1
### dates_ticks = []
### for i in tmp_idx:
### 	a = Results['sel_pft']['dates_win'].flatten()[i]
### 	dates_ticks.append(a)
### for i in range(len(tmp_idx)-1):
### 	fig,ax = plt.subplots(tight_layout = True, figsize = (7,2.3), dpi = 400)
### 	ax.plot(Results['sel_pft']['dates_win'].flatten(), Results['sel_pft']['global_gpp_loss']/10**15, 'r--',
### 			Results['sel_pft']['dates_win'].flatten(),Results['sel_pft']['trend_loss']/10**15,        'k--',
### 			Results['sel_pft']['dates_win'].flatten(),[0]*Results['sel_pft']['dates_win'].flatten().size,'k--',linewidth = .5,alpha = .5)
### 	ax.axvspan(Results['sel_pft']['dates_win'].flatten()[tmp_idx[i]],Results['sel_pft']['dates_win'].flatten()[tmp_idx[i+1]],alpha = 0.4,color = 'skyblue')
### 	ax.axvspan(Results['sel_pft']['dates_win'].flatten()[tmp_idx[5]],Results['sel_pft']['dates_win'].flatten()[tmp_idx[6]],alpha = 0.4,color = 'gray')
### 	ax.set_xticks(dates_ticks) 
### 	ax.annotate('Reference Time Period',xy=(Results['sel_pft']['dates_win'].flatten()[1400],-.5),xytext=(Results['sel_pft']['dates_win'].flatten()[1800],-.55),fontsize = 10,arrowprops=dict(facecolor='black', shrink=0.05,))
### 	for tick in ax.get_xticklabels():
### 		tick.set_rotation(90)
### 	plt.xlabel('Time', fontsize = 14) 
### 	plt.ylabel('Pg C', fontsize = 14)
### 	fig.savefig(in_path+'results/win%s_global_timeseries_per_%.1f.png' %(format(i,'02'),per))
### 	fig.savefig(in_path+'results/win%s_global_timeseries_per_%.1f.pdf' %(format(i,'02'),per))
### 
### 
### 
### 
### """
### # iin progres - - - 
