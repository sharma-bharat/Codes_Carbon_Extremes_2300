# Bharat Sharma
# python 2.7
# Calculations of correlation of triggers

""" the aim of this code is to identify extremes and find correlation of the corresponding drivers (with lags)
   and doing the analysis on consecutiuve time windows
   o	Data for Extremes is a list of triggers of GPP anomalies. Data for Drivers is the anomalies of TS of all drivers except fire. The cumulative lag effects of drivers are considered by averaging the drivers from lag =1 to nth month.
   o	First, the TS of gpp anomalies and driver anomalies are normalized from 0 to 1.
   o	The patch finder function returns a list of triggers of TCE( i.e. gpp ano and dri ano). e.g., in 25 years or 300 months were have 6 TCEs, we will have a list of 6 gpp ano values and 6 driver ano averaged values.
   o	The correlation is done on these qualified values and cc, pv and slope are returned.
   o	This procedure is done for every pixel, time win of 25 years and 0 - 13 months lags.
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

#Inputs:
#------
parser  = argparse.ArgumentParser()
parser.add_argument('--driver_ano'	,	'-dri_a'	, help = "Driver anomalies" 					, type= str		, default= 'n'		) #tmin
parser.add_argument('--variable'   	, 	'-var'    	, help = "Anomalies of carbon cycle variable" 	, type= str		, default= 'gpp'	)
parser.add_argument('--var_in_gC'	, 	'-var_gC'	, help = "Unit:If the variable is gC" 			, type= str		, default= 'gC'		) # gC, ori
parser.add_argument('--model_config', 	'-config' 	, help = "Model Configuration"   				, type= str		, default= 'w_lulcc') # w_lulcc, wo_lulcc, all
parser.add_argument('--thres_type'	, 	'-th_typ' 	, help = "Type of thresholds (independent or misc)", type= str	, default= 'misc'	) # this is mainly for the misc case study filter
parser.add_argument('--ext_type'	, 	'-ext_typ' 	, help = "Type of extreme analysis (pos or neg)", type= str		, default= 'neg'	) #pos/neg

args = parser.parse_args()

#print args
variable	= args.variable 	# variable original
dri_type	= args.driver_ano
th_type		= args.thres_type
ext_type	= args.ext_type

print args
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
if args.var_in_gC == "gC":
	nc_data['var']['wo_lulcc']['gC' ]	= nc4.Dataset(paths['in']['var']['wo_lulcc']	+'cesm1bgc_pftcon_'		+variable+'_gC.nc'	)
	nc_data['var']['w_lulcc' ]['gC' ]	= nc4.Dataset(paths['in']['var']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+variable+'_gC.nc'	)

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
		area		= nc_data['var'][conf][unit].variables['area']
		ano			= nc_data['var'][conf][unit].variables[variable]
		#time_bounds	= nc_data['var'][conf][unit].variables['time_bounds']
		break
	break

##############################################
window 		= 25 #years
win_len     = 12 * window            #number of months in window years
nwin		= int(time.size/win_len) #number of windows

wins    = np.array([format(i,'02' ) for i in range(nwin)])
#number of months you want to use for lags : 0(during) to 12months)
max_lag		=	12 #months
lags    	= np.array([format(i,'02' ) for i in range(max_lag +1)])
dates_ar    = time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates	= [dates_ar[i*win_len] for i in range(nwin)]#list of start dates of 25 year window
end_dates   = [dates_ar[i*win_len+win_len -1] for i in range(nwin)]#list of end dates of the 25 year window

# Calculation of anomalies
idx_dates_win= []   #indicies of time in 30yr windows
dates_win   = []    #sel dates from time variables in win_len windows


#reading driver filenames!
paths['in' ]['dri']				= {}
paths['in' ]['dri']['wo_lulcc'] = '/home/ud4/CESM1-BGC/'+dri_type+'/without_land_use_change/'
paths['in' ]['dri']['w_lulcc' ] = '/home/ud4/CESM1-BGC/'+dri_type+'/with_land_use_change/'

nc_data['dri']					= {}
nc_data['dri']['wo_lulcc']		= nc4.Dataset(paths['in']['dri']['wo_lulcc']	+'cesm1bgc_pftcon_'		+dri_type+'_anomalies.nc'	)
nc_data['dri']['w_lulcc' ]		= nc4.Dataset(paths['in']['dri']['w_lulcc' ]	+'cesm1bgc_with_lulcc_'	+dri_type+'_anomalies.nc'	)
if dri_type == 'fclosscol'		: nc4.Dataset(paths['in']['dri']['wo_lulcc'] +'cesm1bgc_pftcon_'     +dri_type+'_gC.nc'   		) 
if dri_type == 'col_fire_closs'	: nc4.Dataset(paths['in']['dri']['wo_lulcc'] +'cesm1bgc_pftcon_'     +dri_type+'_gC.nc'   		) 
if dri_type == 'spi_gamma_06' 	: nc4.Dataset(paths['in']['dri']['wo_lulcc'] +'cesm1bgc_pftcon_spi_gamma_06.nc'         		)

# Misc case (th_type_ for this code is basically based on the 'keys_wrt_wo_lulcc_pos' 
keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']

# the casename will be added to the filename in order to save files for different cases separately
paths['out']					= {}
paths['out']['dri']				= {}
if th_type == 'misc':
	# the casename will be added to the filename in order to save files for different cases separately
	case_name 	= ext_type+'_'+conf+'_based_pos_wo_lulcc'
	bin_ext 	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_ext_'+case_name+'.nc').variables['gpp_bin_ext']  # the binary matrix of extremes is calculated for each case in MISC category.
	paths['out']['dri']['wo_lulcc'] = '/home/ud4/CESM1-BGC/'+dri_type+'/comparison/pftcon_lulcc/misc/corr_trigger/' + case_name + '/'
	paths['out']['dri']['w_lulcc' ] = '/home/ud4/CESM1-BGC/'+dri_type+'/comparison/pftcon_lulcc/misc/corr_trigger/' + case_name + '/'
else:
	case_name = ''
	paths['out']['dri']['wo_lulcc'] = '/home/ud4/CESM1-BGC/'+dri_type+'/without_land_use_change/corr_triggers/'
	paths['out']['dri']['w_lulcc' ] = '/home/ud4/CESM1-BGC/'+dri_type+'/with_land_use_change/corr_triggers/'


nc_dri 	=   nc_data['dri'][conf].variables[dri_type]  #driver
nc_area =   nc4.Dataset('/home/ud4/CESM1-BGC/prcp/without_land_use_change/cesm1bgc_pftcon_prcp.nc').variables['area']
lf 		=   nc4.Dataset('/home/ud4/CESM1-BGC/prcp/without_land_use_change/cesm1bgc_pftcon_prcp.nc').variables['landfrac']

#for parallel programming: identifying only the non masked locations for better load balancing
lt_ln_mat 	=	create_seq_mat()
land_1d		= 	lt_ln_mat[~lf[...].mask]

# Parallel Computing task distribution
comm		=	MPI.COMM_WORLD
size		=	comm.Get_size() #Number of processors I am asking for
local_rank	=	comm.Get_rank() #Rank of the current processor


if local_rank == 0:
	corr_coeff_mat 	= np.zeros((nwin,len(lags),lat.size,lon.size))
	p_values_mat	= np.zeros((nwin,len(lags),lat.size,lon.size))
	slope_mat		= np.zeros((nwin,len(lags),lat.size,lon.size))
	
	corr_coeff_temp	= np.zeros((nwin,len(lags),lat.size,lon.size))	
	p_values_temp  	= np.zeros((nwin,len(lags),lat.size,lon.size))	
	slope_temp  	= np.zeros((nwin,len(lags),lat.size,lon.size))	

	
	for i in range(1,size):
		print "receiving from %d slave" %(i)
		#try:
		comm.Recv(corr_coeff_temp, 	source=i, tag=00)
		comm.Recv(p_values_temp, 	source=i, tag=01)
		comm.Recv(slope_temp, 		source=i, tag=02)
		print "from slave %d ..."%i, np.nansum(corr_coeff_temp)
		corr_coeff_mat	= corr_coeff_mat 	+ corr_coeff_temp
		p_values_mat	= p_values_mat 		+ p_values_temp
		slope_mat		= slope_mat 		+ slope_temp
		print "................................. received from %d slave" %(i)

	
	#np.savetxt('temp/'+'cc_rank_%s_.csv'%(format(i,'03')), (0,1), delimiter = ',')

	land_masking_ar 	 = np.ma.masked_all((nwin,len(lags),lat.size,lon.size))
	land_masking_ar[:,:] = lf[...]
	
	corr_coeff_mat 	= np.ma.masked_array(corr_coeff_mat	, mask = land_masking_ar.mask)
	p_values_mat 	= np.ma.masked_array(p_values_mat	, mask = land_masking_ar.mask)
	slope_mat 		= np.ma.masked_array(slope_mat		, mask = land_masking_ar.mask)
	print ".............++++++++++++++++++++test"
	#Saving as nc file
	with nc4.Dataset(paths['out']['dri'][conf] + 'corr_coff_cumlag.nc', mode="w") as dset:
		dset        .createDimension( "win" ,size = nwin)
		dset        .createDimension( "lag" ,size = lags.size)
		dset        .createDimension( "lat"	, size = lat.size)
		dset        .createDimension( "lon"	, size = lon.size)
		w   =   dset.createVariable(varname = "win"  ,datatype = float, dimensions = ("win") ,fill_value = np.nan)                       
		t   =   dset.createVariable(varname = "lag"  ,datatype = float, dimensions = ("lag") ,fill_value = np.nan)
		y   =   dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat") ,fill_value = np.nan)
		x   =   dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon") ,fill_value = np.nan)
		z   =   dset.createVariable(varname = "corr_coeff" 	,datatype = float, dimensions = ("win","lag","lat","lon"),fill_value=np.nan)
		zz  =  dset.createVariable(varname = "p_value"     	,datatype = float, dimensions = ("win","lag","lat","lon"),fill_value=np.nan)
		zzz =  dset.createVariable(varname = "slope"     	,datatype = float, dimensions = ("win","lag","lat","lon"),fill_value=np.nan)
		w.axis  =   "T"
		t.axis  =   "T"
		x.axis  =   "X"
		y.axis  =   "Y"
		w[...]  =   wins
		t[...]  =   lags
		x[...]  =   lon[...]
		y[...]  =   lat[...]
		z[...]  =   corr_coeff_mat
		zz[...] =   p_values_mat
		zzz[...]=   slope_mat
		z.missing_value =   np.nan
		z.stardard_name =   "Correlation coefficient of %s_anomalies with gpp negative extreme triggers"%dri_type
		z.setncattr         ("cell_methods",'Stats.linregress correlation coefficient')
		z.units         =   "no units"
		zz.missing_value =   np.nan
		zz.stardard_name =   "p_value for %s anomalies with gpp negative extreme triggers "%dri_type
		zz.setncattr        ("cell_methods",'p-value for testing correlation')
		zz.units        =   "no units"
		zzz.missing_value =   np.nan
		zzz.stardard_name =   "slope for %s anomalies with gpp negative extreme triggers "%(dri_type)
		zzz.setncattr        ("cell_methods",'slope for attribution correlation')
		zzz.units        =   "no units"
		x.units         =   lat.units
		x.missing_value =   np.nan
		x.setncattr         ("long_name",lat.long_name)
		y.units         =   lon.units
		y.missing_value =   np.nan
		y.setncattr         ("long_name",lon.long_name)
		t.units         =   "lag by months"
		w.units         =   "25 year wins"
		w.setncattr         ("long_name","starting from 1850-01-15 i.e. win=0: 1850-01-15 to 1874-12-15")
	
	
	#Saving the table of gobal average and median of corelation coeff
	cc_table_av	= np.zeros((len(lags),nwin))
	pv_table_av	= np.zeros((len(lags),nwin))
	cc_table_md	= np.zeros((len(lags),nwin))
	pv_table_md	= np.zeros((len(lags),nwin))
	
	#Saving table only for CONUS
	top     = 49.345    +.5
	bottom  = 24.743    -.5
	left    = 360-124.78-.5
	right   = 360-66.95 +.5
	cc_table_conus_av	= np.zeros((len(lags),nwin))
	pv_table_conus_av	= np.zeros((len(lags),nwin))
	
#	np.savetxt('temp/'+'table_start_.csv', (0,1), delimiter = ',')
	
	
	for w in range(nwin):
		for l in range(len(lags)):
			cc_table_av[l,w] 		= np.nanmean(corr_coeff_mat[w,l])
			pv_table_av[l,w] 		= np.nanmean(p_values_mat[w,l])
			cc_table_md[l,w] 		= np.nanmedian(corr_coeff_mat[w,l])
			pv_table_md[l,w] 		= np.nanmedian(p_values_mat[w,l])
			cc_table_conus_av[l,w] 	= np.nanmean(corr_coeff_mat[w,l,geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])])
			pv_table_conus_av[l,w] 	= np.nanmean(p_values_mat[w,l, geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])])

	#np.savetxt('temp/'+'table_end_.csv', (0,1), delimiter = ',')
	
	df_cc_av 		=  pd.DataFrame(data = cc_table_av)
	df_pv_av 		=  pd.DataFrame(data = pv_table_av)
	df_cc_md 		=  pd.DataFrame(data = cc_table_md)
	df_pv_md 		=  pd.DataFrame(data = pv_table_md)
	df_cc_conus_av 	=  pd.DataFrame(data = cc_table_conus_av)
	df_pv_conus_av 	=  pd.DataFrame(data = pv_table_conus_av)

	in_yr   = 1850
	win_yr  = [str(in_yr+i*25) + '-'+str(in_yr +(i+1)*25-1)[2:] for i in np.array(wins,dtype =int)]

	col_names 				= win_yr#['win'+ format(i+1,'02') for i in range(nwin)]
	df_cc_av.columns 		= col_names
	df_pv_av.columns 		= col_names
	df_cc_md.columns 		= col_names
	df_pv_md.columns 		= col_names
	df_cc_conus_av.columns 	= col_names
	df_pv_conus_av.columns 	= col_names

	id_names 	= ['lag'+ format(i,'02') for i in range(len(lags))]
	df_cc_av.index 			= id_names
	df_pv_av.index 			= id_names
	df_cc_md.index 			= id_names
	df_pv_md.index 			= id_names
	df_cc_conus_av.index 	= id_names
	df_pv_conus_av.index 	= id_names

	df_cc_av.to_csv 		(paths['out']['dri'][conf] + 'cc_table_av_cumlag.csv'		)
	df_pv_av.to_csv 		(paths['out']['dri'][conf] + 'pv_table_av_cumlag.csv'		)
	df_cc_md.to_csv 		(paths['out']['dri'][conf] + 'cc_table_md_cumlag.csv'		)
	df_pv_md.to_csv 		(paths['out']['dri'][conf] + 'pv_table_md_cumlag.csv'		)
	df_cc_conus_av.to_csv 	(paths['out']['dri'][conf] + 'cc_table_conus_av_cumlag.csv'	)
	df_pv_conus_av.to_csv 	(paths['out']['dri'][conf] + 'pv_table_conus_av_cumlag.csv'	)

else:
	#chunk size or delta
	local_n		=	int(np.ceil(len(land_1d)/(size-1))) #size -1 because we are sending information to the root

	#calculating the range for every parallel process
	begin_idx	=	(local_rank-1)*local_n #local_rank-1 because the rank starts from 1 and we will miss the index 0 if we do not do this ammendment
	end_idx		=	begin_idx+local_n
	
	#Defining empty arrays
	corr_coeff_ar	= 	np.zeros((nwin,len(lags),lat.size,lon.size))
	p_values_ar		= 	np.zeros((nwin,len(lags),lat.size,lon.size))
	slope_ar		= 	np.zeros((nwin,len(lags),lat.size,lon.size))

	if begin_idx <=	len(land_1d)-1:
		if local_n + begin_idx >len(land_1d)-1:
			end_idx = len(land_1d)-1


		loc_pp		=	land_1d[begin_idx:end_idx]		#locations per processor

		#for multiple windows:
		dates_ar	=	time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
		win_len		=	25 *12 #for 25 years
		start_dates	=	[dates_ar[i*win_len] for i in range(nwin)] #list of the start dates of all 25 year windows
		end_dates	=	[dates_ar[i*win_len+win_len-1] for i in range(nwin)] #list of the end dates of all 25 year windows

		bin_ext_loc = np.zeros((win_len,lat.size,lon.size)) #3d array to capture the True binaray anomalies w.r.t. gpp loss events
 
		for win in range(nwin):
			# idx with store all the time indexs b/w start and end dates from the input array
			idx_loc,_	= index_and_dates_slicing(dates_ar,start_dates[win],end_dates[win])
			ano_loc 	= ano[idx_loc[0]:idx_loc[-1]+1,:,:]
			dri_loc		= nc_dri[idx_loc[0]:idx_loc[-1]+1,:,:]
			bin_ext_loc = bin_ext[idx_loc[0]:idx_loc[-1]+1,:,:]
			bin_ext_loc_mask    = np.ma.masked_equal(bin_ext_loc,0)

			for lag in np.asarray(lags, dtype = int):
				#print "calculating on rank_%d_for %s for window %d and lag %d......"%(local_rank,dri_type,win,lag)
				for pixel in loc_pp:
					lt,ln = np.argwhere(lt_ln_mat == pixel)[0]
					if collections.Counter(ano_loc[:,lt,ln])[0] == len(ano_loc[:,lt,ln]):   #if all values are zeros
						corr_coeff_ar 	[win,lag,lt,ln]	= np.nan
						p_values_ar 	[win,lag,lt,ln]	= np.nan
						slope_ar 		[win,lag,lt,ln]	= np.nan
					elif collections.Counter(bin_ext_loc_mask[:,lt,ln].mask)[True]+ lag >= (bin_ext_loc_mask[:,lt,ln].mask).size:
						corr_coeff_ar 	[win,lag,lt,ln]	= np.nan
						p_values_ar 	[win,lag,lt,ln]	= np.nan
						slope_ar 		[win,lag,lt,ln]	= np.nan
					else:
						#np.savetxt('temp/'+'error_rank_win%s_lag%s_lt%s_ln%s.csv'%(str(win),str(lag),str(lt),str(ln)), (0,1), delimiter = ',')
#						print "error_rank_win%s_lag%s_lt%s_ln%s"%(str(win), str(lag), str(lt), str(ln)) 
						triggers = patch_with_gaps_and_eventsize(bin_ar = bin_ext_loc_mask[:,lt,ln], max_gap =2,min_cont_event_size =3, lag=lag)[1]
						if len(triggers) <=2:
							corr_coeff_ar 	[win,lag,lt,ln]	= np.nan
							p_values_ar 	[win,lag,lt,ln]	= np.nan
							slope_ar 		[win,lag,lt,ln]	= np.nan
						
						else:
							print "error_1rank_win%s_lag%s_lt%s_ln%s"%(str(win), str(lag), str(lt), str(ln)) 
							ano_norm	= norm(ano_loc[:,lt,ln])
							dri_norm	= norm(dri_loc[:,lt,ln])
							ano_trigs	= np.array([ano_norm[arg] for arg in triggers])

							if lag == 0:
								dri_trigs 	= np.array([dri_norm[arg] for arg in triggers])
							else:
								dri_lag 	= cumsum_lagged (dri_norm, lag, ignore_t0= True)
								dri_trigs	= np.array([dri_lag[arg] for arg in triggers])
						
							#Calculations of the correlation coffecients, Pvalue and slope of regression line
							print "error_2rank_win%s_lag%s_lt%s_ln%s"%(str(win), str(lag), str(lt), str(ln)) 
							slope,c,corr_coeff,p_value,stderr	= stats.linregress(dri_trigs,ano_trigs)
							if lt == 89 and ln ==228:
						   		print "error_rank_win%s_lag%s_lt%s_ln%s"%(str(win), str(lag), str(lt), str(ln)) 
								print "bharat...............linregress",slope,c,corr_coeff,p_value,stderr
							corr_coeff_ar 	[win,lag,lt,ln]		= corr_coeff
							p_values_ar 	[win,lag,lt,ln]		= p_value
							slope_ar 		[win,lag,lt,ln]		= slope
							if lt == 89 and ln ==228: 
								print "bharat.......array",corr_coeff_ar[:,:,lt,ln]
	
		print "Sending the information from rank : %d" %local_rank
		print "sending cc from rank %d ... "%local_rank ,  np.nansum(corr_coeff_ar)
		comm.Send(corr_coeff_ar	, dest = 0, tag = 00)
		comm.Send(p_values_ar	, dest = 0, tag = 01)
		comm.Send(slope_ar		, dest = 0, tag = 02)
		#df_temp = pd.DataFrame(data = corr_coeff_ar[win,lag])
		#df_temp.to_csv('temp/cc_%d.csv'%local_rank, sep=',')

		print "Data Sent from rank : %d" %local_rank
	
	else:
		print "Sending blank information from rank : %d" %local_rank
		comm.Send(corr_coeff_ar	, dest = 0, tag = 00)
		comm.Send(p_values_ar	, dest = 0, tag = 01)	
		comm.Send(slope_ar		, dest = 0, tag = 02)
		#np.savetxt(in_path_dri+'corr/temp2/'+'cc_win_%s_lag_%s_pixel_%s_cc_%s_pv_%s.csv'%(format(win,'02'),format(lag,'02'),format(pixel,'05'),corr_coeff, p_values), (corr_coeff, p_values), delimiter = ',')
			
