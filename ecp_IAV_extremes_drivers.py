# Bharat Sharma
# python 2.7
# To calculate the IAV of GPP at global and pixel level, the relative IAV i.e. IAV was computed after dividing the anomalies with trend of GPP 
# To calculate the IAV of Climate drivers pixel level 
# to calculate the IAV standard deviation was computed from 1850 to 1874, 1899, 1924 , ... , 2299
import numpy as np              
import pandas as pd
import  matplotlib as mpl
mpl.use('Agg')
from functions import time_dim_dates, index_and_dates_slicing, geo_idx,norm, create_seq_mat
import datetime as dt
import netCDF4 as nc4
import  matplotlib.pyplot as plt
import argparse
from    mpl_toolkits.basemap import Basemap
from    matplotlib import cm
import ssa_copy_2   as ssa
from mpi4py import MPI

parser  = argparse.ArgumentParser()
parser.add_argument ('--lat_idx'	,'-lat_idx'	, help = ' Lat coordinate of the location'							, type = int	, default = 999)
parser.add_argument ('--lon_idx'	,'-lon_idx'	, help = ' Lon coordinate of the location'							, type = int	, default = 999) 
parser.add_argument('--model_config', '-config'	, help = "Model Configuration"   									, type = str	, default = 'wo_lulcc'	) # w_lulcc, wo_lulcc, all
parser.add_argument('--ext_type'    ,   '-ext_typ'  , help = "Type of extreme analysis (pos or neg)"    , type= str, default= 'neg'    )
parser.add_argument('--thres_type'  ,   '-th_typ'   , help = "Type of thresholds (independent or misc)" , type= str, default= 'misc'   )
parser.add_argument('--variables',   '-vars'  , help = "Type all the variables you want to plot( separated by '-' minus sign" , type= str, default= "gpp-prcp-sm-tmax-col_fire_closs")
args = parser.parse_args()

# Running the code
# Enter lat,lon as 999,999 for global scale
# run ecp_IAV_extremes_drivers.py -th_typ misc -vars gpp-prcp-sm-tmax-col_fire_closs -config wo_lulcc -ext_typ neg -lat_idx 102 -lon_idx 26

lat_in	= args.lat_idx
lon_in	= args.lon_idx
conf	= args.model_config
ext_type= args.ext_type
th_type	= args.thres_type
variables_list =  (args.variables).split("-")

""" 
To compute:
-----------
1. to find global IAV of gpp, gpp anomalies and gpp extremes to explain the more threshold of the with LULCC
2. For the amazon discussion compute IAV for per pixel and plot the TS and IAV of the drivers as well
3. not only compute 2 but also compute relative IAV i.e std/mean or something to scale

"""

# Check
a = ['gpp','prcp', 'sm', 'tmax', 'col_fire_closs']
if a == variables_list:
	Variables_name = ['GPP','Prcp', 'Soilmoist',r'$\rm{T}_{\rm{max}}$','Fire']
else:
	Variables_name = variables_list
#reading driver filenames!
nc_data							= {}
nc_data['in']					= {}
nc_data['in' ]['wo_lulcc']		= {}
nc_data['in' ]['w_lulcc']		= {}

nc_data['in']['wo_lulcc']['gpp_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_anomalies_gC.nc'	)['gpp'	]
nc_data['in']['w_lulcc' ]['gpp_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/with_land_use_change/cesm1bgc_with_lulcc_gpp_anomalies_gC.nc')['gpp'	]
nc_data['in']['wo_lulcc']['gpp']		= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_gC.nc'			)['gpp'	]
nc_data['in']['w_lulcc' ]['gpp']		= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/with_land_use_change/cesm1bgc_with_lulcc_gpp_gC.nc'			)['gpp'	]

nc_data['in']['wo_lulcc']['prcp']		= nc4.Dataset('/home/ud4/CESM1-BGC/prcp/without_land_use_change/cesm1bgc_pftcon_prcp.nc'			)['prcp']
nc_data['in']['w_lulcc' ]['prcp']		= nc4.Dataset('/home/ud4/CESM1-BGC/prcp/with_land_use_change/cesm1bgc_with_lulcc_prcp.nc'			)['prcp']
nc_data['in']['wo_lulcc']['sm']			= nc4.Dataset('/home/ud4/CESM1-BGC/sm/without_land_use_change/cesm1bgc_pftcon_sm.nc'				)['sm'	]
nc_data['in']['w_lulcc' ]['sm']			= nc4.Dataset('/home/ud4/CESM1-BGC/sm/with_land_use_change/cesm1bgc_with_lulcc_sm.nc'				)['sm'	]
nc_data['in']['wo_lulcc']['tmax']		= nc4.Dataset('/home/ud4/CESM1-BGC/tmax/without_land_use_change/cesm1bgc_pftcon_tmax.nc'			)['tmax']
nc_data['in']['w_lulcc' ]['tmax']		= nc4.Dataset('/home/ud4/CESM1-BGC/tmax/with_land_use_change/cesm1bgc_with_lulcc_tmax.nc'			)['tmax']
nc_data['in']['wo_lulcc']['tsa']		= nc4.Dataset('/home/ud4/CESM1-BGC/tsa/without_land_use_change/cesm1bgc_pftcon_tsa.nc'				)['tsa'	]
nc_data['in']['w_lulcc' ]['tsa']		= nc4.Dataset('/home/ud4/CESM1-BGC/tsa/with_land_use_change/cesm1bgc_with_lulcc_tsa.nc'				)['tsa'	]

nc_data['in']['wo_lulcc']['prcp_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/prcp/without_land_use_change/cesm1bgc_pftcon_prcp_anomalies.nc'	)['prcp']
nc_data['in']['w_lulcc' ]['prcp_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/prcp/with_land_use_change/cesm1bgc_with_lulcc_prcp_anomalies.nc'	)['prcp']
nc_data['in']['wo_lulcc']['pme_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/pme/without_land_use_change/cesm1bgc_pftcon_pme_anomalies.nc'	)['pme'	]
nc_data['in']['w_lulcc' ]['pme_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/pme/with_land_use_change/cesm1bgc_with_lulcc_pme_anomalies.nc'	)['pme'	]
nc_data['in']['wo_lulcc']['sm_ano']		= nc4.Dataset('/home/ud4/CESM1-BGC/sm/without_land_use_change/cesm1bgc_pftcon_sm_anomalies.nc'		)['sm'	]
nc_data['in']['w_lulcc' ]['sm_ano']		= nc4.Dataset('/home/ud4/CESM1-BGC/sm/with_land_use_change/cesm1bgc_with_lulcc_sm_anomalies.nc'		)['sm'	]
nc_data['in']['wo_lulcc']['tmax_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/tmax/without_land_use_change/cesm1bgc_pftcon_tmax_anomalies.nc'	)['tmax']
nc_data['in']['w_lulcc' ]['tmax_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/tmax/with_land_use_change/cesm1bgc_with_lulcc_tmax_anomalies.nc'	)['tmax']
nc_data['in']['wo_lulcc']['tsa_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/tsa/without_land_use_change/cesm1bgc_pftcon_tsa_anomalies.nc'	)['tsa'	]
nc_data['in']['w_lulcc' ]['tsa_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/tsa/with_land_use_change/cesm1bgc_with_lulcc_tsa_anomalies.nc'	)['tsa'	]
nc_data['in']['wo_lulcc']['tmin_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/tmin/without_land_use_change/cesm1bgc_pftcon_tmin_anomalies.nc'	)['tmin']
nc_data['in']['w_lulcc' ]['tmin_ano']	= nc4.Dataset('/home/ud4/CESM1-BGC/tmin/with_land_use_change/cesm1bgc_with_lulcc_tmin_anomalies.nc'	)['tmin']

iav = {}
iav['global'] = {}
iav['global']['wo_lulcc'] = {}
iav['global']['w_lulcc'] = {}

iav['global']['wo_lulcc']['gpp'] = {}

Win_Start_Years = np.arange(1850,2300,25)  
Time_Windows = [str(1850)+'-'+str(i+24) for i in Win_Start_Years]
lf      	=   nc4.Dataset('/home/ud4/CESM1-BGC/prcp/without_land_use_change/cesm1bgc_pftcon_prcp.nc').variables['landfrac']
area      	=   nc4.Dataset('/home/ud4/CESM1-BGC/prcp/without_land_use_change/cesm1bgc_pftcon_prcp.nc').variables['area']
lt_ln_mat   =   create_seq_mat()
land_1d     =   lt_ln_mat[~lf[...].mask]

# IAV wrt 1850 for every 25 year time period
# ------------------------------------------	i
if (lat_in == 999 and lon_in == 999):
	print " Running the gobal calculations.............."
	# Calculation of the global monthly time series
	gpp_global_mon_wo_lulcc		= np.array([np.sum(nc_data['in']['wo_lulcc']['gpp_ano'][i,:,:]) 	for i in range(nc_data['in']['wo_lulcc']['gpp_ano'].shape[0])])
	gpp_global_mon_w_lulcc		= np.array([np.sum(nc_data['in']['w_lulcc']['gpp_ano'][i,:,:]) 		for i in range(nc_data['in']['w_lulcc']['gpp_ano'].shape[0])])
	
	prcp_global_mon_wo_lulcc	= np.array([np.sum(nc_data['in']['wo_lulcc']['prcp_ano'][i,:,:]) 	for i in range(nc_data['in']['wo_lulcc']['prcp_ano'].shape[0])])
	prcp_global_mon_w_lulcc		= np.array([np.sum(nc_data['in']['w_lulcc']['prcp_ano'][i,:,:]) 	for i in range(nc_data['in']['w_lulcc']['prcp_ano'].shape[0])])

	pme_global_mon_wo_lulcc		= np.array([np.sum(nc_data['in']['wo_lulcc']['pme_ano'][i,:,:]) 	for i in range(nc_data['in']['wo_lulcc']['pme_ano'].shape[0])])
	pme_global_mon_w_lulcc		= np.array([np.sum(nc_data['in']['w_lulcc']['pme_ano'][i,:,:]) 		for i in range(nc_data['in']['w_lulcc']['pme_ano'].shape[0])])

	sm_global_mon_awm_wo_lulcc		= np.array([np.sum(nc_data['in']['wo_lulcc']['sm_ano'][i,:,:]) 		for i in range(nc_data['in']['wo_lulcc']['sm_ano'].shape[0])])
	sm_global_mon_awm_w_lulcc		= np.array([np.sum(nc_data['in']['w_lulcc']['sm_ano'][i,:,:]) 		for i in range(nc_data['in']['w_lulcc']['sm_ano'].shape[0])])

	tmax_global_mon_awm_wo_lulcc	= np.array([np.average(nc_data['in']['wo_lulcc']['tmax_ano'][i,:,:]	, weights = area[...]) 	for i in range(nc_data['in']['wo_lulcc']['tmax_ano'].shape[0])])
	tmax_global_mon_awm_w_lulcc		= np.array([np.average(nc_data['in']['w_lulcc']['tmax_ano'][i,:,:]	, weights = area[...]) 	for i in range(nc_data['in']['w_lulcc']['tmax_ano'].shape[0])])
	
	tsa_global_mon_awm_wo_lulcc		= np.array([np.average(nc_data['in']['wo_lulcc']['tsa_ano'][i,:,:]	, weights = area[...]) 	for i in range(nc_data['in']['wo_lulcc']['tsa_ano'].shape[0])])
	tsa_global_mon_awm_w_lulcc		= np.array([np.average(nc_data['in']['w_lulcc']['tsa_ano'][i,:,:]	, weights = area [...]) for i in range(nc_data['in']['w_lulcc']['tsa_ano'].shape[0])])

	tmin_global_mon_awm_wo_lulcc	= np.array([np.average(nc_data['in']['wo_lulcc']['tmin_ano'][i,:,:]	, weights = area[...]) 	for i in range(nc_data['in']['wo_lulcc']['tmin_ano'].shape[0])])
	tmin_global_mon_awm_w_lulcc		= np.array([np.average(nc_data['in']['w_lulcc']['tmin_ano'][i,:,:]	, weights = area[...]) 	for i in range(nc_data['in']['w_lulcc']['tmin_ano'].shape[0])])

		# IAV is the standard deviation calculated from 1850 up to different time windows of the study
	iav_gpp_global_wo_lulcc	= np.zeros((18, 192, 288))
	iav_gpp_global_w_lulcc		= np.zeros((18, 192, 288))
	iav_prcp_global_wo_lulcc	= np.zeros((18, 192, 288))
	iav_prcp_global_w_lulcc	= np.zeros((18, 192, 288))
	iav_pme_global_wo_lulcc	= np.zeros((18, 192, 288))
	iav_pme_global_w_lulcc		= np.zeros((18, 192, 288))
	iav_tmax_global_wo_lulcc	= np.zeros((18, 192, 288))
	iav_tmax_global_w_lulcc	= np.zeros((18, 192, 288))
	iav_tsa_global_wo_lulcc	= np.zeros((18, 192, 288))
	iav_tsa_global_w_lulcc		= np.zeros((18, 192, 288))
	iav_tmin_global_wo_lulcc	= np.zeros((18, 192, 288))
	iav_tmin_global_w_lulcc	= np.zeros((18, 192, 288))


	iav_gpp_global_wo_lulcc =	[gpp_global_mon_wo_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]
	iav_gpp_global_w_lulcc  = 	[gpp_global_mon_w_lulcc[0:(i+1)*300].std() 		for i in range(len(Win_Start_Years))]

	iav_pme_global_wo_lulcc =	[pme_global_mon_wo_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]
	iav_pme_global_w_lulcc  = 	[pme_global_mon_w_lulcc[0:(i+1)*300].std() 		for i in range(len(Win_Start_Years))]

	iav_prcp_global_wo_lulcc =	[prcp_global_mon_wo_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]
	iav_prcp_global_w_lulcc  = 	[prcp_global_mon_w_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]

	iav_sm_global_wo_lulcc =	[sm_global_mon_awm_wo_lulcc[0:(i+1)*300].std() 		for i in range(len(Win_Start_Years))]
	iav_sm_global_w_lulcc  = 	[sm_global_mon_awm_w_lulcc[0:(i+1)*300].std() 		for i in range(len(Win_Start_Years))]

	iav_tmax_global_wo_lulcc =	[tmax_global_mon_awm_wo_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]
	iav_tmax_global_w_lulcc  = 	[tmax_global_mon_awm_w_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]

	iav_tsa_global_wo_lulcc =	[tsa_global_mon_awm_wo_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]
	iav_tsa_global_w_lulcc  = 	[tsa_global_mon_awm_w_lulcc[0:(i+1)*300].std() 		for i in range(len(Win_Start_Years))]

	iav_tmin_global_wo_lulcc =	[tmin_global_mon_awm_wo_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]
	iav_tmin_global_w_lulcc  = 	[tmin_global_mon_awm_w_lulcc[0:(i+1)*300].std() 	for i in range(len(Win_Start_Years))]

	
	df_global_gpp = pd.DataFrame({'Without LULCC' : iav_gpp_global_wo_lulcc, 'With LULCC':iav_gpp_global_w_lulcc})
	df_global_gpp.index		= Time_Windows
	
	df_global_prcp = pd.DataFrame({'Without LULCC' : iav_prcp_global_wo_lulcc, 'With LULCC':iav_prcp_global_w_lulcc})
	df_global_prcp.index	= Time_Windows

	df_global_pme = pd.DataFrame({'Without LULCC' : iav_pme_global_wo_lulcc, 'With LULCC':iav_pme_global_w_lulcc})
	df_global_pme.index		= Time_Windows

	df_global_sm = pd.DataFrame({'Without LULCC' : iav_sm_global_wo_lulcc, 'With LULCC':iav_sm_global_w_lulcc})
	df_global_sm.index		= Time_Windows

	df_global_tmax = pd.DataFrame({'Without LULCC' : iav_tmax_global_wo_lulcc, 'With LULCC':iav_tmax_global_w_lulcc})
	df_global_tmax.index	= Time_Windows

	df_global_tsa = pd.DataFrame({'Without LULCC' : iav_tsa_global_wo_lulcc, 'With LULCC':iav_tsa_global_w_lulcc})
	df_global_tsa.index		= Time_Windows

	df_global_tmin = pd.DataFrame({'Without LULCC' : iav_tmin_global_wo_lulcc, 'With LULCC':iav_tmin_global_w_lulcc})
	df_global_tmin.index	= Time_Windows

	# Saving the files
	# ----------------
	df_global_gpp	.to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/iav_gpp_global.csv', sep = ',')
	df_global_prcp	.to_csv('/home/ud4/CESM1-BGC/prcp/comparison/pftcon_lulcc/iav_prcp_global.csv'	, sep = ',')
	df_global_pme	.to_csv('/home/ud4/CESM1-BGC/pme/comparison/pftcon_lulcc/iav_pme_global.csv'	, sep = ',')
	df_global_sm	.to_csv('/home/ud4/CESM1-BGC/sm/comparison/pftcon_lulcc/iav_sm_global.csv'		, sep = ',')
	df_global_tmax	.to_csv('/home/ud4/CESM1-BGC/tmax/comparison/pftcon_lulcc/iav_tmax_global.csv'	, sep = ',')
	df_global_tsa	.to_csv('/home/ud4/CESM1-BGC/tsa/comparison/pftcon_lulcc/iav_tsa_global.csv'	, sep = ',')
	df_global_tmin	.to_csv('/home/ud4/CESM1-BGC/tmin/comparison/pftcon_lulcc/iav_tmin_global.csv'	, sep = ',')

	
	print aaa
else:
	# IAV of GPP for a pixel
	# ----------------------
	gpp_px_mon_wo_lulcc = nc_data['in']['wo_lulcc']['gpp_ano'][:,lat_in,lon_in]
	gpp_px_mon_w_lulcc  = nc_data['in']['w_lulcc' ]['gpp_ano'][:,lat_in,lon_in]
	
	iav_gpp_px_wo_lulcc =	[gpp_px_mon_wo_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]
	iav_gpp_px_w_lulcc  = 	[gpp_px_mon_w_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]
	
	df_pixel = pd.DataFrame({'Without LULCC' : iav_gpp_px_wo_lulcc, 'With LULCC':iav_gpp_px_w_lulcc})
	df_pixel.index	= Time_Windows

	df_pixel_norm = pd.DataFrame({'Without LULCC' : norm(iav_gpp_px_wo_lulcc), 'With LULCC': norm(iav_gpp_px_w_lulcc)})
	df_pixel_norm.index	= Time_Windows

	# Rel IAV of GPP for a pixel: IAV devided by trend in GPP
	# -------------------------------------------------------
	trend_px_wo_lulcc = ssa.SSA (nc_data['in']['wo_lulcc']['gpp'][:,lat_in,lon_in],120,[])[0]
	trend_px_w_lulcc = ssa.SSA (nc_data['in']['w_lulcc']['gpp'][:,lat_in,lon_in],120,[])[0]
	
	mean25yr_gpp_trend_wo_lulcc	= [trend_px_wo_lulcc [i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]
	mean25yr_gpp_trend_w_lulcc	= [trend_px_w_lulcc [i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]

	df_pixel_trend_gpp 			= pd.DataFrame({'Without LULCC' : mean25yr_gpp_trend_wo_lulcc, 'With LULCC': mean25yr_gpp_trend_w_lulcc})
	df_pixel_trend_gpp.index	= Time_Windows

	# Relative GPP anomalies to trend
	# -------------------------------
	rel_gpp_ano_wo_lulcc 	= nc_data['in']['wo_lulcc']['gpp_ano'][:,lat_in,lon_in] / trend_px_wo_lulcc
	rel_gpp_ano_w_lulcc  	= nc_data['in']['w_lulcc' ]['gpp_ano'][:,lat_in,lon_in] / trend_px_w_lulcc
	
	# IAV of relative GPP 
	rel_iav_gpp_wo_lulcc	= [rel_gpp_ano_wo_lulcc[0:(i+1)*300].std() for i in range (len(Win_Start_Years))]
	rel_iav_gpp_w_lulcc		= [rel_gpp_ano_w_lulcc[0:(i+1)*300].std()  for i in range (len(Win_Start_Years))]
	
	df_pixel_rel	= pd.DataFrame({'Without LULCC' : rel_iav_gpp_wo_lulcc, 'With LULCC': rel_iav_gpp_w_lulcc})
	df_pixel_rel.index = Time_Windows

	df_pixel_rel_norm		= pd.DataFrame({'Without LULCC' : norm(rel_iav_gpp_wo_lulcc), 'With LULCC': norm(rel_iav_gpp_w_lulcc)})
	df_pixel_rel_norm.index = Time_Windows

	# Saving the files
	# ---------------
	df_pixel_trend_gpp	.to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/trend_gpp_pixel_lat_%d_lon_%d.csv'%(lat_in,lon_in)	, sep = ',')

	df_pixel	.to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/iav_gpp_pixel_lat_%d_lon_%d.csv'%(lat_in,lon_in)			, sep = ',')
	df_pixel_rel.to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/iav_gpp_pixel_rel_trend_lat_%d_lon_%d.csv'%(lat_in,lon_in)	, sep = ',')

	df_pixel_norm	 .to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/norm_iav_gpp_pixel_lat_%d_lon_%d.csv'%(lat_in,lon_in)			, sep = ',')
	df_pixel_rel_norm.to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/norm_iav_gpp_pixel_rel_trend_lat_%d_lon_%d.csv'%(lat_in,lon_in)	, sep = ',')

	# IAV of some of the climate Drivers
	# ----------------------------------
	prcp_px_mon_wo_lulcc = nc_data['in']['wo_lulcc']['prcp_ano'][:,lat_in,lon_in]
	prcp_px_mon_w_lulcc  = nc_data['in']['w_lulcc' ]['prcp_ano'][:,lat_in,lon_in]
	
	iav_prcp_px_wo_lulcc =	[prcp_px_mon_wo_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]
	iav_prcp_px_w_lulcc  = 	[prcp_px_mon_w_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]

	sm_px_mon_wo_lulcc = nc_data['in']['wo_lulcc']['sm_ano'][:,lat_in,lon_in]
	sm_px_mon_w_lulcc  = nc_data['in']['w_lulcc' ]['sm_ano'][:,lat_in,lon_in]
	
	iav_sm_px_wo_lulcc =	[sm_px_mon_wo_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]
	iav_sm_px_w_lulcc  = 	[sm_px_mon_w_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]

	tmax_px_mon_wo_lulcc = nc_data['in']['wo_lulcc']['tmax_ano'][:,lat_in,lon_in]
	tmax_px_mon_w_lulcc  = nc_data['in']['w_lulcc' ]['tmax_ano'][:,lat_in,lon_in]
	
	iav_tmax_px_wo_lulcc =	[tmax_px_mon_wo_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]
	iav_tmax_px_w_lulcc  = 	[tmax_px_mon_w_lulcc[0:(i+1)*300].std() for i in range(len(Win_Start_Years))]

	df_pixel_drivers = pd.DataFrame({'Without LULCC (prcp)' : iav_prcp_px_wo_lulcc, 'With LULCC (prcp)':iav_prcp_px_w_lulcc, 'Without LULCC (sm)' : iav_sm_px_wo_lulcc, 'With LULCC (sm)':iav_sm_px_w_lulcc, 'Without LULCC (tmax)' : iav_tmax_px_wo_lulcc, 'With LULCC (tmax)':iav_tmax_px_w_lulcc})
	df_pixel_drivers.index	= Time_Windows

	#Saving drivers_IAV
	df_pixel_drivers .to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/iav_pixel_drivers_lat_%d_lon_%d.csv'%(lat_in,lon_in) , sep = ',')

	# Normalizing
	# ----------
	df_pixel_drivers_norm = pd.DataFrame({'Without LULCC (prcp)' : norm(iav_prcp_px_wo_lulcc), 'With LULCC (prcp)': norm(iav_prcp_px_w_lulcc), 'Without LULCC (sm)' : norm(iav_sm_px_wo_lulcc), 'With LULCC (sm)': norm(iav_sm_px_w_lulcc), 'Without LULCC (tmax)' : norm(iav_tmax_px_wo_lulcc), 'With LULCC (tmax)': norm(iav_tmax_px_w_lulcc)})
	df_pixel_drivers_norm.index	= Time_Windows

	#Saving drivers_IAV
	df_pixel_drivers_norm .to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/norm_iav_pixel_drivers_lat_%d_lon_%d.csv'%(lat_in,lon_in) , sep = ',')

	# trend of climate drivers
	# -------------------------------------------------------
	prcp_trend_px_wo_lulcc 		= ssa.SSA (nc_data['in']['wo_lulcc']['prcp'][:,lat_in,lon_in],120,[]) [0]
	prcp_trend_px_w_lulcc 		= ssa.SSA (nc_data['in']['w_lulcc']	['prcp'][:,lat_in,lon_in],120,[]) [0]

	sm_trend_px_wo_lulcc 		= ssa.SSA (nc_data['in']['wo_lulcc']['sm'][:,lat_in,lon_in],120,[]) [0]
	sm_trend_px_w_lulcc 		= ssa.SSA (nc_data['in']['w_lulcc']	['sm'][:,lat_in,lon_in],120,[]) [0]

	tmax_trend_px_wo_lulcc 		= ssa.SSA (nc_data['in']['wo_lulcc']['tmax'][:,lat_in,lon_in],120,[]) [0]
	tmax_trend_px_w_lulcc 		= ssa.SSA (nc_data['in']['w_lulcc']	['tmax'][:,lat_in,lon_in],120,[]) [0]

	mean25yr_prcp_trend_px_wo_lulcc =	[prcp_trend_px_wo_lulcc[i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]
	mean25yr_prcp_trend_px_w_lulcc  =	[prcp_trend_px_w_lulcc[i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]

	mean25yr_sm_trend_px_wo_lulcc =	[sm_trend_px_wo_lulcc[i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]
	mean25yr_sm_trend_px_w_lulcc  =	[sm_trend_px_w_lulcc[i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]

	mean25yr_tmax_trend_px_wo_lulcc =	[tmax_trend_px_wo_lulcc[i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]
	mean25yr_tmax_trend_px_w_lulcc  =	[tmax_trend_px_w_lulcc[i:(i+1)*300].mean() for i in range(len(Win_Start_Years))]


	df_pixel_drivers_trend = pd.DataFrame({'Without LULCC (prcp)' : mean25yr_prcp_trend_px_wo_lulcc, 'With LULCC (prcp)':mean25yr_prcp_trend_px_w_lulcc, 'Without LULCC (sm)' : mean25yr_sm_trend_px_wo_lulcc, 'With LULCC (sm)':mean25yr_sm_trend_px_w_lulcc, 'Without LULCC (tmax)' : mean25yr_tmax_trend_px_wo_lulcc, 'With LULCC (tmax)':mean25yr_tmax_trend_px_w_lulcc})
	df_pixel_drivers_trend.index	= Time_Windows
	#Saving drivers_trend
	df_pixel_drivers_trend .to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/pixel_trend_drivers_lat_%d_lon_%d.csv'%(lat_in,lon_in) , sep = ',')

	# Normalizing:
	# -----------
	df_pixel_drivers_trend_norm = pd.DataFrame({'Without LULCC (prcp)' : norm(mean25yr_prcp_trend_px_wo_lulcc), 'With LULCC (prcp)':norm(mean25yr_prcp_trend_px_w_lulcc), 'Without LULCC (sm)' : norm(mean25yr_sm_trend_px_wo_lulcc), 'With LULCC (sm)': norm(mean25yr_sm_trend_px_w_lulcc), 'Without LULCC (tmax)' : norm(mean25yr_tmax_trend_px_wo_lulcc), 'With LULCC (tmax)': norm(mean25yr_tmax_trend_px_w_lulcc)})
	df_pixel_drivers_trend_norm.index = Time_Windows
	#Saving drivers_trend_norm
	df_pixel_drivers_trend_norm .to_csv('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/IAV/pixel_norm_trend_drivers_lat_%d_lon_%d.csv'%(lat_in,lon_in) , sep = ',')





