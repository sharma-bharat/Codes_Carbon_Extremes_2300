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
mpl.use('Agg')

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
#Running: run ecp_dominant_climate_driver_correlation_graphs.py -config wo_lulcc -th_typ misc -ext_typ pos
th_type     = args.thres_type
ext_type    = args.ext_type
conf		= args.model_config

print args
# Misc case (th_type_ for this code is basically based on the 'keys_wrt_wo_lulcc_pos' 
keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
case_name   = ext_type+'_'+conf+'_based_pos_wo_lulcc'

driver_consider = 7

drivers						 = np.array(['prcp','sm','pme','tmax','tsa','tmin','col_fire_closs','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider] 
#drivers						 = np.array(['prcp','sm','pme','tmax','tsa','tmin','fclosscol','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider] 
drivers_names				 = np.array(['Prcp','Soilmoist',r'$P-ET$',r'$T_{max}$',r'$T_{sa}$',r'$T_{min}$','Fire','SPI_gamma_06','SPI_pearson_06','SPEI_gamma_06','SPEI_gamma_12','SPEI_pearson_06'])[:driver_consider]
drivers_code				 = np.array([  10,   20,   30,   40,       50,		60,			70,			80,		90,			100, 110,120])[:driver_consider]
features 			 		 = {}
#features['abr']				 = np.array(['prcp','sm','pme','tmax','tsa', 'tmin','fclosscol','gpp','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider]
features['abr']				 = np.array(['prcp','sm','pme','tmax','tsa', 'tmin','col_fire_closs','gpp','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider]
# The name with which the variables are stored in the nc files:


features['Names']			 	= {}
features['Names']['prcp']	 	= 'prcp'
features['Names']['sm']	 	 	= 'sm'
features['Names']['tmax']	 	= 'tmax'
features['Names']['tmin']	 	= 'tmin'
features['Names']['pme']	 	= 'pme'
#features['Names']['fclosscol']	= 'COL_FIRE_CLOSS'
features['Names']['col_fire_closs']	= 'col_fire_closs'
features['Names']['gpp']	 	= 'gpp'
features['Names']['spi_g06'] 	= 'spi_gamma_06'

paths						= {}
paths['in' ]				= {}
paths['in' ]['gpp'	   ]	= '/home/ud4/CESM1-BGC/gpp/without_land_use_change/'  #just for extraction of time variable
paths['out']				= {}
paths['out']['dri_cc']		= {}

nc_data						= {}
nc_data['gpp' 	  ]			= nc4.Dataset(paths['in']['gpp'    	]   +'cesm1bgc_pftcon_gpp_anomalies_gC.nc'	)
if th_type == 'independent':

	if args.model_config	== 'wo_lulcc' :	
		paths['in' ][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/'
		paths['out'][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/results/'
		nc_data[th_type ]		= nc4.Dataset(paths['in']['wo_lulcc']   +'dominant_driver_correlation.nc'		)
	if args.model_config	== 'w_lulcc' :		
		paths['in' ][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/'
		paths['out'][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/results/'
		nc_data[th_type ]	= nc4.Dataset(paths['in']['w_lulcc' ]   +'dominant_driver_correlation.nc'		)

if th_type  == 'misc' :
	# Following path is for the old attribution method
	#paths['in' ][th_type]           = '/home/ud4/CESM1-BGC/corr_triggers/misc/' + case_name + '/'
	# Following path is for the new attribution method '04/10/2019'
	paths['in' ][th_type]           = '/home/ud4/CESM1-BGC/attr_2019/misc/dom_driver/' + case_name + '/'
	paths['out'][th_type]           = paths['in' ]['misc'] + 'basemaps_dataframes/'
	nc_data[th_type]		=	nc4.Dataset( paths['in' ]['misc'] + 'dominant_driver_correlation.nc' )
	
lat			= nc_data[th_type] .variables['lat' 	  ]
lon			= nc_data[th_type] .variables['lon' 	  ]
ranks		= nc_data[th_type]	.variables['rank'	  ]
wins		= nc_data[th_type]	.variables['win'   	  ]
lags		= nc_data[th_type] .variables['lag' 	  ]
dom_dri_ids	= nc_data[th_type] .variables['dri_id'	  ]
dom_dri_cc	= nc_data[th_type]	.variables['dri_coeff']

time        = nc_data['gpp']['time'] [...]
win_len     = 25*12


dates_ar    = time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates = [dates_ar[i*win_len] for i in range(int(time.size/win_len))] #list of the start dates of all 25 year windows
end_dates   = [dates_ar[i*win_len+win_len-1] for i in range(int(time.size/win_len))] #list of the end dates of all 25 year windows


for win in np.asarray(wins[...], dtype =int):
	for lg in np.asarray(lags[...],dtype = int):
		for rk in np.asarray(ranks[...][0:1], dtype = int):
			counts 			= np.unique(np.ma.masked_equal(np.ma.masked_invalid(dom_dri_ids[rk,win,lg,:,:]),0),return_counts=True)  # there are np.nans and 0's in the array that have to be masked.
			counts_drivers 	= np.array([counts[1][i] for i in range(counts[1].size)])
			#since many drivers were not dominant for most part so only limiting the plot to the relevant ones
			print "counts for dom rank %s and lag %s...:"%(format(rk,'002'), format(lg,'002'))
			tmp_drivers_code	= np.copy(drivers_code)
			for d in counts[0].data:
				tmp_drivers_code = np.ma.masked_equal (tmp_drivers_code, d)
		
			df_counts 		= pd.DataFrame({'Counts':counts_drivers[:-1]})  #the last value corresponds to the masks
			#print "Dataframe", df_counts
			#print "Drivers",drivers, drivers.shape
			#print "mask...shape.",tmp_drivers_code.mask.shape

			df_counts.index = drivers [tmp_drivers_code.mask]
			perc = [round(i*100./sum(df_counts['Counts'].values),2) for i in df_counts['Counts'].values]
			df_counts['percentage']=perc
			#Calculating the mean and std of the climate drivers
			mean_cc	= []
			std_cc 	= []
			for code_id in drivers_code[tmp_drivers_code.mask]:
				#print "code_ID...", code_id
				mean_cc.append(np.ma.mean(dom_dri_cc[rk,win,lg,:,:][~np.ma.masked_not_equal(dom_dri_ids[rk,win,lg,:,:],code_id).mask]))
				std_cc.append(np.ma.std(dom_dri_cc[rk,win,lg,:,:][~np.ma.masked_not_equal(dom_dri_ids[rk,win,lg,:,:],code_id).mask]))
			df_counts['mean_coeff'] = mean_cc
			df_counts['std_coeff'] 	= std_cc

			#print mean_cc
			print 'dataframe_win_%s_lag_%s_and_rank_%s.csv'%(format(win,'02'),format(lg,'02'),format(rk,'02'))
			df_counts .to_csv(paths['out'][th_type]+ 'dataframe_win_%s_lag_%s_and_rank_%s.csv'%(format(win,'02'),format(lg,'02'),format(rk,'02')),sep=',')
		
			del counts, counts_drivers, tmp_drivers_code, df_counts, perc, mean_cc, std_cc

			#Saving table only for CONUS
			top     = 49.345    +.5
			bottom  = 24.743    -.5
			left    = 360-124.78-.5
			right   = 360-66.95 +.5

			counts			= np.unique(np.ma.masked_equal(np.ma.masked_invalid(dom_dri_ids[rk,win,lg,geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])]),0),return_counts=True)
			counts_drivers  = np.array([counts[1][i] for i in range(counts[1].size)])
			tmp_drivers_code= np.copy(drivers_code)
			for d in counts[0].data:
				tmp_drivers_code = np.ma.masked_equal (tmp_drivers_code, d)
			df_counts       = pd.DataFrame({'Counts':counts_drivers[:-1]})
			df_counts.index = drivers [tmp_drivers_code.mask]
			perc = [round(i*100./sum(df_counts['Counts'].values),2) for i in df_counts['Counts'].values]
			df_counts['percentage']=perc
			#Calculating the mean and std of the climate drivers
			mean_cc = []
			std_cc  = []
			for code_id in drivers_code[tmp_drivers_code.mask]:
		#	print "code_ID...", code_id
				mean_cc.append(np.ma.mean	(dom_dri_cc[rk,win,lg,geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])][~np.ma.masked_not_equal(dom_dri_ids[rk,win,lg,geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])],code_id).mask]))
				std_cc.append(np.ma.std		(dom_dri_cc[rk,win,lg,geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])][~np.ma.masked_not_equal(dom_dri_ids[rk,win,lg,geo_idx(bottom,lat[...]):geo_idx(top,lat[...]),geo_idx(left,lon[...]):geo_idx(right,lon[...])],code_id).mask]))
			df_counts['mean_coeff'] = mean_cc
			df_counts['std_coeff']  = std_cc
			df_counts .to_csv(paths['out'][th_type] + 'dataframe_CONUS_win_%s_lag_%s_and_rank_%s.csv'%(format(win,'02'),format(lg,'02'),format(rk,'02')),sep=',')


			from    mpl_toolkits.basemap import Basemap
			from    matplotlib import cm
			import matplotlib.patches as patches
	
			fig,ax  = plt.subplots(figsize = (6,2.6),tight_layout=True,dpi=500)
	
			bmap    = Basemap(  projection  =   'eck4',
							lon_0       =   0,
							resolution  =   'c',   
							)
			LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij') 														###################
			#LAT,LON = np.meshgrid(lat[...][60], lon[...],indexing ='ij')
			cmap 	= plt.get_cmap('rainbow', drivers_code.size)
			ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_equal(np.ma.masked_invalid(dom_dri_ids[rk,win,lg,:,:]),0),latlon=True,cmap=cmap,vmin = 5, vmax = drivers_code[-1]+5)
			cbar 	= plt.colorbar(ax ,ticks = range(drivers_code[0],drivers_code[-1]+1,10))
			cbar	.ax.set_yticklabels(drivers_names)
			cbar	.ax.set_ylabel('Drivers',fontsize =12)
			bmap	.drawparallels(np.arange(-90., 90., 30.),fontsize=10, linewidth = .2)
			bmap	.drawmeridians(np.arange(0., 360., 60.),fontsize=10, linewidth = .2)
			bmap	.drawcoastlines(linewidth = .2)
			text_win = []
			for i in start_dates:
				a = i.year
				b = a+24
				c = str(a)+'-'+str(b)[2:]
				text_win.append(c)
	
			plt.title (text_win[win])
			fig.savefig(paths['out'][th_type] + 'basemap_win_%s_lag_%s_dom_rank_%s.pdf'%(format(win,'02'),format(lg,'02'),format(rk,'02')))	
			#fig.savefig(paths['out'][th_type] + case_name+ 'basemap_win_%s_lag_%s_dom_rank_%s.png'%(format(win,'02'),format(lg,'02'),format(rk,'02')))	



