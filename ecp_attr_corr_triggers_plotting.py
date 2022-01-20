#	Bharat Sharma (bharat.sharma.neu@gmail.com)
# python 2.7
# plotting percent distribution of dominant drivers
# --------------------------------------------------

import numpy as np              
import pandas as pd
import  matplotlib as mpl
mpl.use('Agg')

import  matplotlib.pyplot as plt
import argparse

parser  = argparse.ArgumentParser()
parser.add_argument ('--cum_lag'	,'-lag'		, help = 'cum lag months? (multiple lag optional) use , to add multiple'		, type = str	, default = '01,02,03'	) 
parser.add_argument ('--dom_rank'	,'-rank'	, help = 'dominant driver rank'													, type = str	, default = '00'		)
parser.add_argument('--model_config', '-config'	, help = "Model Configuration"   												, type = str	, default = 'wo_lulcc'	) # w_lulcc, wo_lulcc, all
args = parser.parse_args()

""" remarks:
1.	lag	: cummulative lag effects are calculated for the months 00 to 12.
2.	rank: rank are calculated from 00 to 05: 00 is the most dominant driver
"""

lg = (args.cum_lag).split(',')
rk = args.dom_rank
if args.model_config	== 'wo_lulcc' :		
	configs	= ['wo_lulcc']
	conf 	= 'wo_lulcc'
if args.model_config	== 'w_lulcc' :		
	configs	= ['w_lulcc']
	conf	= 'w_lulcc'
if args.model_config	== 'all' :		
	configs	= ['wo_lulcc','w_lulcc']

paths						= {}
paths['in' ]				= {}
paths['in' ]['wo_lulcc'] 	= '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/results/'
paths['in' ]['w_lulcc' ] 	= '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/results/'
paths['in' ]['gpp'	   ]	= '/home/ud4/CESM1-BGC/gpp/without_land_use_change/'
paths['out']				= {}
paths['out']['wo_lulcc'] 	= '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/per_dom/'
paths['out']['w_lulcc' ] 	= '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/per_dom/'

wins 	= [format(i,'02') for i in range(18)]
lags 	= [format(i,'02') for i in range(13)]
ranks 	= [format(i,'02') for i in range(6)]

filenames={}
#format>>> filenames [wins] [lags] [ranks]
for w in wins:
	filenames[w] = {}
	for l in lags[1:7]:
		filenames[w][l] = {}
		for r in ranks[:4]:
			filenames[w][l][r] = 'dataframe_win_%s_lag_%s_and_rank_%s.csv'%(w,l,r) 

in_yr 	= 1850
win_yr 	= [str(in_yr+i*25) + '-'+str(in_yr +(i+1)*25-1)[2:] for i in range(18)]

driver_names = ['Prcp','Soilmoist',r'$P-ET$',r'$T_{max}$',r'$T_{sa}$',r'$T_{min}$','Fire']
#driver_names = ['Prcp', 'Soiloist','PME','Tmax','Tsa', 'Tmin', 'Fire']

if len(lg)==1:

	data_percent = np.zeros((len(win_yr), len(driver_names)))
	print "data_percent shape: ", data_percent.shape
	print "data shape", np.transpose(pd.read_csv(paths['in'][conf]+filenames[w][lg][rk])['percentage']).shape

	df = pd.DataFrame( data_percent , index = win_yr, columns = driver_names) #main dataframe

	for w in wins:
		data = pd.read_csv(paths['in'][conf]+filenames[w][lg][rk])
		drivers = data.iloc[:,0]
		data_df	= pd.DataFrame(data = pd.read_csv(paths['in'][conf]+filenames[w][lg][rk])['percentage'].values.reshape(1,len(drivers)), index = [win_yr[int(w)]], columns = drivers)  #dataframe for a particuar windw
		for idx,dr in enumerate (drivers):
			df.loc[data_df.index,driver_names[idx]] = data_df[dr]


	fig,ax	= plt.subplots(figsize = (8,5), tight_layout=True, dpi=400)
	for i in range(len(driver_names)):
		ax.plot( range(18), df.iloc[:,i], label = driver_names[i], linewidth = 2)
	ax.set_xticks(range(18))
	ax.set_xticklabels(df.index,fontsize =12)
	for tick in ax.get_xticklabels():
		tick.set_rotation(90)
	ax.set_xlabel('Time ($25-yr)$ wins',fontsize =16)
	ax.legend(loc = 'upper center',ncol=len(driver_names), bbox_to_anchor=(.5,1.08),frameon =False) 
	ax.set_ylabel('Percent of dominant driver', fontsize = 16)
	fig.savefig(paths['out'][conf]+'percent_dominance_lag_%s_rank_%s.pdf'%(lg,rk))
	fig.savefig(paths['out'][conf]+'percent_dominance_lag_%s_rank_%s.png'%(lg,rk))

if len(lg)>1:
	data_lag = {}
	for LAG in lg:
		data_percent = np.zeros((len(win_yr), len(driver_names)))
		print "data_percent shape: ", data_percent.shape
		print "data shape", np.transpose(pd.read_csv(paths['in'][conf]+filenames[w][LAG][rk])['percentage']).shape

		df = pd.DataFrame( data_percent , index = win_yr, columns = driver_names) #main dataframe

		for w in wins:
			data = pd.read_csv(paths['in'][conf]+filenames[w][LAG][rk])
			drivers = data.iloc[:,0]
			data_df	= pd.DataFrame(data = pd.read_csv(paths['in'][conf]+filenames[w][LAG][rk])['percentage'].values.reshape(1,len(drivers)), index = [win_yr[int(w)]], columns = drivers)  #dataframe for a particuar windw
			for idx,dr in enumerate (drivers):
				df.loc[data_df.index,driver_names[idx]] = data_df[dr]
		data_lag[LAG] = df
	
#Plotting subplots
	color_list 		= ['b','b','b','g','g','g','r']
	linestyle_list 	= [':','-','--','-',':','--','-',]
	fig,ax  = plt.subplots(nrows=3,ncols= 1,gridspec_kw = {'wspace':0, 'hspace':0.02},tight_layout = True, figsize = (6.5,7), dpi = 400)
	ax      = ax.ravel()
	for lag_idx,LAG in enumerate(lg):
		for dri_idx in range(len(driver_names)):
			ax[lag_idx].plot( range(18), data_lag[LAG].iloc[:,dri_idx], label = driver_names[dri_idx], color=color_list[dri_idx],linestyle=linestyle_list[dri_idx], linewidth = 1)
		ax[lag_idx].set_xticks(range(18))
		ax[lag_idx].tick_params(axis="x",direction="in")
		ax[lag_idx].set_xticklabels([])
		ax[lag_idx].set_ylabel("Lag: %s"%(LAG))
		ax[lag_idx].set_ylim([0,32])
	ax[lag_idx].set_xticklabels(df.index,fontsize =9)
	for tick in ax[lag_idx].get_xticklabels():
		tick.set_rotation(90)
	ax[lag_idx].set_xlabel('Time ($25-yr)$ wins',fontsize =14)
	fig.text(0.03, 0.5, 'Percent Area under Extreme Events', va='center', ha='center', rotation='vertical', fontsize=14)
	ax[0].legend(loc = 'upper center',ncol=len(driver_names), bbox_to_anchor=(.44,1.15),frameon =False,fontsize=9,handletextpad=0.1)
	plt.gcf().subplots_adjust(bottom=0.1)
	fig.savefig(paths['out'][conf]+'percent_dominance_multilag_%s_rank_%s.pdf'%('-'.join(lg),rk))
				 



		


print "Done .................... percent_dominance_lag_%s_rank_%s"%(lg,rk)
