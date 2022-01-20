# Bharat Sharma
# python 2.7
# 

import  matplotlib as mpl
mpl.use('Agg')
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab as plt
from functions import time_dim_dates, index_and_dates_slicing, geo_idx, mpi_local_and_global_index, create_seq_mat, cumsum_lagged,patch_with_gaps_and_eventsize, norm
import netCDF4 as nc4
import datetime as dt
import argparse
import pandas as pd

driver_consider = 7

drivers						 = np.array(['prcp','sm','pme','tmax','tsa','tmin','col_fire_closs','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider] 
#drivers						 = np.array(['prcp','sm','pme','tmax','tsa','tmin','fclosscol','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider] 
drivers_names				 = np.array(['Prcp','Soilmoist',r'$P-ET$',r'$T_{max}$',r'$T_{sa}$',r'$T_{min}$','Fire','SPI_gamma_06','SPI_pearson_06','SPEI_gamma_06','SPEI_gamma_12','SPEI_pearson_06'])[:driver_consider]
drivers_code				 = np.array([  10,   20,   30,   40,       50,		60,			70,			80,		90,			100, 110,120 ])[:driver_consider]
features 			 		 = {}
#features['abr']				 = np.array(['prcp','sm','pme','tmax','tsa', 'tmin','fclosscol','gpp','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider]
features['abr']				 = np.array(['prcp','sm','pme','tmax','tsa', 'tmin','col_fire_closs','gpp','spi_g06','spi_p06','spei_g06','spei_g12','spei_p06'])[:driver_consider]

# The name with which the variables are stored in the nc files:
features['Names']			 	= {}
features['Names']['prcp']	 	= 'prcp'
features['Names']['sm']	 	 	= 'sm'
features['Names']['pme']	 	= 'pme'
features['Names']['tmax']	 	= 'tmax'
features['Names']['tsa']	 	= 'tsa'
features['Names']['tmin']	 	= 'tmin'
features['Names']['col_fire_closs']	= 'col_fire_closs'
features['Names']['gpp']	 	= 'gpp'
features['Names']['spi_g06'] 	= 'spi_gamma_06'

parser  = argparse.ArgumentParser()
parser.add_argument('--model_config', '-config' 	, help = "Model Configuration"   						, type= str, default= 'w_lulcc'	) # w_lulcc, wo_lulcc, all
parser.add_argument('--plot_win'	, '-pwin'		, help = "which time period to plot? 2000-24:win 06" 	, type= int, default=  6		) # 0 to 17
parser.add_argument('--plot_lag'	, '-plag'		, help = "which (cum) lag in months to plot? " 			, type= int, default=  1		) # 0 to 
parser.add_argument('--thres_type'  , '-th_typ'   	, help = "Type of thresholds (independent or misc)" 	, type= str, default= 'misc'    ) # this is mainly for the misc case study filter
parser.add_argument('--ext_type'    , '-ext_typ'  	, help = "Type of extreme analysis (pos or neg)"		, type= str, default= 'neg'     )
args = parser.parse_args()

#Running: run ecp_dom_driver_RGB_tce.py -config w_lulcc -th_typ misc -ext_typ neg -pwin 6 -plag 1

th_type     = args.thres_type
ext_type    = args.ext_type
conf        = args.model_config

#temp
print args
# Misc case (th_type_ for this code is basically based on the 'keys_wrt_wo_lulcc_pos' 
keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
case_name   = ext_type+'_'+conf+'_based_pos_wo_lulcc'

paths						= {}
paths['in' ]				= {}
paths['in' ]['gpp'	   ]	= '/home/ud4/CESM1-BGC/gpp/without_land_use_change/'  #just for extraction of time variable
paths['out']				= {}
paths['out']['dri_cc']		= {}

nc_data						= {}
nc_data['gpp' 	  ]			= nc4.Dataset(paths['in']['gpp'    	]   +'cesm1bgc_pftcon_gpp_anomalies_gC.nc'	)
#if th_type == 'independent':
#	if args.model_config	== 'wo_lulcc' :	
#		paths['in' ][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/'
#		paths['out'][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/rgb_plots/'
#		paths['out']['dri_cc'][th_type] = '/home/ud4/CESM1-BGC/corr_triggers/without_land_use_change/dri_cc_plot/'
#		nc_data[th_type ]		= nc4.Dataset(paths['in']['wo_lulcc']   +'dominant_driver_correlation.nc'		)
#	if args.model_config	== 'w_lulcc' :		
#		paths['in' ][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/'
#		paths['out'][th_type] 			= '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/rgb_plots/'
#		paths['out']['dri_cc'][th_type] = '/home/ud4/CESM1-BGC/corr_triggers/with_land_use_change/dri_cc_plot/'
#		nc_data[th_type ]	= nc4.Dataset(paths['in']['w_lulcc' ]   +'dominant_driver_correlation.nc'		)
if th_type  == 'misc' :
	paths['in' ][th_type]           = '/home/ud4/CESM1-BGC/attr_2019/misc/dom_driver/'+ case_name + '/'
	paths['out'][th_type]           = paths['in' ]['misc'] + 'rgb_plots/'
	paths['out']['dri_cc'][th_type] = paths['in' ]['misc'] + 'dri_cc_plot/'
	nc_data[th_type]		=	nc4.Dataset( paths['in' ]['misc'] + 'dominant_driver_correlation.nc' )
	
lat			= nc_data[th_type] .variables['lat' 	  ]
lon			= nc_data[th_type] .variables['lon' 	  ]
ranks		= nc_data[th_type]	.variables['rank'	  ]
wins		= nc_data[th_type]	.variables['win'   	  ]
lags		= nc_data[th_type] .variables['lag' 	  ]
dom_dri_ids	= nc_data[th_type] .variables['dri_id'	  ]
dom_dri_cc	= nc_data[th_type]	.variables['dri_coeff']

time 		= nc_data['gpp']['time'] [...]
win_len     = 25*12

dates_ar    = time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates = [dates_ar[i*win_len] for i in range(int(time.size/win_len))] #list of the start dates of all 25 year windows
end_dates   = [dates_ar[i*win_len+win_len-1] for i in range(int(time.size/win_len))] #list of the end dates of all 25 year windows

pwin 		= args.plot_win
plag		= args.plot_lag

""" RGB dic:
	R(fire) : Fire 
	B(water) : max(Prcp, PME,SM)
	G(temp) : max(Tmax, Tmin, Tsa)
	 u'10:prcp,20:sm,30:pme,40:tmax,50:tsa,60:tmin,70:fclosscol,'
"""
rgb_val			= {}
rgb_val['fire'] = np.ma.masked_not_equal(dom_dri_ids[:,pwin,plag,:,:],	70	 )     	#Fire id is 70
rgb_val['water']= np.ma.masked_outside	(dom_dri_ids[:,pwin,plag,:,:],	9 ,31)		# Water includes Prcp, P-ET, SM and those ids are 10,20,30
rgb_val['temp']	= np.ma.masked_outside	(dom_dri_ids[:,pwin,plag,:,:],	39,61)		# Temp includes Tmax,Tsa and Tmin and those ids are 40,50,60
rgb_val['r']	=  np.ma.max(np.ma.abs(np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:],mask = rgb_val['fire'] [...].mask)), axis = 0)
rgb_val['g']	=  np.ma.max(np.ma.abs(np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:],mask = rgb_val['temp'] [...].mask)), axis = 0)
rgb_val['b']	=  np.ma.max(np.ma.abs(np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:],mask = rgb_val['water'] [...].mask)), axis = 0)

rgb_val['color']		= {}
rgb_val['color']['r']	= 'fire'
rgb_val['color']['g']	= 'temp'
rgb_val['color']['b']	= 'water'

#print a
		
rgb_val['rgb']	= np.ma.asarray([rgb_val['r'],rgb_val['g'],rgb_val['b']])  	# the shape right now is 3,192,288 
rgb_val['rgb']	= np.rollaxis(rgb_val['rgb'],0,3)							# the shape right now is 192,288,3 and this is what is needed for RGB plotting

_all_mask		= rgb_val['rgb'].mask.all(axis=2) #'all' tests for the condition where all values are 'True'
rgb_val['rgb'].data[np.where(rgb_val['rgb'].mask)] = 0 # Replacing all the mask values with 0; places where there is a valid 'r' and 'g' and masked 'b' will replace 'b' with 0.
_all_mask		= _all_mask[...,np.newaxis]*np.ones(3,dtype=bool)
rgb_val['rgb'].data[np.where(_all_mask)]	= 1

# New: Count and sign of dominant drivers
# ---------------------------------------
tmp 		= np.ma.asarray	( [rgb_val['r'],rgb_val['g'],rgb_val['b']])
lt_ln_mat	= create_seq_mat( lat.size, lon.size)	# creating a 
mask_1d 	= lt_ln_mat[~np.sum	( tmp,axis = 0) .mask]
"""
Creating a DataFrame to store the dominant drivers and their correlation sign

1. We created a database will of location and the dominant drivers contribution to RGB map, their CC and sign
2. The saved entries are: 
	-	location (lat,lon)
	- 	Dominant Driver Ids  ('10:prcp,20:sm,30:pme,40:tmax,50:tsa,60:tmin,70:fclosscol,')
	-	Dominant Driver CC
	-	Compound Drivers in order : Fire-Temp-Water (RGB)
		o	The Positive correlation is assigned a number 4
		o	The Negative correlation is assigned a number 2
		o	The Non-significant driver pool is assigned a number 9
3. Pools Drivers case for negative TCE events:
	-	299 : Pure Fire		fire
	-	929	: Pure Temp		hot
	-	994	: Pure Water	dry
	-	949	:				cold
	-	992	:				wet
	- 	924	:				hot & dry
	-	944	:				cold & dry
	-	922	:				hot & wet
	-	224	:				hot & dry & fire
	-	244 :               cold & dry & fire
	-	222 :               hot & wet & fire
	-	229	:				hot & not dry & fire
	-	249 :				cold & not dry & fire
	-	294	:				not hot & dry & fire
"""
Compound_Pools_Names	= {}
Compound_Pools_Names [299] = 'fire'
Compound_Pools_Names [929] = 'hot'
Compound_Pools_Names [994] = 'dry'
Compound_Pools_Names [949] = 'cold'
Compound_Pools_Names [992] = 'wet'
Compound_Pools_Names [924] = 'hot & dry'
Compound_Pools_Names [944] = 'cold & dry'
Compound_Pools_Names [922] = 'hot & wet'
Compound_Pools_Names [224] = 'hot & dry & fire'
Compound_Pools_Names [244] = 'cold & dry & fire'
Compound_Pools_Names [222] = 'hot & wet & fire'
Compound_Pools_Names [229] = 'hot & not-dry & fire'
Compound_Pools_Names [249] = 'cold & not-dry & fire'
Compound_Pools_Names [294] = 'not-hot & dry & fire'


df_null = 	np.ma.masked_all((mask_1d.size,9))
df_pool	=	pd.DataFrame(data = df_null, columns = ['lat','lon','Fire','Water','Temp','F_val','W_val','T_val','Pools'])
"""
we will create a database will of location and the dominant drivers contribution to RGB map, their CC and sign
"""

# for a particular lag and win
for idx, pixel in enumerate(mask_1d):
	lt,ln	= np.argwhere(lt_ln_mat == pixel)[0]
	print lt,ln
	# Fire
	f_tmp	= np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:],mask = rgb_val['fire'] [...].mask)[:,lt,ln]
	if f_tmp.max() == f_tmp.min():
		f_val = f_tmp.max()
		if f_val > 0: f_sign = 'pos'
		else: f_sign = 'neg'
	elif f_tmp.max() >= abs(f_tmp.min()):
		f_val = f_tmp.max()
		f_sign = 'pos'
	elif f_tmp.max() < abs(f_tmp.min()):
		f_val = f_tmp.min()
		f_sign = 'neg'
	else:
		f_val = None;	f_sign = None;	f_id = None;	f = None

	if f_sign == 'pos': 
		f = 70
		FS = '4'
	elif f_sign == 'neg': 
		f = -70
		FS = '2'
	else:
		FS = '9'

	# temp
	t_tmp   = np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:],mask = rgb_val['temp'] [...].mask)[:,lt,ln]
	t_ids   = np.ma.masked_array(data=dom_dri_ids[:,pwin,plag,:,:],mask = rgb_val['temp'] [...].mask)[:,lt,ln]

	if t_tmp.max() == t_tmp.min():
		t_val = t_tmp.max()
		if t_val > 0: t_sign = 'pos'
		else: t_sign = 'neg'
	elif t_tmp.max() >= abs(t_tmp.min()):
		t_val = t_tmp.max()
		t_sign = 'pos'
	elif t_tmp.max() < abs(t_tmp.min()):
		t_val = t_tmp.min()
		t_sign = 'neg'
	else:
		t_val = None;	t_sign = None;	t_id = None;	t = None
	
	if t_val != None:
		t_id = np.argwhere(t_tmp==t_val)[0]
	
	if t_sign == 'pos': 
		t = int(t_ids[t_id])
		TS = '4'
	elif t_sign == 'neg': 
		t = -int(t_ids[t_id])
		TS = '2'
	else:
		TS = '9'

	# water
	w_tmp   = np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:],mask = rgb_val['water'] [...].mask)[:,lt,ln]
	w_ids   = np.ma.masked_array(data=dom_dri_ids[:,pwin,plag,:,:],mask = rgb_val['water'] [...].mask)[:,lt,ln]

	if w_tmp.max() == w_tmp.min():
		w_val = w_tmp.max()
		if w_val > 0: w_sign = 'pos'
		else: w_sign = 'neg'
	elif w_tmp.max() >= abs(w_tmp.min()):
		w_val = w_tmp.max()
		w_sign = 'pos'
	elif w_tmp.max() < abs(w_tmp.min()):
		w_val = w_tmp.min()
		w_sign = 'neg'
	else:
		w_val = None;	w_sign = None;	w_id = None;	w = None

	if w_val != None:
		w_id = np.argwhere(w_tmp==w_val)[0]
	
	if w_sign == 'pos': 
		w = int(w_ids[w_id])
		WS = '4'
	elif w_sign == 'neg': 
		w = -int(w_ids[w_id])
		WS = '2'
	else:
		WS = '9'
	pool_sign	= FS+TS+WS

	df_pool.iloc[idx,:] = np.array([lt,ln,f,w,t,f_val,w_val,t_val,int(pool_sign)])
	del lt,ln,f,w,t,f_val,w_val,t_val,pool_sign


print (" Total Fire  +/-\t%d/%d"%(np.ma.sum(np.ma.array(df_pool['Fire'])>0),np.sum(np.ma.array(df_pool['Fire'])<0)))
print (" Total Water +/-\t%d/%d"%(np.ma.sum(np.ma.array(df_pool['Water'])>0),np.sum(np.ma.array(df_pool['Water'])<0)))
print (" Total Temp  +/-\t%d/%d"%(np.ma.sum(np.ma.array(df_pool['Temp'])>0),np.sum(np.ma.array(df_pool['Temp'])<0)))
	

compound_pools_all	= np.asarray(df_pool['Pools'], dtype =int) # list of all the compound pools present/ obsereved
compound_pools_counts_id, compound_pools_count_tot = np.unique (compound_pools_all, return_counts = True)

# Filtering 95 percentile of total compound pools
cp_class 	= []  # Compound Pools Names/Classification
cp_total 	= []  # Compound Pools Total Number
cp_frac		= []  # Compound Pools Fraction of Total
cp_ids		= []
compound_pools_dict_4plot = {}
for c_idx, c_id in enumerate (compound_pools_counts_id):
	th_5 = compound_pools_all.size*0.05
	if compound_pools_count_tot [c_idx] > th_5:
		if c_id in Compound_Pools_Names.keys():
			c_name = Compound_Pools_Names[c_id]
		else:
	 		c_name = str(c_id)
		cp_class .append (c_name)
		cp_ids	 .append (c_id)
		cp_total .append (compound_pools_count_tot [c_idx])
		cp_frac	 .append (compound_pools_count_tot [c_idx] / float(compound_pools_all.size))

# Saving it as a csv file:
dict_cp = {}
dict_cp ['class'] 	= cp_class
dict_cp ['IDs']		= cp_ids
dict_cp ['total']	= cp_total
dict_cp ['frac']	= cp_frac

df_cp	= pd.DataFrame(dict_cp)
df_cp = df_cp.sort_values(by = ['frac'], ascending = False)
df_cp.to_csv(paths['out'][th_type]+'compound_driver_frac_win_%s_lag_%s.csv'%(format(pwin,'02'),format(plag,'02')), index=False)



# New End -------------------------------





#The basemap for some reason, is not able to change the center of the map hence the data matrix have to be flipped mannually
data			= np.zeros(rgb_val['rgb'].data.shape)
data[:,:144,:]	= rgb_val['rgb'].data[:,144:,:]
data[:,144:,:]	= rgb_val['rgb'].data[:,:144,:]

# Ploting the RGB Map
# -------------------
# if you want to add a circular rgb color wheel? say 'yes'
rgb_wheel = 'no'

# -------------------
fig,ax  = plt.subplots(figsize = (6,2.8),tight_layout=True,dpi=500)
bmap = Basemap(projection='eck4',
		    	lon_0= 0,
				resolution = 'c')
#if the lats/lons are not +1/+1 of shape data, the colors are screwed
lons 	= np.linspace(0,360,data.shape[1]+1)
lats	= np.linspace(-90, 90,data.shape[0]+1)
LON,LAT = np.meshgrid(lons, lats)
ax = bmap.pcolormesh(LON,LAT,data[...,0],color=data.reshape((-1,3)),latlon=True,linewidth=0)
ax.set_array(None)
bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
bmap    .drawcoastlines(linewidth = .25,color='lightgrey')

if rgb_wheel == 'yes':
	display_axes = fig.add_axes([0.1,0.3,0.2,0.2], projection='polar')
	display_axes._direction = 2*np.pi
	norm = mpl.colors.Normalize(0.0, 2*np.pi)
	quant_steps = 2056
	cb = mpl.colorbar.ColorbarBase(display_axes, cmap=cm.get_cmap('hsv',quant_steps),
									norm=norm,
									orientation='horizontal')
	cb.outline.set_visible(False) 
	#for j, lab in enumerate(['R','B','G']):
	cb.ax.text(0.98, 1.05, 'Fire', fontsize=6)
	cb.ax.text(.35, 1.8, 'Temp', fontsize=6, rotation=120)
	cb.ax.text(.63, 1.50, 'Water',fontsize=6, rotation=60)
	
	#cb.ax.set_yticklabels(["a","b","c","d","e","f"])
	display_axes.set_axis_off()

fig.savefig( paths['out'][th_type]+'dom_dri_compound_RGB_win_%s_lag_%s.pdf'%(format(pwin,'02'),format(plag,'02'))) 
plt.close(fig)
print "Done", paths['out'][th_type]+'dom_dri_compound_RGB_win_%s_lag_%s'%(format(pwin,'02'),format(plag,'02'))


# Spatial plot of individual driver correlatons
# --------------------------------------------
# In case you want to plot individual driver correlations : ind_plot = 'yes' 
ind_plot = 'no'

if ind_plot == 'yes':
	#if pwin == 6:
	for idx,dri in enumerate (drivers_names):

		fig1,ax  = plt.subplots(figsize = (6,2.8),tight_layout=True,dpi=500)
		bmap = Basemap(projection='eck4',
						lon_0= 0,
						resolution = 'c')
		#if the lats/lons are not +1/+1 of shape data, the colors are screwed
		
		LON,LAT = np.meshgrid(lon[...],lat[...])
		plot_data = np.ma.max(np.ma.masked_array(data=dom_dri_cc[:,pwin,plag,:,:], mask = np.ma.masked_not_equal(dom_dri_ids[:,pwin,plag,:,:], drivers_code[idx]) .mask),axis = 0)
		ax  	= bmap.pcolormesh(LON,LAT,plot_data,cmap= cm.PuOr,latlon=True,linewidth=0,vmax=1,vmin=-1)	
		bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
		bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
		bmap    .drawcoastlines(linewidth = .25,color='lightgrey')
		cbar    = plt.colorbar(ax)
		plt.title("%s"%dri)
		fig1.savefig( paths['out']['dri_cc'][th_type]+'win_%s/dri_win_%s_lag_%s_%s.pdf'%(format(pwin,'02'),format(pwin,'02'),format(plag,'02'),drivers[idx])) 
		plt.close(fig1)


"""
colors = [1,2,3]
with nc4.Dataset(paths['out'][conf]+'dom_dri_rgb_values_win_%s_lag_%s.nc'%(format(pwin,'02'),format(plag,'02')),mode="w") as dset:
	dset        .createDimension( "lat" ,size = lat.size)
	dset        .createDimension( "lon" ,size = lon.size)
	dset        .createDimension( "colors" ,size = len(colors))
	t   =   dset.createVariable(varname = "colors"        ,datatype = float   , dimensions = ("colors")             ,fill_value = -999)
	y   =   dset.createVariable(varname = "lat"         ,datatype = float   , dimensions = ("lat")              ,fill_value = -999)
	x   =   dset.createVariable(varname = "lon"         ,datatype = float   , dimensions = ("lon")              ,fill_value = -999)
	z   =   dset.createVariable(varname = "dri_id"      ,datatype = float   , dimensions = ("colors","lat","lon") ,fill_value = -999)
	x.axis  =   "X"
	y.axis  =   "Y"
	t.axis  =   "T"
	t[...]  =   colors
	x[...]  =   lon[...]
	y[...]  =   lat[...]
	z[...]  =   rgb_val['rgb']
	x.units         =   lat.units
	x.missing_value =   -999
	x.setncattr         ("long_name",lat.long_name)
	y.units         =   lon.units
	y.missing_value =   -999
	y.setncattr         ("long_name",lon.long_name)
	z.missing_value =   -999
	z.stardard_name =   "RGB Colors"
	z.setncattr         ("long_names","R,G,B")

"""
