# Bharat Sharma
# python 2.7
# The aim is to create a stacked graphs of gpp and climate drivers to identify the climate that led to an extreme event at that location
# The TCE events are also shaded and  we could focus on a particular TCE and see the detais
import numpy as np              
import pandas as pd
import  matplotlib as mpl
mpl.use('Agg')
from functions import time_dim_dates, index_and_dates_slicing, geo_idx
import datetime as dt
import netCDF4 as nc4
import  matplotlib.pyplot as plt
import argparse
from    mpl_toolkits.basemap import Basemap
from    matplotlib import cm

parser  = argparse.ArgumentParser()
parser.add_argument ('--lat_idx'	,'-lat_idx'	, help = ' Lat coordinate of the location'										, type = float	, default = 105	) 
parser.add_argument ('--lon_idx'	,'-lon_idx'	, help = ' Lon coordinate of the location'										, type = float	, default = 29	) 
#parser.add_argument ('--cum_lag'	,'-lag'		, help = 'cum lag months? (multiple lag optional) use , to add multiple'		, type = str	, default = '01,02,03'	) 
#parser.add_argument ('--dom_rank'	,'-rank'	, help = 'dominant driver rank'													, type = str	, default = '00'		)
parser.add_argument('--model_config', '-config'	, help = "Model Configuration"   												, type = str	, default = 'wo_lulcc'	) # w_lulcc, wo_lulcc, all
parser.add_argument('--ext_type'    ,   '-ext_typ'  , help = "Type of extreme analysis (pos or neg)"    , type= str, default= 'neg'    )
parser.add_argument('--thres_type'  ,   '-th_typ'   , help = "Type of thresholds (independent or misc)" , type= str, default= 'misc'   )
parser.add_argument('--variables',   '-vars'  , help = "Type all the variables you want to plot( separated by '-' minus sign" , type= str, default= "prcp-sm-tmax-col_fire_closs")
parser.add_argument('--plot_win'    , '-pwin'       , help = "which time period to plot? 2000-24:win 06"    , type= int, default=  6   )
parser.add_argument('--plot_ano'    , '-pano'       , help = "Do you want to plot the anomalies also"    , type= str, default=  'n'   )
args = parser.parse_args()

# Running the code
# run ecp_shaded_extremes_drivers.py -pwin 6 -pano y -th_typ misc -vars prcp-sm-tmax-col_fire_closs -lat_idx 102 -lon_idx 26 -config wo_lulcc -ext_typ neg
pwin    = args.plot_win
lat_in	= args.lat_idx
lon_in	= args.lon_idx
conf	= args.model_config
ext_type= args.ext_type
th_type	= args.thres_type
variables_list =  (args.variables).split("-")
plt_ano	= args.plot_ano

#reading driver filenames!
nc_data							= {}
nc_data['in']						= {}
nc_data['in' ][conf]				= {}

if conf in ['wo_lulcc','pftcon']:
	for var in variables_list : 
		nc_data['in' ][conf][var]	 = nc4.Dataset('/home/ud4/CESM1-BGC/'+var+'/without_land_use_change/cesm1bgc_pftcon_'+ var+'.nc')[var]
		if plt_ano in ['y','Y','yes']:
			nc_data['in' ][conf][var+'_ano']	 = nc4.Dataset('/home/ud4/CESM1-BGC/'+var+'/without_land_use_change/cesm1bgc_pftcon_'+ var+'_anomalies.nc')[var]
if conf in ['w_lulcc','lulcc']:
	for var in variables_list : 
		nc_data['in' ][conf][var]	 = nc4.Dataset('/home/ud4/CESM1-BGC/'+var+'/with_land_use_change/cesm1bgc_with_lulcc_'+ var+'.nc')[var]
		if plt_ano in ['y','Y','yes']:
			nc_data['in' ][conf][var+'_ano']	 = nc4.Dataset('/home/ud4/CESM1-BGC/'+var+'/without_land_use_change/cesm1bgc_with_lulcc'+ var+'_anomalies.nc')[var]

if conf in ['wo_lulcc','pftcon']:
	nc_gpp		= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_gC.nc')['gpp']
	nc_gpp_ano	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_anomalies_gC.nc')['gpp']
if conf in ['w_lulcc','lulcc']:
	nc_gpp		= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_with_lulcc_gpp_gC.nc')['gpp']
	nc_gpp_ano	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_with_lulcc_gpp_anomalies_gC.nc')['gpp']


time	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_gC.nc')['time'][...]
lat		= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_gC.nc')['lat' ][...]
lon		= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/without_land_use_change/cesm1bgc_pftcon_gpp_gC.nc')['lon' ][...]


win_len     = 25*12
dates_ar    = time_dim_dates(base_date=dt.date(1850,01,01), total_timestamps=time.size)
start_dates = [dates_ar[i*win_len] for i in range(int(time.size/win_len))] #list of the start dates of all 25 year windows
end_dates   = [dates_ar[i*win_len+win_len-1] for i in range(int(time.size/win_len))] #list of the end dates of all 25 year windows



keys_wrt_wo_lulcc_pos = ['neg_w_lulcc_based_pos_wo_lulcc', 'pos_w_lulcc_based_pos_wo_lulcc','neg_wo_lulcc_based_pos_wo_lulcc','pos_wo_lulcc_based_pos_wo_lulcc']
case_name   = ext_type+'_'+conf+'_based_pos_wo_lulcc'

#binary of 1s for extreme events for 2000-24
bin_tce_1s	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_1s_'+case_name+'.nc')['gpp_TCE_1s'] [0,:,lat_in,lon_in]
bin_tce_01s	= nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_01s_'+case_name+'.nc')['gpp_TCE_01s'] [0,:,lat_in,lon_in]
if ext_type == 'neg':
	case_name_alt   = 'pos'+'_'+conf+'_based_pos_wo_lulcc'
	bin_tce_1s_alt  = nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_1s_'+case_name_alt+'.nc')['gpp_TCE_1s'] [0,:,lat_in,lon_in]
	bin_tce_01s_alt = nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_01s_'+case_name_alt+'.nc')['gpp_TCE_01s'] [0,:,lat_in,lon_in]
if ext_type == 'pos':
	case_name_alt   = 'neg'+'_'+conf+'_based_pos_wo_lulcc'
	bin_tce_1s_alt  = nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_1s_'+case_name_alt+'.nc')['gpp_TCE_1s'] [0,:,lat_in,lon_in]
	bin_tce_01s_alt = nc4.Dataset('/home/ud4/CESM1-BGC/gpp/comparison/pftcon_lulcc/nc_files/bin_TCE_01s_'+case_name_alt+'.nc')['gpp_TCE_01s'] [0,:,lat_in,lon_in]

from    scipy   import  ndimage
larray,narray   = ndimage.label(bin_tce_1s,structure = np.ones(3))
locations       = ndimage.find_objects(larray)

ext_locs = []
for loc in locations:
	start 	= dates_ar[pwin*win_len:(pwin+1)*win_len][loc[0]][0]
	end 	= dates_ar[pwin*win_len:(pwin+1)*win_len][loc[0]][-1]+dt.timedelta(1)
	ext_locs.append((start,end))


drivers_names                = np.array(['Prcp','Soilmoist',r'$P-ET$',r'$T_{max}$',r'$T_{sa}$',r'$T_{min}$','Fire'])
drivers_code                 = np.array([  10,   20,   30,   40,       50,      60,         70])
if th_type  == 'misc' :
	path_in           = '/home/ud4/CESM1-BGC/corr_triggers/misc/' + case_name + '/'
	dom_dri        	  =   nc4.Dataset( path_in + 'dominant_driver_correlation.nc' )

#printing the dominant driver names and correlation coefficients:
# top 3 ranks
# lags from 0-3 months
tmp_dri_id = dom_dri['dri_id'][0:3,pwin,0:4,lat_in,lon_in]
tmp_dri_cc = dom_dri['dri_coeff'][0:3,pwin,0:4,lat_in,lon_in]

print "rank 1 : Most dominant driver"# and lag 0 is ignored"
for i in range(tmp_dri_id.shape[0]):
	for j in range(tmp_dri_id.shape[1])[1:]:
		if tmp_dri_id.mask[i,j] ==False:
			print "rank:",i+1, "  lag:",j,"  driver:" , drivers_names[np.where(drivers_code==tmp_dri_id[i,j])[0]], " Corr Coeff:", "%.2f"%tmp_dri_cc[i,j]


fig,ax  = plt.subplots(nrows=(len(variables_list)+2),ncols= 1,gridspec_kw = {'wspace':0, 'hspace':0.05},tight_layout = True, figsize = (10,7), dpi = 400)
fig.suptitle ('TS of the Drivers')
ax      = ax.ravel()
for idx, var in enumerate (variables_list):
	ax[idx] .plot(dates_ar[pwin*win_len:(pwin+1)*win_len], nc_data['in'][conf][var][pwin*win_len:(pwin+1)*win_len, lat_in, lon_in], label = var)
	print var
	ax[idx].set_xticklabels([])
	ax[idx].tick_params(axis="x",direction="in")
	ax[idx].set_ylabel(" %s"%(var.upper()))
	for (start,end) in ext_locs:
		ax[idx].axvspan(start,end,alpha = .3, color = 'red')

ax[-2].plot(dates_ar[pwin*win_len:(pwin+1)*win_len], nc_gpp[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in], label = 'gpp')
ax[-2].set_ylabel("GPP")
for (start,end) in ext_locs:
	ax[-2].axvspan(start,end,alpha = .3, color = 'red')
ax[-2].set_xticklabels([])
ax[-2].tick_params(axis="x",direction="in")

ax[-1].plot(dates_ar[pwin*win_len:(pwin+1)*win_len], nc_gpp_ano[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in], label = 'gpp')
ax[-1].set_xlabel('Time',fontsize =14)
ax[-1].set_ylabel("GPP_ANO")
for (start,end) in ext_locs:
	ax[-1].axvspan(start,end,alpha = .3, color = 'red')
#ax[0].legend(loc = 'upper center',ncol=len(variables_list)+1, bbox_to_anchor=(.1,1.15),frameon =False,fontsize=9,handletextpad=0.1)
fig.savefig('/home/ud4/CESM1-BGC/corr_triggers/misc/'+case_name+'/shaded_TCE/win_%d.pdf'%int(pwin))

if plt_ano in ['y','Y','yes']:
	fig_ano,ax  = plt.subplots(nrows=(len(variables_list)+2),ncols= 1,gridspec_kw = {'wspace':0, 'hspace':0.05},tight_layout = True, figsize = (10,7), dpi = 400)
	fig_ano.suptitle ('Anomalies of the Driversi: %s'%case_name)
	ax      = ax.ravel()
	for idx, var in enumerate (variables_list):
		ax[idx] .plot(dates_ar[pwin*win_len:(pwin+1)*win_len], nc_data['in'][conf][var+'_ano'][pwin*win_len:(pwin+1)*win_len, lat_in, lon_in], label = var)
		print var
		ax[idx].set_xticklabels([])
		ax[idx].tick_params(axis="x",direction="in")
		ax[idx].set_ylabel(" %s"%(var.upper()))
		for (start,end) in ext_locs:
			ax[idx].axvspan(start,end,alpha = .3, color = 'red')

	ax[-2].plot(dates_ar[pwin*win_len:(pwin+1)*win_len], nc_gpp[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in], label = 'gpp')
	ax[-2].set_ylabel("GPP")
	for (start,end) in ext_locs:
		ax[-2].axvspan(start,end,alpha = .3, color = 'red')
	ax[-2].set_xticklabels([])
	ax[-2].tick_params(axis="x",direction="in")

	ax[-1].plot(dates_ar[pwin*win_len:(pwin+1)*win_len], nc_gpp_ano[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in], label = 'gpp')
	ax[-1].set_xlabel('Time',fontsize =14)
	ax[-1].set_ylabel("GPP_ANO")
	for (start,end) in ext_locs:
		ax[-1].axvspan(start,end,alpha = .3, color = 'red')
	#ax[0].legend(loc = 'upper center',ncol=len(variables_list)+1, bbox_to_anchor=(.1,1.15),frameon =False,fontsize=9,handletextpad=0.1)
	fig_ano.savefig('/home/ud4/CESM1-BGC/corr_triggers/misc/'+case_name+'/shaded_TCE/win_%d_ano.pdf'%int(pwin))

#ploting the location on the map
loc_on_map = np.ma.masked_all((lat.size,lon.size))
loc_on_map[int(lat_in),int(lon_in)] = 1
lat_deg = lat[int(lat_in)]
lon_deg = lon[int(lon_in)]


fig_map,ax = plt.subplots(figsize = (5,2.8),tight_layout=True,dpi=500)
bmap    = Basemap(  projection  =   'eck4',
					lon_0       =   0.,
					resolution  =   'c')
x,y 	= bmap(lon_deg,lat_deg)
bmap	. plot(x,y,'bo',markersize=10)
LAT,LON = np.meshgrid(lat[...], lon[...],indexing ='ij')
ax      = bmap.pcolormesh(LON,LAT,np.ma.masked_invalid(loc_on_map),latlon=True,cmap= cm.afmhot_r)#,vmax= ymax, vmin= ymin)
bmap    .drawparallels(np.arange(-90., 90., 30.),fontsize=14, linewidth = .2)
bmap    .drawmeridians(np.arange(0., 360., 60.),fontsize=14, linewidth = .2)
bmap    .drawcoastlines(linewidth = .2)
plt.title("2000-24; The location!")
fig_map.savefig('/home/ud4/CESM1-BGC/corr_triggers/misc/'+case_name+'/shaded_TCE/win_%d_map.pdf'%int(pwin))
plt.close(fig_map)

# Focusing on one TCE
#--------------------
while True:
	focus_TCE = raw_input('Do you want to plot any of the TCE? (left most : 0 and then +1 towards right) ')
	if focus_TCE in ['n','N','']:
		break
	else:
		focus_TCE = int (focus_TCE)

		tmp_locations = [] # modifying the locations of the extreme events to show drivers upto 12 months prior the TCE and shade will correspond to real time of extreme
		for loc in locations:
			tmp_locations.append(slice(loc[0].start-12, loc[0].stop+2))
		
		fig_focus,ax  = plt.subplots(nrows=(len(variables_list)+2),ncols= 1,gridspec_kw = {'wspace':0, 'hspace':0.05},tight_layout = True, figsize = (10,7), dpi = 400)
		fig_focus.suptitle ('TS of the Drivers Focus TCE')
		ax      = ax.ravel()
		for idx, var in enumerate (variables_list):
			ax[idx] .plot(dates_ar[pwin*win_len:(pwin+1)*win_len][tmp_locations[focus_TCE]], nc_data['in'][conf][var][pwin*win_len:(pwin+1)*win_len, lat_in, lon_in][tmp_locations[focus_TCE]], label = var)
			print var
			ax[idx].set_xticklabels([])
			ax[idx].tick_params(axis="x",direction="in")
			ax[idx].set_ylabel(" %s"%(var.upper()))
			for (start,end) in ext_locs[focus_TCE:focus_TCE+1]:
				ax[idx].axvspan(start,end,alpha = .3, color = 'red')
		ax[-2].plot(dates_ar[pwin*win_len:(pwin+1)*win_len][tmp_locations[focus_TCE]], nc_gpp[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in][tmp_locations[focus_TCE]], label = 'gpp')
		ax[-2].set_ylabel("GPP")
		for (start,end) in ext_locs [focus_TCE:focus_TCE+1]:
			ax[-2].axvspan(start,end,alpha = .3, color = 'red')
		ax[-2].set_xticklabels([])
		ax[-2].tick_params(axis="x",direction="in")

		ax[-1].plot(dates_ar[pwin*win_len:(pwin+1)*win_len][tmp_locations[focus_TCE]], nc_gpp_ano[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in][tmp_locations[focus_TCE]], label = 'gpp')
		ax[-1].set_xlabel('Time',fontsize =14)
		ax[-1].set_ylabel("GPP_ANO")
		for (start,end) in ext_locs[focus_TCE:focus_TCE+1]:
			ax[-1].axvspan(start,end,alpha = .3, color = 'red')
		#ax[0].legend(loc = 'upper center',ncol=len(variables_list)+1, bbox_to_anchor=(.1,1.15),frameon =False,fontsize=9,handletextpad=0.1)
		fig_focus.savefig('/home/ud4/CESM1-BGC/corr_triggers/misc/'+case_name+'/shaded_TCE/win_%d_focus_TCE.pdf'%int(pwin))
	
		# Anomalies of the TCE
		# --------------------

		fig_ano_focus,ax  = plt.subplots(nrows=(len(variables_list)+2),ncols= 1,gridspec_kw = {'wspace':0, 'hspace':0.05},tight_layout = True, figsize = (10,7), dpi = 400)
		fig_ano_focus.suptitle ('Anomalies of the Drivers Focus TCE')
		ax      = ax.ravel()
		for idx, var in enumerate (variables_list):
			ax[idx] .plot(dates_ar[pwin*win_len:(pwin+1)*win_len][tmp_locations[focus_TCE]], nc_data['in'][conf][var+'_ano'][pwin*win_len:(pwin+1)*win_len, lat_in, lon_in][tmp_locations[focus_TCE]], label = var)
			print var ,'_ano'
			ax[idx].set_xticklabels([])
			ax[idx].tick_params(axis="x",direction="in")
			ax[idx].set_ylabel(" %s"%(var.upper()))
			for (start,end) in ext_locs[focus_TCE:focus_TCE+1]:
				ax[idx].axvspan(start,end,alpha = .3, color = 'red')
		ax[-2].plot(dates_ar[pwin*win_len:(pwin+1)*win_len][tmp_locations[focus_TCE]], nc_gpp[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in][tmp_locations[focus_TCE]], label = 'gpp')
		ax[-2].set_ylabel("GPP")
		for (start,end) in ext_locs [focus_TCE:focus_TCE+1]:
			ax[-2].axvspan(start,end,alpha = .3, color = 'red')
		ax[-2].set_xticklabels([])
		ax[-2].tick_params(axis="x",direction="in")

		ax[-1].plot(dates_ar[pwin*win_len:(pwin+1)*win_len][tmp_locations[focus_TCE]], nc_gpp_ano[pwin*win_len:(pwin+1)*win_len, lat_in, lon_in][tmp_locations[focus_TCE]], label = 'gpp')
		ax[-1].set_xlabel('Time',fontsize =14)
		ax[-1].set_ylabel("GPP_ANO")
		for (start,end) in ext_locs[focus_TCE:focus_TCE+1]:
			ax[-1].axvspan(start,end,alpha = .3, color = 'red')
		#ax[0].legend(loc = 'upper center',ncol=len(variables_list)+1, bbox_to_anchor=(.1,1.15),frameon =False,fontsize=9,handletextpad=0.1)
		fig_ano_focus.savefig('/home/ud4/CESM1-BGC/corr_triggers/misc/'+case_name+'/shaded_TCE/win_%d_focus_ano_TCE.pdf'%int(pwin))

