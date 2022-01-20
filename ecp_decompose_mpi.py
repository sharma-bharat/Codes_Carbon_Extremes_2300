""" Time series decomposition using parallel computing using the SSA approach usign the ssa function in  ssa.py"""

import ssa 	        as ssa
import numpy 		as np
import pandas 		as pd
import glob,os
import netCDF4 		as nc4
from mpi4py 		import MPI
import sys
import collections
import datetime 	as dt
from calendar 		import monthrange
#Calling my functions
sys.path.append('/home/ud4/repos/Codes_Carbon_Extremes_2300/')
from functions 		import time_dim_dates, index_and_dates_slicing, mpi_local_and_global_index, geo_idx
import argparse

#Inputs:
#-------
parser  = argparse.ArgumentParser()
parser.add_argument('--variable',   	'-var'	,   help = "Variable to decompose" , type= str, default= 'n')
parser.add_argument('--model_config',   '-config', 	help = "Model Configuration"   , type= str, default= 'w_lulcc')
parser.add_argument('--var_in_gC',		'-var_gC',	help = "If the variable is gC" , type= str, default= 'n')
args = parser.parse_args()

# enter the variable name and the corresponding file name that you want to decompose
variable= args.variable
if args.model_config == 'wo_lulcc': 	#without landuse and landcover change
	in_path = '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
	if args.var_in_gC =='n':
		ncin	=	nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'.nc')
	else: 
		ncin   	=   nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'_gC.nc')
	ncin_area   =   nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'.nc')
elif args.model_config == 'w_lulcc': 		#with landuse and land cover change
	in_path = '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'
	if args.var_in_gC =='n':
		ncin    =   nc4.Dataset(in_path+'cesm1bgc_with_lulcc_'+variable+'.nc')
	else:
		ncin    =   nc4.Dataset(in_path+'cesm1bgc_with_lulcc_'+variable+'_gC.nc')
	ncin_area   =   nc4.Dataset(in_path+'cesm1bgc_with_lulcc_'+variable+'.nc')	

#filename= 	glob.glob(in_path+"*gpp.nc")

ga_names=   ncin.ncattrs() 
time    =   ncin.variables['time']
lat     =   ncin.variables['lat']
lon     =   ncin.variables['lon']
#nc_gpp	=	ncin.variables['GPP']
nc_var	=	ncin.variables[variable]
nc_area	=	ncin_area.variables['area']

#SSA tool
window	=	120 	#for trend removal of 10 years and above

# Parallel Computing task distribution
comm		=	MPI.COMM_WORLD
size		=	comm.Get_size() #Number of processors I am asking for
local_rank	=	comm.Get_rank() #Rank of the current processor
#chunk size or delta
local_n		=	lat.size/size
#calculating the range for every parallel process
begin_idx	=	local_rank*local_n
end_idx		=	begin_idx+local_n

resids		=	np.zeros((time.size,local_n,lon.size))
lats_pp		=	lat[...][begin_idx:end_idx]		#lons per processor
loc_idx_n_lats_pp	=	mpi_local_and_global_index(begin_idx = begin_idx,local_n=local_n) #local idx and global lon index per processor

for k,i in loc_idx_n_lats_pp:
	for j in range(lon.size):
		tmp_ts	=	nc_var[:,i,j]
		#counter	=	collections.Counter(tmp_ts)	

		if nc_area[...].mask[i,j]	== 	True: 			#where ocean/masked print np.nan array
			resids[:,k,j] 			= 	np.array([np.nan]*time.size)
			#print ("1")
		elif collections.Counter(tmp_ts)[0]	== 	len(tmp_ts): 	#where values are all zero print all zeros too
			resids[:,k,j]			=	np.array([0]*time.size)
			#print ("2")
		else:											#else run the decomposition 
			y						= 	np.array(tmp_ts)
			resids[:,k,j]			= 	ssa.GetResidual(y,window)	
			#print ("3")
	print k,i
	#with nc4.Dataset(in_path+'temp3_copy/cesm1gbc_pftcon_gpp_anomalies_'+format(local_rank,'03')+'.nc',mode="w") as dset:
	if  args.model_config == 'wo_lulcc':
		if args.var_in_gC =='n':
			out_filename = in_path+'temp/cesm1bgc_pftcon_'+variable+'_anomalies_'+format(local_rank,'03')+'.nc'
		else: 
			out_filename = in_path+'temp/cesm1bgc_pftcon_'+variable+'_anomalies_gC_'+format(local_rank,'03')+'.nc'
	elif  args.model_config == 'w_lulcc':	
		if args.var_in_gC =='n':
			out_filename = in_path+'temp/cesm1bgc_with_lulcc_'+variable+'_anomalies_'+format(local_rank,'03')+'.nc'
		else:
			out_filename = in_path+'temp/cesm1bgc_with_lulcc_'+variable+'_anomalies_gC_'+format(local_rank,'03')+'.nc'
	
	with nc4.Dataset(out_filename,mode="w") as dset:
		#for ga_name in  ga_names:
			#gval=   ncin.__dict__[ga_name]
			#dset    .setncattr(ga_name,gval)
		dset        .createDimension( "time",size = time.size)
		dset        .createDimension( "lat"	,size = local_n)
		dset        .createDimension( "lon"	,size = lon.size) 
		t  	=	dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time") ,fill_value = 1.e+36)
		y   =	dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat")  ,fill_value = 1.e+36)
		x   =	dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon")  ,fill_value = 1.e+36)
		#z   =	dset.createVariable(varname = "gpp"  ,datatype = float, dimensions = ("time","lat","lon"),fill_value=np.nan)
		z   =	dset.createVariable(varname = variable  ,datatype = float, dimensions = ("time","lat","lon"),fill_value=1.e+36)
		t.axis  =   "T"
		x.axis  =   "X"
		y.axis  =   "Y"
		t[...]  =  	time[...]
		x[...]  =   lon[...]
		y[...]  =   lat[i]
		z[...]  =   resids
		z.missing_value =   1.e+36
		z.stardard_name =   variable+" anomalies using the ssa"
		z.setncattr			("cell_methods",nc_var.cell_methods)
		z.units 		=   ncin.variables[variable].units
		x.units			=	lat.units
		x.missing_value	=	1.e+36
		x.setncattr			("long_name",lat.long_name)
		y.units			=	lon.units
		y.missing_value =   1.e+36
		y.setncattr         ("long_name",lon.long_name)
		t.units			= 	time.units
		t.setncattr			("calendar",time.calendar)
		t.setncattr			("long_name",time.long_name)

