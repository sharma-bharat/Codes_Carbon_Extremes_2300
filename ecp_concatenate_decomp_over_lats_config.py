#python 2
#whereever the variable var is used, it is equal to the new varibale that is entered in the command line while running
import netCDF4 as nc4
import numpy as np
import sys
import argparse

#Inputs:
#-------
parser  = argparse.ArgumentParser()
parser.add_argument('--variable',       '-var'  ,    help = "Variable to decompose" , type= str, default= 'n')
parser.add_argument('--model_config',   '-mod_conf', help = "Model Configuration"   , type= str, default= 'w_lulcc')
parser.add_argument('--var_in_gC',      '-var_gC',   help = "If the variable is gC" , type= str, default= 'n')
args = parser.parse_args()

#concatenates all the files for the entered variable
variable= args.variable

if args.model_config == 'wo_lulcc': 	#without landuse and landcover change
	in_path = '/home/ud4/CESM1-BGC/'+variable+'/without_land_use_change/'
	if args.var_in_gC =='n':
		nc_read	= nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'.nc')
		text	= "temp/cesm1bgc_pftcon_"+variable+"_anomalies_"
	else: 
		nc_read   	= nc4.Dataset(in_path+'cesm1bgc_pftcon_'+variable+'_gC.nc')
		text	= "temp/cesm1bgc_pftcon_"+variable+"_anomalies_gC_"
elif args.model_config == 'w_lulcc': 		#with landuse and land cover change
	in_path = '/home/ud4/CESM1-BGC/'+variable+'/with_land_use_change/'
	if args.var_in_gC =='n':
		nc_read    =   nc4.Dataset(in_path+'cesm1bgc_with_lulcc_'+variable+'.nc')
		text    = "temp/cesm1bgc_with_lulcc_"+variable+"_anomalies_"
	else:
		nc_read    =   nc4.Dataset(in_path+'cesm1bgc_with_lulcc_'+variable+'_gC.nc')
		text    = "temp/cesm1bgc_with_lulcc_"+variable+"_anomalies_gC_"


ll		= [format(i,'003') for i in range(192)]
filenames = [text+ll[i]+'.nc' for i in range(len(ll))]

anomalies = np.zeros((5412,192,288))
for i in range(len(filenames)):
	print (in_path+filenames[i])
	anomalies[:,i,:] = nc4.Dataset(in_path+filenames[i],'r').variables[variable][:,0,:]
"""
#using the mask
lf_2d	= nc_read.variables['landfrac'][...]
lf_3d	= np.ma.masked_all(anomalies.shape)
for i in range(lf_3d.shape[0]):
	lf_3d[i] = lf_2f
anomalies	= np.ma.masked_array(anomalies, mask = lf_3d.mask)
"""

ga_names  	= nc_read.ncattrs()
var			= nc_read.variables[variable]
lat			= nc_read.variables['lat']
lon			= nc_read.variables['lon']
time		= nc_read.variables['time']

#saving with mask
anomalies	= np.ma.masked_array(anomalies, mask = var[...].mask)
#Saving the NC file	
if  args.model_config == 'wo_lulcc':
	if args.var_in_gC =='n':
		out_filename = in_path+'cesm1bgc_pftcon_'+variable+'_anomalies.nc'
	else: 
		out_filename = in_path+'cesm1bgc_pftcon_'+variable+'_anomalies_gC.nc'
elif  args.model_config == 'w_lulcc':	
	if args.var_in_gC =='n':
		out_filename = in_path+'cesm1bgc_with_lulcc_'+variable+'_anomalies.nc'
	else:
		out_filename = in_path+'cesm1bgc_with_lulcc_'+variable+'_anomalies_gC.nc'


with nc4.Dataset(out_filename,mode="w") as dset:
	for ga_name in  ga_names:
	    gval=   nc_read.__dict__[ga_name]
	    dset    .setncattr(ga_name,gval)
	dset        .createDimension( "time",size = time.size)
	dset        .createDimension( "lat"	,size = lat.size)
	dset        .createDimension( "lon"	,size = lon.size) 
	t  	=	dset.createVariable(varname = "time" ,datatype = float, dimensions = ("time") ,fill_value = 1.e+36)
	y   =	dset.createVariable(varname = "lat"  ,datatype = float, dimensions = ("lat")  ,fill_value = 1.e+36)
	x   =	dset.createVariable(varname = "lon"  ,datatype = float, dimensions = ("lon")  ,fill_value = 1.e+36)
	z   =	dset.createVariable(varname = variable  ,datatype = float, dimensions = ("time","lat","lon"),fill_value= 1.e+36)
	t.axis  =   "T"
	x.axis  =   "X"
	y.axis  =   "Y"
	t[...]  =  	time[...]
	x[...]  =   lon[...]
	y[...]  =   lat[...]
	z[...]  =   anomalies
	z.missing_value =   1.e+36
	z.standard_name =   variable+" anomalies using the ssa"
	z.setncattr			("cell_methods",var.cell_methods)
	z.units 		=   var.units
	x.units			=	lat.units
	x.missing_value	=	1.e+36
	x.setncattr			("long_name",lat.long_name)
	y.units			=	lon.units
	y.missing_value =   np.nan
	y.setncattr         ("long_name",lon.long_name)
	t.units			= 	time.units
	t.setncattr			("calendar",time.calendar)
	t.setncattr			("long_name",time.long_name)

