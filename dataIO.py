#/usr/bin/python
# data input and output
import netCDF4 as nc
import parameters as param
from swe_rk4 import swe_rk4 as swe
import varGlobal as var
import time

def output_netcdf():
	da=nc.Dataset(param.fname,'w',format='NETCDF4')
	da.createDimension('time',param.nt)
	da.createDimension('lon',param.nx)
	da.createDimension('lat',param.ny)
	da.createVariable("time",'f',("time"))
	da.createVariable("lon",'f',("lon"))
	da.createVariable("lat",'f',("lat"))
	da.variables['lat'][:]=param.lat
	da.variables['lon'][:]=param.lon
	da.variables['time'][:]=range(0,param.nt)

	da.createVariable('height','f8',('time','lat','lon'))
	da.createVariable('uwnd','f8',('time','lat','lon'))
	da.createVariable('vwnd','f8',('time','lat','lon'))

  # change dimension from Lon-Lat-Time to Time-Lat-Lon
	h = np.swapaxes(var.h,0,2)
  u = np.swapaxes(var.u,0,2)
  v = np.swapaxes(var.v,0,2)
	da.variables['height'][:] = var.h
	da.variables['uwnd'][:] = var.u
	da.variables['vwnd'][:] = var.v
	da.close()


