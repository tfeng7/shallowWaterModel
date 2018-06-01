#/usr/bin/python
# initialze the variables for the model
import parameters as param
import numpy as np
import math
import varGlobal as var

def init():
  # meshgrid
	lon2d, lat2d =  np.meshgrid (param.lon,param.lat)
	f = lon2d

	var._init()
	h = var.getVar("h",0)
	np.set_printoptions(threshold='nan')
	h[:,:] = param.h0 - param.perturb * np.exp(-(lat2d**2.0+lon2d**2.0))
	print("max h = ",np.amax(h))
	print("min h = ",np.amin(h))

	# cal coriolis parameter
	for i1 in range(0,param.ny):
	    f[:, i1] = 2.0*param.omega*math.sin(param.pi*param.lat[i1]/param.ly)

	var.setVar(h,"h",0)
	var.setVar(f,"f",0)

	print("Initialization finished.")
