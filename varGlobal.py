#/usr/bin/python
# Global var
import numpy as np
import parameters as param
from decimal import *
import math

# f ~ coriolis matrix
def _init():
	global u
	global v
	global h
	global f
	u = np.zeros([param.nx, param.ny, param.nt], dtype = Decimal, order = 'F')
	v = np.zeros([param.nx, param.ny, param.nt], dtype = Decimal, order = 'F')
	h = np.zeros([param.nx, param.ny, param.nt], dtype = Decimal, order = 'F')
	f = np.zeros([param.nx, param.ny], dtype = Decimal, order = 'F')
	for i1 in range(0,param.ny):
		f[:, i1] = 2.0*param.omega*math.sin(param.pi*param.lat[i1]/param.ly)

def setVar(x,varname,t):
	global u
	global v
	global h
	global f
	if varname == "u":
		u[:,:,t] = x[:,:]
	elif varname == "v":
		v[:,:,t] = x[:,:]
	elif varname == "h":
		h[:,:,t] = x[:,:]
	elif varname == "f":
		f[:,:] = x[:,:]
	else:
		print("ERROR!!!Input varname is not \"u\" or \"v\" or \"h\" or \"f\"!!!Check your code.")

def getVar(varname,t):
	global u
	global v
	global h
	global f

	if varname is "u":
		x = u[:,:,t]
	elif varname is "v":
		x = v[:,:,t]
	elif varname is "h":
		x = h[:,:,t]

	return(x)


