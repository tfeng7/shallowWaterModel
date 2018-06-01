#/usr/bin/python
# solver module
import parameters as param
import initialization as init
import numpy as np
import varGlobal as var
import pdb
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

def __init__():
	global ku1,kv1,kh1,ku2,kv2,kh2,ku3,kv3,kh3,ku4,kv4,kh4
	nx  = param.nx
	ny  = param.ny
	uo  = np.zeros([nx, ny], dtype = float, order = 'F')
	vo  = np.zeros([nx, ny], dtype = float, order = 'F')
	ho  = np.zeros([nx, ny], dtype = float, order = 'F')
	ku1 = np.zeros([nx, ny], dtype = float, order = 'F')
	kv1 = np.zeros([nx, ny], dtype = float, order = 'F')
	kh1 = np.zeros([nx, ny], dtype = float, order = 'F')
	ku2 = np.zeros([nx, ny], dtype = float, order = 'F')
	kv2 = np.zeros([nx, ny], dtype = float, order = 'F')
	kh2 = np.zeros([nx, ny], dtype = float, order = 'F')
	ku3 = np.zeros([nx, ny], dtype = float, order = 'F')
	kv3 = np.zeros([nx, ny], dtype = float, order = 'F')
	kh3 = np.zeros([nx, ny], dtype = float, order = 'F')
	ku4 = np.zeros([nx, ny], dtype = float, order = 'F')
	kv4 = np.zeros([nx, ny], dtype = float, order = 'F')
	kh4 = np.zeros([nx, ny], dtype = float, order = 'F')

def boundary_periodic(i2,nx,xb,xc,xf):
	# the boundary_periodic can be modified to a matrix instead the loop in swe_solve
	if i2 == 0:
		xb = nx-1
		xc = i2
		xf = i2+1
	elif i2 == nx-1:
		xb = i2-1
		xc = i2
		xf = 0
	else:
		xb = i2-1
		xc = i2
		xf = i2+1

	return xb, xc, xf

def swe_solve(uin, vin, hin, uout, vout, hout):
	# calculate gradient
	global t
	xb = 0
	xc = 0
	xf = 0
	yb = 0
	yc = 0
	yf = 0
	ny = param.ny
	nx = param.nx
	dx = param.dx
	dy = param.dy
	vd = param.vd
	nu = param.nu
	f = var.f
	dx20 = 1.0/(2.0*dx) # here might be an issue in original code
	dy20 = 1.0/(2.0*dy) # should be 1.0/(2.0*dy)
	gdx20 = param.g * dx20
	gdy20 = param.g * dy20
	# diffusion terms
	nudxx = nu/(dx*dx)
	nudyy = nu/(dy*dy)
	ku1 = np.zeros([nx, ny], dtype = float, order = 'F')
	kv1 = np.zeros([nx, ny], dtype = float, order = 'F')
	kh1 = np.zeros([nx, ny], dtype = float, order = 'F')

  # precalculate terms h*u and h*v
	huin = hin * uin
	hvin = hin * vin

  # use numpy.gradient to simplify and boost the gradient calculation
  # use 2-order center finite diff instead of 1-order
  # use Arakawa-C grid
	if param.loopMethod == 1:  # loop
		for j2 in range(0,ny):
			# boundary_periodic(j2,ny,yb,yc,yf)
			yb,yc,yf = boundary_periodic(j2,ny,yb,yc,yf)
			for i2 in range(0,nx):
				#if j2==50 and i2 == 50:
				#	pdb.set_trace()
				xb,xc,xf = boundary_periodic(i2,nx,xb,xc,xf)
				phupx = (huin[xf,yc]-huin[xb,yc])*dx20
				# phupy = (huin[xc,yf]-huin[xc,yb])*dy20
				# phvpx = (hvin[xf,yc]-hvin[xb,yc])*dx20
				phvpy = (hvin[xc,yf]-hvin[xc,yb])*dy20
				upupx = uin[xc,yc] * (uin[xf,yc]-uin[xb,yc])*dx20
				vpupy = vin[xc,yc] * (uin[xc,yf]-uin[xc,yb])*dy20
				upvpx = uin[xc,yc] * (vin[xf,yc]-vin[xb,yc])*dx20
				vpvpy = vin[xc,yc] * (vin[xc,yf]-vin[xc,yb])*dy20
				fv   = f[xc,yc]*vin[xc,yc]
				fu   = f[xc,yc]*uin[xc,yc]
				gphpx = (hin[xf,yc]-hin[xb,yc])*gdx20
				gphpy = (hin[xc,yf]-hin[xc,yb])*gdy20
				bu = vd*uin[xc,yc]
				bv = vd*vin[xc,yc]
				pupxx = (uin[xf,yc]-2.0*uin[xc,yc]+uin[xb,yc])*nudxx
				pupyy = (uin[xc,yf]-2.0*uin[xc,yc]+uin[xc,yb])*nudyy
				pvpxx = (vin[xf,yc]-2.0*vin[xc,yc]+vin[xb,yc])*nudxx
				pvpyy = (vin[xc,yf]-2.0*vin[xc,yc]+vin[xc,yb])*nudyy

				phpt = - (phupx + phvpy)
				pupt = - (upupx + vpupy) + fv - gphpx - bu + (pupxx + pupyy)
				pvpt = - (upvpx + vpvpy) - fu - gphpy - bv + (pvpxx + pvpyy)
				ku1[i2,j2] = pupt
				kv1[i2,j2] = pvpt
				kh1[i2,j2] = phpt

	elif param.loopMethod == 2:  # martix manipulation
				phupx = (huin[2:nx,1:ny-1]-huin[0:(nx-2),1:(ny-1)])*dx20
				phupy = (huin[1:(nx-1),2:ny]-huin[1:(nx-1),0:(ny-2)])*dy20
				phvpx = (hvin[2:nx,1:(ny-1)]-hvin[0:(nx-2),1:(ny-1)])*dx20
				phvpy = (hvin[1:(nx-1),2:ny]-hvin[1:(nx-1),0:(ny-2)])*dy20
				upupx = uin[1:(nx-1),1:(ny-1)] * (uin[2:nx,1:(ny-1)]-uin[0:(nx-2),1:(ny-1)])*dx20
				vpupy = vin[1:(nx-1),1:(ny-1)] * (uin[1:(nx-1),2:ny]-uin[1:(nx-1),0:(ny-2)])*dy20
				upvpx = uin[1:(nx-1),1:(ny-1)] * (vin[2:nx,1:(ny-1)]-vin[0:(nx-2),1:(ny-1)])*dx20
				vpvpy = vin[1:(nx-1),1:(ny-1)] * (vin[1:(nx-1),2:ny]-vin[1:(nx-1),0:(ny-2)])*dy20
				fv   = f[1:(nx-1),1:(ny-1)]*vin[1:(nx-1),1:(ny-1)]
				fu   = f[1:(nx-1),1:(ny-1)]*uin[1:(nx-1),1:(ny-1)]
				gphpx = (hin[2:nx,1:(ny-1)]-hin[0:(nx-2),1:(ny-1)])*gdx20
				gphpy = (hin[1:(nx-1),2:ny]-hin[1:(nx-1),0:(ny-2)])*gdy20
				bu = vd*uin[1:(nx-1),1:(ny-1)]
				bv = vd*vin[1:(nx-1),1:(ny-1)]
				pupxx = (uin[2:nx,1:(ny-1)]-2.0*uin[1:(nx-1),1:(ny-1)]+uin[0:(nx-2),1:(ny-1)])*nudxx
				pupyy = (uin[1:(nx-1),2:ny]-2.0*uin[1:(nx-1),1:(ny-1)]+uin[1:(nx-1),0:(ny-2)])*nudyy
				pvpxx = (vin[2:nx,1:(ny-1)]-2.0*vin[1:(nx-1),1:(ny-1)]+vin[0:(nx-2),1:(ny-1)])*nudxx
				pvpyy = (vin[1:(nx-1),2:ny]-2.0*vin[1:(nx-1),1:(ny-1)]+vin[1:(nx-1),0:(ny-2)])*nudyy

				kh1[1:nx-1,1:ny-1] = - (phupx + phvpy)
				ku1[1:nx-1,1:ny-1] = - (upupx + vpupy) + fv - gphpx - bu + (pupxx + pupyy)
				kv1[1:nx-1,1:ny-1] = - (upvpx + vpvpy) - fu - gphpy - bv + (pvpxx + pvpyy)

				#deal with boundaries
				for j2 in [0,ny-1]:
					yb,yc,yf = boundary_periodic(j2,ny,yb,yc,yf)
					for i2 in range(0,nx):
						# pdb.set_trace()
						xb,xc,xf = boundary_periodic(i2,nx,xb,xc,xf)
						phupx = (huin[xf,yc]-huin[xb,yc])*dx20
						phupy = (huin[xc,yf]-huin[xc,yb])*dy20
						phvpx = (hvin[xf,yc]-hvin[xb,yc])*dx20
						phvpy = (hvin[xc,yf]-hvin[xc,yb])*dy20
						upupx = uin[xc,yc] * (uin[xf,yc]-uin[xb,yc])*dx20
						vpupy = vin[xc,yc] * (uin[xc,yf]-uin[xc,yb])*dy20
						upvpx = uin[xc,yc] * (vin[xf,yc]-vin[xb,yc])*dx20
						vpvpy = vin[xc,yc] * (vin[xc,yf]-vin[xc,yb])*dy20
						fv   = f[xc,yc]*vin[xc,yc]
						fu   = f[xc,yc]*uin[xc,yc]
						gphpx = (hin[xf,yc]-hin[xb,yc])*gdx20
						gphpy = (hin[xc,yf]-hin[xc,yb])*gdy20
						bu = vd*uin[xc,yc]
						bv = vd*vin[xc,yc]
						pupxx = (uin[xf,yc]-2.0*uin[xc,yc]+uin[xb,yc])*nudxx
						pupyy = (uin[xc,yf]-2.0*uin[xc,yc]+uin[xc,yb])*nudyy
						pvpxx = (vin[xf,yc]-2.0*vin[xc,yc]+vin[xb,yc])*nudxx
						pvpyy = (vin[xc,yf]-2.0*vin[xc,yc]+vin[xc,yb])*nudyy

						phpt = - (phupx + phvpy)
						pupt = - (upupx + vpupy) + fv - gphpx - bu + (pupxx + pupyy)
						pvpt = - (upvpx + vpvpy) - fu - gphpy - bv + (pvpxx + pvpyy)
						ku1[i2,j2] = pupt
						kv1[i2,j2] = pvpt
						kh1[i2,j2] = phpt

				for i2 in [0,nx-1]:
					# boundary_periodic(j2,ny,yb,yc,yf)
					xb,xc,xf = boundary_periodic(i2,nx,xb,xc,xf)
					for j2 in range(1,ny-1):
						# if t>30:
						# 	pdb.set_trace()
						yb,yc,yf = boundary_periodic(j2,ny,yb,yc,yf)
						phupx = (huin[xf,yc]-huin[xb,yc])*dx20
						phupy = (huin[xc,yf]-huin[xc,yb])*dy20
						phvpx = (hvin[xf,yc]-hvin[xb,yc])*dx20
						phvpy = (hvin[xc,yf]-hvin[xc,yb])*dy20
						upupx = uin[xc,yc] * (uin[xf,yc]-uin[xb,yc])*dx20
						vpupy = vin[xc,yc] * (uin[xc,yf]-uin[xc,yb])*dy20
						upvpx = uin[xc,yc] * (vin[xf,yc]-vin[xb,yc])*dx20
						vpvpy = vin[xc,yc] * (vin[xc,yf]-vin[xc,yb])*dy20
						fv   = f[xc,yc]*vin[xc,yc]
						fu   = f[xc,yc]*uin[xc,yc]
						gphpx = (hin[xf,yc]-hin[xb,yc])*gdx20
						gphpy = (hin[xc,yf]-hin[xc,yb])*gdy20
						bu = vd*uin[xc,yc]
						bv = vd*vin[xc,yc]
						pupxx = (uin[xf,yc]-2.0*uin[xc,yc]+uin[xb,yc])*nudxx
						pupyy = (uin[xc,yf]-2.0*uin[xc,yc]+uin[xc,yb])*nudyy
						pvpxx = (vin[xf,yc]-2.0*vin[xc,yc]+vin[xb,yc])*nudxx
						pvpyy = (vin[xc,yf]-2.0*vin[xc,yc]+vin[xc,yb])*nudyy

						phpt = - (phupx + phvpy)
						pupt = - (upupx + vpupy) + fv - gphpx - bu + (pupxx + pupyy)
						pvpt = - (upvpx + vpvpy) - fu - gphpy - bv + (pvpxx + pvpyy)
						ku1[i2,j2] = pupt
						kv1[i2,j2] = pvpt
						kh1[i2,j2] = phpt

	#pdb.set_trace()
	if t>30:
		np.set_printoptions(threshold=np.inf)
		#print kh1
	return ku1, kv1, kh1

def swe_rk4(ts):
	'''
	 t timestep
	'''
	global ku1,kv1,kh1,ku2,kv2,kh2,ku3,kv3,kh3,ku4,kv4,kh4
	global t
	t = ts
	nx = param.nx
	ny = param.ny


	dt = param.dt
	dt2 = dt / 2.0
	dt6 = dt / 6.0

	uo = var.getVar("u",ts-1)
	vo = var.getVar("v",ts-1)
	ho = var.getVar("h",ts-1)

  # rk4 solve
	ku1, kv1, kh1 = swe_solve(uo        ,         vo,         ho,ku1, kv1, kh1 )
	ku2, kv2, kh2 = swe_solve(uo+ku1*dt2, vo+kv1*dt2, ho+kh1*dt2,ku2, kv2, kh2 )
	ku3, kv3, kh3 = swe_solve(uo+ku2*dt2, vo+kv2*dt2, ho+kh2*dt2,ku3, kv3, kh3 )
	ku4, kv4, kh4 = swe_solve(uo+ku3*dt , vo+kv3*dt , ho+kh3*dt ,ku4, kv4, kh4 )

	un = uo + (ku1 + 2.0*ku2 + 2.0*ku3 + ku4) * dt6
	vn = vo + (kv1 + 2.0*kv2 + 2.0*kv3 + kv4) * dt6
	hn = ho + (kh1 + 2.0*kh2 + 2.0*kh3 + kh4) * dt6

	# u[:,:,ts+1] = un
	# v[:,:,ts+1] = vn
	# h[:,:,ts+1] = hn
	var.setVar(un,"u",ts)
	var.setVar(vn,"v",ts)
	var.setVar(hn,"h",ts)

	if param.iplot == 1:
		plt.clf()
		fig = plt.figure()
		#ax = fig.gca(projection='3d')
		#surf = ax.plot_surface(param.lon, param.lat, ku1)
		plt.contourf(param.lon,param.lat,hn, cmap='coolwarm')
		plt.colorbar()
		plt.pause(1e-5)
		plt.close(fig)

	print("U min/max = ", np.amin(var.u[:,:,ts]), "/",np.amax(var.u[:,:,ts]))

