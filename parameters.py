#/usr/bin/python
# static parameters
import numpy as np

nx = 101
ny = 101

xstart = -5 # meter
xend = 5
ystart = -5
yend = 5

dt = 0.0003 # dt < 1/300 *min(dx,dy) following the CFL criterion
tstart = 0.
tend = 100.


iplot = 1        # 1: plot in every timestep; 2: off
loopMethod = 2   # 1: integrate using loop (slow) 2: integrate using matrix arithmetic
perturb = 0.01   # 1: initial pertubation factor
fname = "data_small_rotation.nc" # output filename

# static parameters,!!!DO NOT CHANGE!!!if you don't know the exactly meaning
lx = xend - xstart
ly = yend - ystart
lt = tend - tstart
dx = lx/float(nx-1)
dy = ly/float(ny-1)
nt = lt/float(dt-1) + 1

g = 9.8065
pi = 4.0*tan(1.0)
h0 = 1.0

omega = 1e-5 # rotation parameter earth~7.2921 × 10-5 rad s-1
vd = 1e-3    # viscous drag coeffient of atmosphere ~ 0.005 (over water), 0.015 (over grass) , ocean over the seabed~0.0025
nu = 1.5e-5  # kinematic viscosity of the atmosphere ~ 1.460x10^-5 m^2 s-1 , Ocean~ 1.0 × 10−6 m^2·s^-1

lat =  np.arange(ystart, yend+dy/2, dy)
lon =  np.arange(xstart, xend+dx/2, dx)