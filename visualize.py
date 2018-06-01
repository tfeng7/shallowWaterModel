# Program to visualize the data

import numpy as np
from netCDF4 import Dataset as dt
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

fname = 'data_small_rotation.nc'
fp = dt(fname, 'r')

# print(fp.variables)

x = np.array(fp.variables['lon'][:])
y = np.array(fp.variables['lat'][:])
t = np.array(fp.variables['time'][:])
h = np.array(fp.variables['height'][:])
u = np.array(fp.variables['uwnd'][:])
v = np.array(fp.variables['vwnd'][:])

uv = (u**2.0 + v**2.0)**0.5
h = np.swapaxes(h,0,2)
uv = np.swapaxes(uv,0,2)

# regridding data for smoother values

X, Y = np.meshgrid(x,y)

# X, Y = np.mgrid[0:10:1001j, 0:10:1001j]
# H = np.zeros((10001,1001,1001))
'''
for i in range(0, 10001):
    H[i,:,:] = griddata(x, y, h[i,:,:], (X, Y), method='linear')
'''

fig = plt.figure()


ax = fig.gca(projection='3d')

##plot_args = {'rstride': 1, 'cstride': 1, 'facecolor':
##             'blue', 'linewidth': 0.01, 'antialiased': True,
##             'shade': True}

plot_args = {'rstride': 1, 'cstride': 1, 'cmap':
             'Blues', 'linewidth': 0.01, 'antialiased': True,
             'shade': True}

for i in range(0, 10001, 10):
    ax.clear()
    surf = ax.plot_surface(X, Y, h[i,:,:], **plot_args)
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    ax.set_zlabel('h')
    ax.set_zlim([0.99, 1.01])
    plt.pause(1e-5)


#matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

#vmax = np.amax(h)
#vmin = np.amin(h)

#contours = np.linspace(vmin, vmax, 11)
# for i in range(0, 10001, 10):
#    plt.clf()
#    # plt.contourf(x,y,h[i,:,:], cmap='coolwarm')
#    plt.contourf(x,y,h[i,:,:], cmap='coolwarm')
#    plt.xlabel('x axis')
#    plt.ylabel('y axis')
#    plt.suptitle('t = {} s'.format(i*0.1))
#    # plt.colorbar(ticks=contours)
#    plt.colorbar()
#    #plt.contour(x,y,h[i,:,:], colors='black')
#    plt.contour(y,x,uv[i,:,:], colors='black')
#    plt.title('Wave motion')
#    plt.pause(1e-3)

