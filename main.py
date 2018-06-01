#!/usr/bin/python
# Study to make a shallow water model following a video of youtube.
# Video: https://www.youtube.com/watch?v=7_kcM8Yfv3A
import parameters as param
import initialization as init
import dataIO as dio
import swe_rk4
import varGlobal as var

def main():
	'''
	i: Iterating variable
	nt: Maximum timesteps
	'''
	print("Reading parameters and variables")

	init.init()
	print("Initialization finished. Start calculation.")

	#print(var.h[:,:,0])
	for i in range(1,param.nt):
		print("Integrating for timesstep t = "+str(i)+" of "+str(param.nt-1))
		swe_rk4.__init__()
		swe_rk4.swe_rk4(i)

	dio.output_netcdf()

if __name__ == '__main__':
	main()
