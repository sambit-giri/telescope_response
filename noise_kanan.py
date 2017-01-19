import numpy as np
import c2raytools as c2t
#from telescope_functions import *
from Ctelescope import *

def kelvin_jansky_conversion(ncells, z, boxsize=None):
	if not boxsize: boxsize = c2t.conv.LB
	KB_SI       = 1.38e-23
	janskytowatt= 1e-26
	dist_z      = c2t.z_to_cdist(z)
	boxsize_pp  = boxsize/dist_z				 #in rad	
	omega_pixel = boxsize_pp**2/ncells**2
	omega_total = boxsize_pp**2.0
	c_light_SI  = 2.99792458e+8                              #in m
	mktomujy_nuc= 2.0*KB_SI/c_light_SI/c_light_SI/janskytowatt*((c2t.z_to_nu(z)*1e6)**2.0)*1e3
	con_sol     = mktomujy_nuc*omega_pixel
	return con_sol

def jansky_2_kelvin(array, z, boxsize=None):
	ncells  = array.shape[0]
	con_sol = kelvin_jansky_conversion(ncells, z, boxsize=boxsize)	
	return  array/con_sol

def kelvin_2_jansky(array, z, boxsize=None):
	ncells  = array.shape[0]
	con_sol = kelvin_jansky_conversion(ncells, z, boxsize=boxsize)	
	return  array*con_sol

def noise_map(ncells, z, depth_mhz, obs_time=1000, filename=None, boxsize=None, total_int_time=4., int_time=10., declination=30.):
	"""
	Parameter
	z: 		   Redshift.
	ncells:    The grid size.
	depth_mhz: The bandwidth in MHz.
	obs_time:  The observation time in hours.
	filename:  The path to the file containing the telescope configuration.	
	
	Return
	noise_map: A 2D slice of the interferometric noise at that frequency.
	"""
	if not filename: filename = 'input/lati_long_cor_final.dat'
	uv_map, N_ant  = daily_observation(z, ncells, filename, total_int_time=total_int_time, int_time=int_time, boxsize=boxsize, declination=declination)
	sigma, rms_noi = kanan_noise_image_ska(z, uv_map, depth_mhz, obs_time, N_ant_ska=N_ant)
	noise_real = np.random.normal(loc=0.0, scale=rms_noi, size=(ncells, ncells))
	noise_imag = np.random.normal(loc=0.0, scale=rms_noi, size=(ncells, ncells))
	noise_arr  = noise_real + 1.j*noise_imag
	noise_four = apply_uv_response(noise_arr, uv_map)
	noise_map  = np.fft.ifft2(noise_four)
	return np.real(noise_map)

def telescope_response_on_image(array, z, depth_mhz, obs_time=1000, filename=None, boxsize=None, total_int_time=4., int_time=10., declination=30.):
	if not filename: filename = 'input/lati_long_cor_final.dat'
	uv_map, N_ant  = daily_observation(z, ncells, filename, obs_time, total_int_time=total_int_time, int_time=int_time, boxsize=boxsize, declination=declination)
	img_arr  = np.fft.fft2(array)
	img_four = apply_uv_response(img_arr, uv_map)
	img_map  = np.fft.ifft2(img_four)
	return np.real(img_map)






		

