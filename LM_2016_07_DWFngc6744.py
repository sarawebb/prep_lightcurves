import numpy as np
import matplotlib as mpl
from astropy.table import Table, Column, join 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
import math
import os
import glob
import sys
from sortedcontainers import SortedDict
import datetime as dt
import imageio
import os
from PIL import Image
from matplotlib.colors import LogNorm
from astropy.nddata.utils import Cutout2D
from astropy import units as u
import datetime as dt 
import astropy.units as u
from astroML.crossmatch import crossmatch_angular 
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import glob
import csv
import scipy as sp
from scipy import stats
from scipy.optimize import curve_fit
from math import isclose

field = 'DWF_ngc6744_slew'
year = '2016'
month = '07'
sorted_image_path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/'+year+'_'+month+'_'+field+'_sorted.ascii'

source_lists = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/source_lists/'+ field +'/'+year+'/'+ month+'/'

start_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year+'/'+month+'/'+ field+'/g_band/single/'

image_names = np.loadtxt(sorted_image_path, unpack = True, skiprows = 1, usecols =[0], dtype = str )
image_mjds = np.loadtxt(sorted_image_path, unpack =True, skiprows= 1, usecols = [1])

sorted_table = Table()
sorted_table['filenames'] = image_names
sorted_table['mjds'] = image_mjds 
print(sorted_table)

im_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year +'/' + month +'/'+ field +'/g_band/single/*'

im_paths = glob.glob(im_path)

def func(x, a, b, c):
	return np.float64(a) * np.exp(np.float64(b) * np.float64(x)) + np.float64(c)

	
def powerlaw(x, a, b, c):
	return c + x**a * b
	
filenames = []
limiting_mags = []
	
for row in sorted_table: 
	image_path = str(start_path + row['filenames']) 
	for i in im_paths: 
		#print(i)
		if i == image_path[:len(image_path)-5]:
			fuckingwork = os.path.join(i ,'final_source_cats')
			os.chdir(fuckingwork)
			for filename in os.listdir('.'):
				if filename.endswith('cat_CORRECTED.ascii'):
					print(filename)
					mags = []
					errors = []
					filenames.append(filename)
					cwd = os.getcwd()
					cat_to_open = os.path.join(cwd, filename)
					RA, DEC, g_mag_auto, g_mag_auto_err, g_mag_aper, g_mag_aper_error = np.loadtxt( cat_to_open, unpack =True, skiprows = 1)
					
					#print(len(RA))
					
					for k, j  in zip(g_mag_auto, g_mag_auto_err): 
						if k <= 25: 
							mags.append(k)
							errors.append(j) 
					
					
					#plt.plot(mags, errors, 'r.')
					
					
					xx = mags 
					yy = errors 
					
					new_binned_x = []
					new_binned_y = [] 
					new_min_x = []
					new_min_y = []
					
					for i in np.arange(12, 22, 0.1 ): 
						tempx = []
						tempy = []
						for x, y in zip(xx, yy): 
							if i - 0.1 <= x <= i+0.1: 
								tempx.append(x)
								tempy.append(y)
								
						median_x, blah = sp.stats.mode(tempx)
						median_y, blah = sp.stats.mode(tempy)
						
						
						if 12 <= median_x <= 22: 
							new_binned_x.append(median_x[0])
							new_binned_y.append(median_y[0])
						else: 
							pass
							
													
					#plt.plot(xx, yy, 'b.')	
					plt.plot(xx, yy, 'g.')
					plt.plot(new_binned_x, new_binned_y, 'r*')
					#print(new_binned_x)
					#print(new_binned_y)
					plt.xlabel('g mag')
					plt.ylabel('g mag error')		
					plt.title(filename)					
					
					
					lim_save_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+year+'/'+month+'/'+field +'/g_band/single/' + filename[:26] + '/limiting_mags/'
					plot_save_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+year+'/'+month+'/'+field +'/g_band/single/' + filename[:26] + '/limiting_mags/check_plots/'
					print(lim_save_path)
					if not os.path.exists(lim_save_path):
						os.makedirs(lim_save_path, 0o755)
					else: 
						pass 
						
					if not os.path.exists(plot_save_path):
						os.makedirs(plot_save_path)
					else:
						pass
					
					
					
					
					#yn = new_binned_y + 0.2*np.random.normal(size=len(new_binned_x))
					#popt, pcov = curve_fit(func, new_binned_x, new_binned_y, maxfev=100000)
					
					yn = new_binned_y + 0.2*np.random.normal(size=len(new_binned_x))
					popt, pcov = curve_fit(powerlaw, new_binned_x, new_binned_y, maxfev=100000)
					
					
					
					testing = powerlaw(xx, *popt) # exponational 
					
					#p = np.polyfit(xx, yy, 1) #polynomal 
					
					#ext_line = func(xx, *popt)
					#print(len(ext_line))
					#print(len(xx))
					
					x = np.array(xx)
					
					
					#sp.stats.expon.fit(new_binned_x, new_binned_y)
					
					test = np.sort(testing)
					test2 = np.sort(xx)
					plt.plot(test2, powerlaw(test2, *popt), 'b-')
					plt.ylim(-0.1, 0.4)
					plt.savefig(plot_save_path + filename + '_magerr_vs_mag_plot.png', overwrite = True)
					#plt.show()
					
					a, b, c = popt 
					#print(a, b, c)
					limiting = ((0.198 - c)/(b))**(1/(a))
					print(limiting)
					
					#print(limiting)
					limiting_mags.append(limiting)
					
					'''
					counts, bins, bars = plt.hist(mags, bins = 30) 
					
					x = bins[3:len(bins)-5]
					y = counts[3:len(counts)-4]
					
					xx = bins[:len(bins)-1]
					yy = counts 
					
					yn = y + 0.2*np.random.normal(size=len(x))
					
					
					x = np.array(x)
					popt, pcov = curve_fit(func, x, yn, maxfev=100000)
					
					
					#plt.plot(x, yn, 'ko', label ='work god darn it')
					#plt.plot(x, func(x, *popt), 'r-', label = 'WORK MOTHERFUCKER')
					#plt.legend 
					
					#plt.yscale('log', nonposy='clip')
					
					
					ext_line = func(bins, *popt)
					
					xp = np.linspace(min(bins), max(bins), 100)
					p = np.poly1d(np.polyfit(xx, yy, 30))
					plt.plot( xp, p(xp), '-')
					plt.ylim(0, max(counts))
					plt.show()
					#plt.plot(bins,ext_line, 'b-')
					#plt.close()
					#plt.plot(bins[:len(bins)-1],counts, 'r.')
					#plt.yscale('log', nonposy='clip')
					print(max(counts))
					
					
					for i, j, k in zip(counts[len(counts)-5:], ext_line[len(ext_line)-5:], bins[len(bins)-5:]):
						diff = []
						a = 0.5*j
						if isclose(a, i, abs_tol=30) == True: 
							print(i, a, k)
							b = a -i 
							diff.append(b)
							limiting_mags.append(k)
						else: 
							pass 
						#else: 
							#print('NOTHING FOUND')
							#a = 99
							#limiting_mags.append(a)
					print('HERE')
					#print(limiting_mags)
					#print(len(limiting_mags))
					#print(len(filenames))
			
				'''
#print(filenames)
#print(limiting_mags)
info_table = Table()
info_table['filenames']	= filenames
info_table['mags'] = limiting_mags
#print(info_table['filenames'])
output = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/field_limiting_mags/' + field + '_' + year + '_' + month + '_limitingmags.ascii'
	
info_table.write(output, format = 'ascii', overwrite = True)
						
