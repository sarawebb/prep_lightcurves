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

field = 'Dusty11'
year = '2018'
month = '06'

filesave_path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/source_lists/' + field + '/' + year + '/' + month + '/'

if not os.path.exists(filesave_path): 
	os.makedirs(filesave_path, 0o755)

field_source_RA = []
field_source_DEC = []
field_source_g_mag_AUTO= [] 
field_source_g_mag_err_AUTO = [] 
field_source_g_mag_APER = []
field_source_g_mag_err_APER = []

sorted_image_path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/'+ year +'_' + month + '_' +field +'_sorted.ascii'
print(sorted_image_path)
image_names = np.loadtxt(sorted_image_path, unpack = True, skiprows = 1, usecols =[0], dtype = str )
image_mjds = np.loadtxt(sorted_image_path, unpack =True, skiprows= 1, usecols = [1])

start_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year +'/' + month +'/'+ field +'/g_band/single/'
print(start_path)

sorted_table = Table()
sorted_table['filenames'] = image_names
sorted_table['mjds'] = image_mjds 

print(sorted_table[:1])



im_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year +'/' + month +'/'+ field +'/g_band/single/*'

im_paths = glob.glob(im_path)

for row in sorted_table: 
	image_path = str(start_path + row['filenames']) 
	for i in im_paths: 
		#print(i)
		if i == image_path[:len(image_path)-5]  :
			fuckingwork = os.path.join(i ,'final_source_cats')
			#print(fuckingwork)
			#print(fuckingwork)
			#print('yass')
			#os.getcwd()
			os.chdir(fuckingwork)
			
			#print(filename)
			image_info_ra = []
			image_info_dec = []
			image_info_g_mag_auto = []
			image_info_g_mag_auto_err = []
			image_info_g_mag_aper = []
			image_info_g_mag_aper_err = []
			image_info_mjd = []
			for filename in os.listdir('.'):
			
				
			 	
				if filename.endswith('cat_CORRECTED.ascii'): 
					cwd = os.getcwd()
					cat_to_open = os.path.join(cwd, filename)
					#print(cat_to_open)
					
					
					RA, DEC, g_mag_auto, g_mag_auto_err, g_mag_aper, g_mag_aper_error = np.loadtxt( cat_to_open, unpack =True, skiprows = 1)
					
					
					for k in RA: 
						image_info_ra.append(k)
					
					for k in DEC: 
						image_info_dec.append(k)	
					
						
					for k in g_mag_auto: 
						image_info_g_mag_auto.append(k)
					
					for k in g_mag_auto_err: 
						image_info_g_mag_auto_err.append(k)
						
					for k in g_mag_aper: 
						image_info_g_mag_aper.append(k)
				
					for k in g_mag_aper_error: 
						image_info_g_mag_aper_err.append(k)
			

			for k in range(len(image_info_ra)):
				image_info_mjd.append(row['mjds'])

			info_table = Table()
			info_table['RA'] = image_info_ra
			info_table['DEC'] = image_info_dec
			info_table['g_mag_auto'] = image_info_g_mag_auto
			info_table['g_mag_auto_err'] = image_info_g_mag_auto_err
			info_table['g_mag_aper'] = image_info_g_mag_aper
			info_table['g_mag_aper_err'] = image_info_g_mag_aper_err
			info_table['MJD'] = image_info_mjd
			
			
			
			print(info_table)
			
			
			output = filesave_path + row['filenames'] + '_'+ year +'_'+ month + '_'+ field + '.ascii'
			print(output)
			info_table.write(output, format = 'ascii', overwrite =True) 
		
			
