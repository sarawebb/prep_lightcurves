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


field = '4hr'
year = '2015'
month = '01'

source_list = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/total_field_sources/fields/' + year + '_' + month +'_' + field + '_full_field_sources.ascii'

source_ra = np.loadtxt(source_list, unpack = True, skiprows = 1, usecols =[0] )
source_dec = np.loadtxt(source_list, unpack =True, skiprows= 1, usecols = [1])

source_table = Table()
source_table['RA'] = source_ra 
source_table['DEC'] = source_dec
print(source_table) 

start_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year+'/'+month+'/'+ field+'/g_band/single/'
sorted_image_path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/'+year+'_'+month+'_'+field+'_sorted.ascii'

image_names = np.loadtxt(sorted_image_path, unpack = True, skiprows = 1, usecols =[0], dtype = str )
image_mjds = np.loadtxt(sorted_image_path, unpack =True, skiprows= 1, usecols = [1])

sorted_table = Table()
sorted_table['filenames'] = image_names
sorted_table['mjds'] = image_mjds 
print(sorted_table)

im_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year +'/' + month +'/'+ field +'/g_band/single/*'

im_paths = glob.glob(im_path)

for row in source_table: 
	SOR_X = np.empty((1, 2), dtype = np.float64)
	SOR_X[:,0] = row['RA']
	SOR_X[:,1] = row['DEC']
	print(SOR_X)
	
	mags = []
	magerrs = []
	filenames =[]
	
	for row1 in sorted_table:
		image_path = str(start_path + row1['filenames'])
		for i in im_paths: 
			if i == image_path[:len(image_path)-5]:
				enter_cats = os.path.join(i, 'final_source_cats')
				os.chdir(enter_cats)
				filenames.append(row1['filenames']
				
				checker = []
				for filename in os.listdir('.'):
					if filename.endswith('cat_CORRECTED.ascii'):
						cwd = os.getcwd()
						cat_to_open = os.path.join(cwd, filename)
						RA, DEC, g_mag_auto, g_mag_auto_err, g_mag_aper, g_mag_aper_error = np.loadtxt( cat_to_open, unpack =True, skiprows = 1)
						
						CAT_X = np.empty((len(RA), 4), dtype = np.float64)
						CAT_X[:, 0] = RA 
						CAT_X[:, 1] = DEC
						CAT_X[:, 2] = g_mag_auto
						CAT_X[:, 3] = g_mag_auto_err 
						
						max_radius = 1/3600 #.5 arc second
						dist_between, ind_row = crossmatch_angular(SOR_X, CAT_X, max_radius)
						match = ~np.isinf(dist_between)
						
						if len(match) != 0: 
							match_table = Table()
							match_table['matched_true_false'] = match 
							match_table['matched_ID'] = ind_row
							match_table['matched_source_RA'] = row['RA']
							match_table['match_source_DEC'] = row['DEC']
							
							cat_match_true = []
							cat_row_match = []
							
							for row2 in match_table: 
								if row2['matched_true_false'] == True: 
									cat_match_true.append(row2['matched_true_false'])
									cat_row_match.append(row2['matched_ID']) 
									
							cat_ra = []
							cat_dec = []
							cat_g_band = []
							cat_g_err = []
							
							for j in cat_row_match: 
								cat_ra.append(CAT_X[j, 0])
								cat_dec.append(CAT_X[j, 1])
								cat_g_band.append(CAT_X[j, 2])
								cat_g_err.append(CAT_X[j, 3])
							
							 
							if len(cat_g_band) == 1: 
								mags.append(cat_g_band)
								magarrs.append(cat_g_err)
								checker.append(cat_g_band)
								
						if len(match) == 0: 
							
						 		
							
						
