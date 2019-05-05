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


field = ''
year = ''
month = ''

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
	SOR_X[1,1] = row['RA']
	SOR_X[1,2] = row['DEC']
	print(SOR_X)
	for row1 in sorted_table:
		image_path = str(start_path + row['filenames'])
		for i in im_paths: 
			if i == image_path[:len(image_path)-5]:
				enter_cats = os.path.join(i, 'final_source_cats')
				os.chdir(enter_cats)
				for filename in os.listdir('.'):
					if filename.endswith('cat_CORRECTED.ascii'):
						cwd.os.getcwd()
						cat_to_open = os.path.join(cwd, filename)
						RA, DEC, g_mag_auto, g_mag_auto_err, g_mag_aper, g_mag_aper_error = np.loadtxt( cat_to_open, unpack =True, skiprows = 1)
						
						CAT_X = np.empty((len(RA), 4), dtype = np.float64)
						CAT_X[:, 0] = RA 
						CAT_X[:, 1] = DEC
						CAT_X[:, 2] = g_mag_auto
						CAT_X[:, 3] = g_mag_auto_err 
						
						
