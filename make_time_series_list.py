import numpy as np
import matplotlib as mpl
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
from astropy.table import Table, Column, join 
from astropy import units as u
import datetime as dt 
import glob

year = '2018'
month = '06'
field = 'Dusty11'
path = '/mnt/dwf/archive_NOAO_data/data/'+ year +'/'+ month +'/'+ field +'/g_band/single/'


filenames = []
mjd_images = []

for filename in os.listdir(path):
	if filename.endswith('.fits'): 
		hdulist = fits.open(path + '/'+ filename)
		head = hdulist[0].header
		mjd = head['MJD-OBS'] 
		filenames.append(filename)
		mjd_images.append(mjd)
		

print(filenames, mjd_images)
		
info = Table()
info['filename'] = filenames 
info['MJD'] = mjd_images 

output = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/' + year + '_' + month + '_' + field + '.ascii' 
info.write(output, format='ascii', overwrite = True)

