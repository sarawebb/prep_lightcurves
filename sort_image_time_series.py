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


path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/'

for i in os.listdir(path): 
	if i.endswith('.ascii'): 
		
		try: 
			filename = np.loadtxt(path + i, unpack = True, skiprows = 1, usecols =[0], dtype = str)
			mjd = np.loadtxt(path + i, unpack =True, skiprows= 1, usecols = [1])
		
			mydict = SortedDict()
			for j in range(len(filename)): 
				mydict[mjd[j]] = filename[j], mjd[j]
		
		 	
			#print(mydict)
			
		
			filenames_sorted = []
			mjd_sorted = []
			
			for key, values in mydict.items():
				filenames_sorted.append(values[0])
				mjd_sorted.append(values[1])
				
			#print(filenames_sorted) 
			#print(mjd_sorted)
				
			info = Table()
			info['filename'] = filenames_sorted 
			info['MJD'] = mjd_sorted
			print(info)
			output = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/'+ i[:len(i)-6] +'_sorted.ascii' 
			print(output)
			info.write(output, format='ascii', overwrite = True)
 
		except: 
			pass
