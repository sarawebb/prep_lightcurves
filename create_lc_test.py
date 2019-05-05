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

## TESTING FOR ONE FIELD IN 4hr field 2015 ###


sorted_image_path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/2015_01_4hr_sorted.ascii'
#cat_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/*/final_source_cats/' 
#cat_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/*/' 

#cat_paths = glob.glob(cat_path)

source_lists = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/source_lists'

start_path = '/mnt/dwf/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/'

image_names = np.loadtxt(sorted_image_path, unpack = True, skiprows = 1, usecols =[0], dtype = str )
image_mjds = np.loadtxt(sorted_image_path, unpack =True, skiprows= 1, usecols = [1])

sorted_table = Table()
sorted_table['filenames'] = image_names
sorted_table['mjds'] = image_mjds 


field_source_RA = []
field_source_DEC = []
field_source_g_mag_AUTO= [] 
field_source_g_mag_err_AUTO = [] 
field_source_g_mag_APER = []
field_source_g_mag_err_APER = []
field_source_MJD = []

firstfile = image_names[0]
otherfiles = image_names[1:]
#print(otherfiles)


firstfound = False


while firstfound == False: 
	for filename in os.listdir(source_lists):
		if filename.startswith(str(firstfile)):
			print('FIRST: ' + filename)
			ra, dec, g_mag_auto, g_mag_auto_err, g_mag_aper, g_mag_aper_err, mjd = np.loadtxt(source_lists + '/' +filename, unpack = True, skiprows=1)
			for k in ra:
				field_source_RA.append(k)
			for k in dec: 
				field_source_DEC.append(k)
			for k in g_mag_auto: 
				field_source_g_mag_AUTO.append(k)
			for k in g_mag_auto_err: 
				field_source_g_mag_err_AUTO.append(k)
			for k in g_mag_aper:
				field_source_g_mag_APER.append(k)
			for k in g_mag_aper_err: 
				field_source_g_mag_err_APER.append(k)
			for k in mjd:
				field_source_MJD.append(k)
			firstfound = True
	
	for name in otherfiles: 
		for filename in os.listdir(source_lists):
			if filename.startswith(name): 
				print(filename) 
				first_cat = np.empty((len(field_source_RA), 7), dtype = np.float64)
				first_cat[:,0] = field_source_RA
				first_cat[:,1] = field_source_DEC
				first_cat[:,2] = field_source_g_mag_AUTO
				first_cat[:,3] = field_source_g_mag_err_AUTO
				first_cat[:,4] = field_source_g_mag_APER
				first_cat[:,5] = field_source_g_mag_err_APER
				first_cat[:,6] = field_source_MJD
				
				first_cat_ra_dec = np.empty((len(field_source_RA), 2), dtype = np.float64) 
				first_cat_ra_dec[:, 0] = field_source_RA 
				first_cat_ra_dec[:, 1] = field_source_DEC
				
				ra, dec, g_mag_auto, g_mag_auto_err, g_mag_aper, g_mag_aper_err, mjd = np.loadtxt(source_lists + '/' +filename, unpack = True, skiprows=1)
			
				new_ra = []
				new_dec = []
				new_gmagauto = []
				new_gmagautoerr = []
				new_gmagaper = []
				new_gmagapererr = []
				new_mjd = []
				
				for k in ra:
					new_ra.append(k)
				for k in dec: 
					new_dec.append(k)
				for k in g_mag_auto: 
					new_gmagauto.append(k)
				for k in g_mag_auto_err: 
					new_gmagautoerr.append(k)
				for k in g_mag_aper:
					new_gmagaper.append(k)
				for k in g_mag_aper_err: 
					new_gmagapererr.append(k)
				for k in mjd:
					new_mjd.append(k)
				
				#print(len(new_ra)) 
				new_cat = np.zeros((int(len(new_ra)), 7), dtype = np.float64)
				new_cat[:,0] = new_ra
				new_cat[:,1] = new_dec
				new_cat[:,2] = new_gmagauto
				new_cat[:,3] = new_gmagautoerr
				new_cat[:,4] = new_gmagaper
				new_cat[:,5] = new_gmagapererr
				new_cat[:,6] = new_mjd 
				
				new_cat_ra_dec = np.zeros((int(len(new_ra)), 2), dtype = np.float64)
				new_cat_ra_dec[:,0] = new_ra
				new_cat_ra_dec[:,1] = new_dec
				
				print(new_cat[:,0])
				print(first_cat[:,0])
				max_radius = .5/3600 #1 arc second 
				print('first_len: ' + str(len(first_cat)))
				print('new_len: ' + str(len(new_cat)))
				
				dist_between, ind_row = crossmatch_angular(first_cat_ra_dec, new_cat_ra_dec, max_radius)
				match =~np.isinf(dist_between)
				for item in match: 
					if item == True: 
						print(item)
				#print(match)
				
				'''
				elif len(new_cat) <= len(first_cat):
					('LESS THEN')
					dist_between, ind_row = crossmatch_angular( new_cat, first_cat)
					match =~np.isinf(dist_between)
					print(match)
					
					match_table = Table()
					match_table['matched_true_false'] = match 
					match_table['match_ID'] = ind_row
					match_table['matched_firstcat_RA'] = first_cat[:,0]
					match_table['matched_firstcat_DEC'] = first_cat[:,1] 
					match_table['matched_firstcat_gmag_auto'] = first_cat[:,2] 
					match_table['matched_firstcat_gmag_auto_err'] = first_cat[:,3]
					match_table['matched_firstcat_gmag_aper'] = first_cat[:,4]
					match_table['matched_firstcat_gmag_aper_err'] = first_cat[:,5]
					match_table['matched_firstcat_MJD'] = first_cat[:,6]
					
					new_cat_match_true = []
					new_cat_match_false = []
					new_cat_ID_true = []
					new_cat_ID_false = []
					first_cat_RA_true = []
					first_cat_DEC_true = []
					first_cat_gmag_auto_true = []
					first_cat_gmag_auto_err_true = []
					first_cat_gmag_aper_true = []
					first_cat_gmag_aper_err_true = []
					first_cat_mjd_true = []
					
					first_cat_RA_false= []
					first_cat_DEC_false = []
					first_cat_gmag_auto_false = []
					first_cat_gmag_auto_err_false = []
					first_cat_gmag_aper_false = []
					first_cat_gmag_aper_err_false = []
					first_cat_mjd_false = []
					
					for row in match_table: 
						if row['matched_true_false'] == True:
							 new_cat_match_true.append(row['matched_true_false'])
							 new_cat_ID_true.append(row['match_ID'])
							 first_cat_RA_true.append(row['matched_firstcat_RA'])
							 first_cat_DEC_true.append(row['matched_firstcat_DEC'])
							 first_cat_gmag_auto_true.append(row['matched_firstcat_gmag_auto'])
							 first_cat_gmag_auto_err_true.append(row['matched_firstcat_gmag_auto_err'])
							 first_cat_gmag_aper_true.append(row['matched_firstcat_gmag_aper'])
							 first_cat_gmag_aper_err_true.append(row['matched_firstcat_gmag_aper_err'])
							 first_cat_mjd_true.append(row['matched_firstcat_MJD'])
							 
					for row in match_table: 
						if row['matched_true_false'] == False:
							 new_cat_match_false.append(row['matched_true_false'])
							 new_cat_ID_false.append(row['match_ID'])
							 first_cat_RA_false.append(row['matched_firstcat_RA'])
							 first_cat_DEC_false.append(row['matched_firstcat_DEC'])
							 first_cat_gmag_auto_false.append(row['matched_firstcat_gmag_auto'])
							 first_cat_gmag_auto_err_false.append(row['matched_firstcat_gmag_auto_err'])
							 first_cat_gmag_aper_false.append(row['matched_firstcat_gmag_aper'])
							 first_cat_gmag_aper_err_false.append(row['matched_firstcat_gmag_aper_err'])
							 first_cat_mjd_false.append(row['matched_firstcat_MJD'])
					
					new_cat_RA_true  = []
					new_cat_DEC_true  = []
					new_cat_gmag_auto_true  = []
					new_cat_gmag_auto_err_true  = []
					new_cat_gmag_aper_true  = []
					new_cat_gmag_aper_err_true  = []
					new_cat_mjd_true = []
					
					for j in new_cat_match_true: 
						new_cat_RA_true.append(new_cat[j,0])
						new_cat_DEC_true.append(new_cat[j,1])
						new_cat_gmag_auto_true.append(new_cat[j, 2])
						new_cat_gmag_auto_err_true.append(new_cat[j, 3])
						new_cat_gmag_aper_true.append(new_cat[j, 4])
						new_cat_gmag_aper_err_true.append(new_cat[j, 5])
						new_cat_mjd_true.append(new_cat[j, 5])
						
					new_cat_RA_false = []
					new_cat_DEC_false  = []
					new_cat_gmag_auto_false  = []
					new_cat_gmag_auto_err_false  = []
					new_cat_gmag_aper_false  = []
					new_cat_gmag_aper_err_false  = []
					new_cat_mjd_false = []
					
					for k in new_cat_match_false: 
						new_cat_RA_false.append(new_cat[k,0])
						new_cat_DEC_false.append(new_cat[k,1])
						new_cat_gmag_auto_false.append(new_cat[k, 2])
						new_cat_gmag_auto_err_false.append(new_cat[k, 3])
						new_cat_gmag_aper_false.append(new_cat[k, 4])
						new_cat_gmag_aper_err_false.append(new_cat[k, 5])
						new_cat_mjd_false.append(new_cat[k, 5])
						
					
					match_true_table = Table () 
					match_true_table['RA'] = new_cat_RA_true
					match_true_table['DEC'] = new_cat_DEC_true
					match_true_table['gmag_auto'] = new_cat_gmag_auto_true
					match_true_table['gmag_auto_err'] = new_cat_gmag_auto_err_true	
					match_true_table['gmag_aper'] = new_cat_gmag_aper_true
					match_true_table['gmag_aper_err'] = new_cat_gmag_aper_err_true
					match_true_table['mjd'] = new_cat_mjd_true
					
					match_false_table = Table()
					match_false_table['RA'] = new_cat_RA_false
					match_false_table['DEC'] = new_cat_DEC_false
					match_false_table['gmag_auto'] = new_cat_gmag_auto_false
					match_false_table['gmag_auto_err'] = new_cat_gmag_auto_err_false
					match_false_table['gmag_aper'] = new_cat_gmag_aper_false
					match_false_table['gmag_aper_err'] = new_cat_gmag_aper_err_false
					match_false_table['mjd'] = new_cat_mjd_false
					
					print('match_true_table: ' + str(match_true_table))
					print('match_false_table: ' + str(match_false_table))
					
					'''
