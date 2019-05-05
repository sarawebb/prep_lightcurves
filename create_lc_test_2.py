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

field = 'Centaurus'
year = '2018'
month = '06'
sorted_image_path = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/image_mjd_lists/'+year+'_'+month+'_'+field+'_sorted.ascii'


source_lists = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/source_lists/'+ field +'/'+year+'/'+ month+'/'

start_path = '/mnt/dwf/archive_NOAO_data/data_outputs/'+ year+'/'+month+'/'+ field+'/g_band/single/'

image_names = np.loadtxt(sorted_image_path, unpack = True, skiprows = 1, usecols =[0], dtype = str )
image_mjds = np.loadtxt(sorted_image_path, unpack =True, skiprows= 1, usecols = [1])

sorted_table = Table()
sorted_table['filenames'] = image_names
sorted_table['mjds'] = image_mjds 
print(sorted_table)

field_source_RA = []
field_source_DEC = []
field_source_g_mag_AUTO= [] 
field_source_g_mag_err_AUTO = [] 
field_source_g_mag_APER = []
field_source_g_mag_err_APER = []
field_source_MJD = []

firstfile = image_names[0]
otherfiles = image_names[1:]
otherfiles_mjds = image_mjds[1:]
print(firstfile)


firstfound = False
#print(os.listdir(source_lists))

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
			#print(field_source_RA[1:10])
			
	for name, mjds in zip(otherfiles, otherfiles_mjds): 
		for filename in os.listdir(source_lists):
			if filename.startswith(name): 
				
				print(source_lists + '/' + filename) 
				
				#print(field_source_RA[1:10])
				try: 
					print('length of field_source_RA: ' + str(len(field_source_RA)))
					print('length of field_source_DEC: ' + str(len(field_source_DEC)))
			
				
					testing = np.array(field_source_RA)
					testing_dec = np.array(field_source_DEC)
				
					first_cat = np.empty((len(field_source_RA), 2), dtype = np.float64)
					first_cat[:, 0] = field_source_RA 
					first_cat[:, 1] = field_source_DEC 
				
				
				
					#print('FIRST CAT TEST: ' +  str(first_cat)) 
					#print('TEST: ' + str(first_cat[1,0]))
			
				
					#print(first_cat)
			
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
					#print('unpacked SE CAT')
				
					#print('NEW RA: ' + str(len(new_ra))) 
					new_cat = np.empty((int(len(new_ra)), 2), dtype = np.float64)
					new_cat[:,0] = new_ra
					new_cat[:,1] = new_dec
				
					#print('TEST NEW: ' + str(new_cat[1,0]))
				
					#new_cat[:,2] = new_gmagauto
					#new_cat[:,3] = new_gmagautoerr
					#new_cat[:,4] = new_gmagaper
					#new_cat[:,5] = new_gmagapererr
					#new_cat[:,6] = new_mjd 

					max_radius = .5/3600 #.5 arc second 
					#print('first_len: ' + str(len(first_cat)))
					#print('new_len: ' + str(len(new_cat)))
				
					if len(first_cat) >= len(new_cat): 
						print('GREATER THEN')
						dist_between, ind_row = crossmatch_angular(first_cat, new_cat, max_radius)
						match =~np.isinf(dist_between)
					
						#print('first_cat: ' + str(len(first_cat)))
						#print('new_cat: ' + str(len(new_cat)))
						#print('match: ' + str(len(match)))
					
					
						match_table = Table()
						match_table['matched_true_false'] = match 
						match_table['match_ID'] = ind_row
						match_table['matched_firstcat_RA'] = first_cat[:,0]
						match_table['matched_firstcat_DEC'] = first_cat[:,1] 
				
				
				
						new_cat_line_ID = range(len(match))
						#print('LENGTH OF MATCH: ' + str(len(match)))
						match_table['match_table_false_id'] = new_cat_line_ID 
					
					
						new_cat_match_true = []
						new_cat_match_false = []
						new_cat_ID_true = []
						new_cat_ID_false = []
						first_cat_RA_true = []
						first_cat_DEC_true = []
				
				
						first_cat_RA_false= []
						first_cat_DEC_false = []
					
						print('LENGTH TRUE: ' + str(len(first_cat_match_true)) )
						print('LENGTH FALSE: ' + str(len(first_cat_match_false))) 
						print('LENGTH ID TRUE: ' + str(len(first_cat_ID_true)))
						print('LENGTH ID FALSE: ' + str(len(first_cat_ID_false)))
						print('LENGTH OF NEWCAT: ' + str(len(new_cat)))
						print('LENGTH OF FIRSTCAT: ' + str(len(first_cat)))
				
						for row in match_table: 
							if row['matched_true_false'] == True:
						 		new_cat_match_true.append(row['matched_true_false'])
						 		new_cat_ID_true.append(row['match_ID'])
						 		first_cat_RA_true.append(row['matched_firstcat_RA'])
						 		first_cat_DEC_true.append(row['matched_firstcat_DEC'])
						 	
						#print('new_cat_match_true:' + str(len(new_cat_match_true))) 
					
						#print('CHECK: ' + str(first_cat_RA_true))
					
					
					
						for row in match_table: 
							if row['matched_true_false'] == False:
								new_cat_match_false.append(row['matched_true_false'])
								new_cat_ID_false.append(row['match_table_false_id'])
								first_cat_RA_false.append(row['matched_firstcat_RA'])
								first_cat_DEC_false.append(row['matched_firstcat_DEC'])
						
					
						#print('new_cat_match_false:' + str(len(new_cat_match_false))) 
					
						new_cat_RA_true  = []
						new_cat_DEC_true  = []
						new_cat_gmag_auto_true  = []
						new_cat_gmag_auto_err_true  = []
						new_cat_gmag_aper_true  = []
						new_cat_gmag_aper_err_true  = []
						new_cat_mjd_true = []
				
						for j in new_cat_ID_true: 
							new_cat_RA_true.append(new_cat[j,0])
							new_cat_DEC_true.append(new_cat[j,1])
					
					
					
						new_cat_RA_false = []
						new_cat_DEC_false  = []
						new_cat_gmag_auto_false  = []
						new_cat_gmag_auto_err_false  = []
						new_cat_gmag_aper_false  = []
						new_cat_gmag_aper_err_false  = []
						new_cat_mjd_false = []
				
						for k in new_cat_ID_false: 
							new_cat_RA_false.append(first_cat[k,0])
							new_cat_DEC_false.append(first_cat[k,1])
					
					
				
						match_true_table = Table () 
						match_true_table['RA'] = new_cat_RA_true
						match_true_table['DEC'] = new_cat_DEC_true
					
				
						match_false_table = Table()
						match_false_table['RA'] = new_cat_RA_false
						match_false_table['DEC'] = new_cat_DEC_false
					
						print('TRUE: ' + str(len(match_true_table)))
						print('FALSE: ' + str(len(match_false_table)))
					
					
					
					elif len(first_cat) <= len(new_cat):
						print('LESS THEN')
						dist_between, ind_row = crossmatch_angular( new_cat, first_cat, max_radius)
						match =~np.isinf(dist_between)
						#print(match)
				
						#print('first_cat: ' + str(len(first_cat)))
						#print('new_cat: ' + str(len(new_cat)))
						#print('match: ' + str(len(match)))
					
						#print('NEW_CAT_TEST: ' + str(new_cat)) 
					
						match_table = Table()
						match_table['matched_true_false'] = match 
						match_table['match_ID'] = ind_row
						match_table['matched_firstcat_RA'] = new_cat[:,0]
						match_table['matched_firstcat_DEC'] = new_cat[:,1] 
					
						new_cat_line_ID = range(len(match))
						#print('LENGTH OF MATCH: ' + str(len(match)))
						match_table['match_table_false_id'] = new_cat_line_ID 
						#print(match_table)
					
						first_cat_match_true = []
						first_cat_match_false = []
						first_cat_ID_true = []
						first_cat_ID_false = []
						new_cat_RA_true = []
						new_cat_DEC_true = []
						new_cat_RA_false = []
						new_cat_DEC_false = []
					
						#print(match_table['matched_true_false'])
						for row in match_table: 
							if row['matched_true_false'] == True:
						 		first_cat_match_true.append(row['matched_true_false'])
						 		first_cat_ID_true.append(row['match_ID'])
						 		new_cat_RA_true.append(row['matched_firstcat_RA'])
						 		new_cat_DEC_true.append(row['matched_firstcat_DEC'])
						#print('CHECK: ' + str(new_cat_RA_true) )	
						for row in match_table: 
							if row['matched_true_false'] == False:
								first_cat_match_false.append(row['matched_true_false'])
								first_cat_ID_false.append(row['match_table_false_id'])
								new_cat_RA_false.append(row['matched_firstcat_RA'])
								new_cat_DEC_false.append(row['matched_firstcat_DEC'])
							
						#print('RA OF FALSE: ' + str(new_cat_RA_false))
					
						print('LENGTH TRUE: ' + str(len(first_cat_match_true)) )
						print('LENGTH FALSE: ' + str(len(first_cat_match_false))) 
						print('LENGTH ID TRUE: ' + str(len(first_cat_ID_true)))
						print('LENGTH ID FALSE: ' + str(len(first_cat_ID_false)))
						print('LENGTH OF NEWCAT: ' + str(len(new_cat)))
						print('LENGTH OF FIRSTCAT: ' + str(len(first_cat)))
						#print(first_cat_ID_false)
						first_cat_RA_true  = []
						first_cat_DEC_true  = []
				
						for j in first_cat_ID_true: 
							#print('J: ' + str(j))
							RA_test = first_cat[j, 0]
							#print('RA_test: ' + str(RA_test))
							DEC_test = first_cat[j, 1]
							first_cat_RA_true.append(RA_test)
							first_cat_DEC_true.append(DEC_test)
					
						#print('CHECK: ' + str(first_cat_RA_true[1,0]))
					
						first_cat_RA_false = []
						first_cat_DEC_false  = []
					
						# The sources in the new cat that weren't found in the first catalog! 
					
						for k in first_cat_ID_false: 
							#print('K: ' + str(k))
							RA = new_cat[k, 0]
							DEC = new_cat[k, 1]
							first_cat_RA_false.append(RA)
							first_cat_DEC_false.append(DEC)
					
						#print('LENGTH FALSE: ' + str(first_cat_RA_false))
						#print('LENGTH TRUE: ' + str(first_cat_RA_true))
						match_true_table = Table () 
						match_true_table['RA'] = first_cat_RA_true
						match_true_table['DEC'] = first_cat_DEC_true
					
						match_false_table = Table()
						match_false_table['RA'] = first_cat_RA_false
						match_false_table['DEC'] = first_cat_DEC_false
					
						#print('TRUE: ' + str(len(match_true_table)))
						#print('FALSE: ' + str(len(match_false_table)))
					
					master_RA = []
					master_DEC = []
					master_gmag_auto = []
					master_mjd = []
				
				
				
					for u in match_true_table['RA']:
						master_RA.append(u) 
					
					for o in match_false_table['RA']: 
						master_RA.append(o) 
					
					for p in match_true_table['DEC']:
						master_DEC.append(p)
					
					for l in match_false_table['DEC']: 
						master_DEC.append(l)
				
				
					print('MASTER RA: ' + str(len(master_RA)))
					print('MASTER DEC: ' + str(len(master_DEC)))
				
				
					field_source_RA = master_RA 
					field_source_DEC = master_DEC 
				
					print(field_source_RA[1:10])
				
				
				
				except: 
					print('PASSED')
	
	final_table = Table()
	final_table['RA'] = field_source_RA
	final_table['DEC'] = field_source_DEC
	
	
	
	output = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/total_field_sources/fields/' + year + '_' + month + '_' +field + '_full_field_sources.ascii'
	final_table.write(output, format ='ascii', overwrite = True)
	
