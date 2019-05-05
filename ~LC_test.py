
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


field = 'Centaurus'
year = '2018'
month = '06'

source_list = '/mnt/dwf/archive_NOAO_data/scripts/create_lc/total_field_sources/fields/' + year + '_' + month + '_' field +'_full_field_sources.ascii''

image lists
