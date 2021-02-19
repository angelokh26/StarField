

"""
Created on Mon Sep 16 22:35:16 2019

This program is used to locate and identify stars in a given field. 
It is assumed that the given files are calibrated properly.

@author: Angelo Hollett - angelokh26@gmail.com
"""

import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u

#================data access===================
tbl = ascii.read("field1textclean.txt")              # Read Field 1 data from text
tbl2 = ascii.read("field2text.txt")                  # Read Field 2 data from text
#==============================================

# field array and centre conversion to galactic Coordinates
# ------------------------------------------------------------------------------#
#Field 1
field1_c_icrs = SkyCoord(ra=184*u.degree, dec=32*u.degree, frame='icrs')
field1_gal = field1_c_icrs.galactic

#Field 2
field2_c_icrs = SkyCoord(ra=318*u.degree, dec=0.5*u.degree, frame='icrs')
field2_gal = field2_c_icrs.galactic

#======Compute the list of coordinates for each table and convert to Galactic Coordinates=========

field1_c = SkyCoord(ra=tbl["ra"]*u.degree, dec= tbl["dec"]*u.degree, frame='icrs')        #Field 1, Table 1
field2_c = SkyCoord(ra=tbl2["ra"]*u.degree, dec= tbl2["dec"]*u.degree, frame='icrs')      #Field 2, Table 2

field1_galactic = field1_c.galactic
field2_galactic = field2_c.galactic

#=============================Print out this info===================================
print ("============================================================================")
print ("The full table array of coordinates in Galactic coordinates is: ")
print ("============================================================================")
print field1_galactic


print ("============================================================================")
print ("The centre of field 1 in Galactic coordinates is: ") 
print field1_gal
print ("============================================================================")
print ("The centre of field 2 in Galactic coordinates is: ") 
print field2_gal
print ("============================================================================")

#Creating arrays of RA and DEC values (not used in final figures)
ra = coord.Angle(tbl["ra"]*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(tbl["dec"]*u.degree)

#============Create arrays of field psf magnitudes in R band=========================
field1_mag = np.array(tbl["psfMag_r"])
field2_mag = np.array(tbl2["psfMag_r"])

#===========calculate number of stars in each field brighter than r=21.5=============

# Field 1
count1 = 0
k = 21.5
for i in field1_mag : 
    if i < k : 
        count1 = count1 + 1

#field 2
count2 = 0
for i in field2_mag : 
    if i < k : 
        count2 = count2 + 1

#====================================================================================

#================Print out information regarding each field==========================
#-------------------------total stars in each field----------------------------------
#Field 1
field1_length = len(field1_galactic_l_array)
print ("=============================================================================")
print ("The number of stars in field 1 is:")
print field1_length

#Field 2
field2_length = len(field2_galactic_l_array)
print ("=============================================================================")
print ("The number of stars in field 2 is:")
print field2_length        

#--------------------Stars brighter than mag r = 21.5--------------------------------
print ("============================================================================")
print ("The number of stars brighter than r = 21.5 in field 1 is: ")
print (str(count))

print ("============================================================================")
print ("The number of stars brighter than r = 21.5 in field 2 is: ")
print (str(count2))

#========================Create the figure and subplots=============================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[10,20])
fig.subplots_adjust(hspace=0.3)

#=========================Field 1 centering and scaling==============================
# coordinate l
field1_galactic_l_array = np.array(field1_galactic.l)
field1_galactic_l_centered = field1_galactic_l_array - 178.499
field1_galactic_l_center_corr = np.multiply(field1_galactic_l_centered,0.1585)         #correct by factor cosb

#coordinate b
field1_galactic_b_array = np.array(field1_galactic.b)
field1_galactic_b_centered = field1_galactic_b_array - 80.888

#====================Field 2 centering and scaling====================================
# coordinate l
field2_galactic_l_array = np.array(field2_galactic.l)
field2_galactic_l_centered = field2_galactic_l_array - 51.218
field2_galactic_l_center_corr = np.multiply(field2_galactic_l_centered,0.8612)         #correct for factor cosb

#coordinate b
field2_galactic_b_array = np.array(field2_galactic.b)
field2_galactic_b_centered = field2_galactic_b_array + 30.548
#======================================================================================

#=========================Plot device -- Field 1======================================
#ax1.scatter(field1_galactic.l, field1_galactic.b, color = 'darkcyan', s=3)                       #Plot uncentered
ax1.scatter(field1_galactic_l_center_corr, field1_galactic_b_centered,
            color = 'k', s=10)                                                                    #plot centered
#=====================================================================================

#=========================Plot device -- Field 2======================================
#ax2.scatter(field2_galactic.l, field2_galactic.b, color = 'darkcyan', s=3)                       #Plot uncentered
ax2.scatter(field2_galactic_l_center_corr, field2_galactic_b_centered,
            color = 'k', s=3)   #Plot centered
# ====================================================================================


#======================Plot labels etc==========================================
ax1.set_title("Field 1 - Center: (l,b) = (178.499, 80.888)", fontsize=20)
ax1.set_xlabel("Relative Galactic Longitude (L) (Degrees)", fontsize=15)
ax1.set_ylabel("Relative Galactic Latitude (b) (Degrees)",fontsize=15)

ax2.set_title("Field 2 - Center: (l,b) = (51.218, -30.548)", fontsize=20)
ax2.set_xlabel("Relative Galactic Longitude (L) (Degrees)", fontsize=15)
ax2.set_ylabel("Relative Galactic Latitude (b) (Degrees)",fontsize=15)
# =============================================================================