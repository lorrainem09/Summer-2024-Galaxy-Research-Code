#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 13:25:14 2024

@author: lorrainemarcelin
"""

import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.io import fits

from astropy.coordinates import SkyCoord


from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70 * u.km/u.s/u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3) 


hdulist = fits.open('jw_o3_candidates_VI_v20240115v2.fits')
center_coords= Table(hdulist[1].data) #putting data to astropy table
good_gals = center_coords[np.where(center_coords['Score_VI']==3)]
good_coords = SkyCoord(ra=good_gals['ra'], dec=good_gals['dec'], unit=u.degree)

XX=(good_gals["ra"])# getting list of X and Y coordinates from table
YY =(good_gals["dec"])


#Getting redshift of galaxies

SED_properties = Table.read('OIII_SED_properties_022624_delayedtau.dat', format='ascii')
redshift = SED_properties['z']
separation=[]
galaxy_pair=[]
restframe_velocity=[]

#Calculating separation distance between closest galaxy
galaxy=0 #keep track of first galaxy
flag_array = np.zeros(len(good_coords)) # keep track of galaxies that have a match pair
threshold = 30/3600 # 30 arcsec threshold
flag_value = 0

for gal, (RA1,DEC1) in enumerate(zip(XX,YY)):

    c1 = SkyCoord(RA1*u.deg, DEC1*u.deg, frame='icrs')
    
    if flag_array[gal] != 0:
        print("galaxy %i has already been tagged, moving on."%gal)
        continue
        
    angular_separation = c1.separation(good_coords).degree
    match = np.where((angular_separation>0)&(angular_separation<=threshold))[0] # find the closest match that is not itself and is within the distance threshold
 
    if len(match) == 0: # no matches, just continue
        continue
    elif len(match) > 0:
        # if there's more than 2 galaxies in the merger,'match' array might have length>1, so just take the minimum
        #print(gal,'match =',match)
        match_separation = []
        match_index=0
        for index in match: #creating array of the separation values of the matched galaxies
            match_separation.append(float(angular_separation[index]))
            match_index += 1
        match_separation = np.array(match_separation)
        match_gal = np.where(angular_separation==np.min(match_separation)) #closest galaxy number
        
        pair = match[np.where(match==match_gal[0])]#finding minimum matched angular separation
        
        z1 = redshift[gal]
        z2 = redshift[pair]
        zmed = np.median([z1, z2[0]])
        dA = cosmo.angular_diameter_distance(zmed).to('kpc').value * np.pi/180 # this has units of kpc per degree
        physical_sep = angular_separation[pair] * dA 
        #print(physical_sep)
        
        # also take the redshift difference, equation 2
        c = 3.0e5 # speed of light in km/s
        delta_v = c * (np.abs(z1-z2)) / (1 + zmed)#Calculating restframe velocity between galaxy pairs
        restframe_velocity.append(delta_v)
        #print(delta_v)
        
        if physical_sep <= 50 and delta_v <= 315 or physical_sep<=100 and delta_v<=100: # within 30 kpc and 5000 km/s, you can change these if you want
            flag_value +=1 # label the pair - at the end, flag_array will contain integers that indicate which galaxies are matched
            flag_array[pair] = flag_value
            flag_array[gal] = flag_value # flag_array will now tie together galaxy C1 & its match with a unique number
            galaxy_pair.append(gal)
            galaxy_pair.append(pair)
            print(gal, pair)
            


#Calculating merger pair fraction
merger_fraction=len(np.where(flag_array!=0)[0])/len(good_coords)/2 #divide by 2 so you’re not double counting



