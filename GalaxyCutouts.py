#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:58:10 2024

@author: lorrainemarcelin
"""
#import library
import numpy as np
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
#reading fits file with center coordinates

hdulist = fits.open('jw_o3_candidates_VI_v20240115v2.fits')
center_coords= Table(hdulist[1].data) #putting data to astropy table
good_coords = center_coords[np.where(center_coords['Score_VI']==3)]

XX=(good_coords["ra"])# getting list of X and Y coordinates from table
YY =(good_coords["dec"])

#Determining which overdensity each cutout is in
is_overdensity = np.zeros(len(good_coords)) # initialize an array of zeros; if it’s in an overdensity then we change it to a nonzero value
protocluster = good_coords[np.where(good_coords['ZO3_2D']>=6.5)]
is_overdensity[np.where(good_coords['ZO3_2D']>=6.5)] = 1
overdensity1 =  good_coords[np.where((good_coords['ZO3_2D']>=5.35)&(good_coords['ZO3_2D']<=5.45))]
is_overdensity[np.where((good_coords['ZO3_2D']>=5.35)&(good_coords['ZO3_2D']<=5.45))] = 2
overdensity2 =  good_coords[np.where((good_coords['ZO3_2D']>=6.15)&(good_coords['ZO3_2D']<=6.3))]
is_overdensity[np.where((good_coords['ZO3_2D']>=6.15)&(good_coords['ZO3_2D']<=6.3))] = 3
field = good_coords[np.where(is_overdensity==0)] 

'''
galaxy_file = fits.open ('J0305M3150_mosaic_F356W_psfmatch.fits')#importing fits file of image
#print(galaxy_file.info())
image= galaxy_file[0].data/galaxy_file[0].header['PHOTMJSR']*galaxy_file[0].header['XPOSURE']#getting image data and converting it to counts
w = wcs.WCS(galaxy_file[0].header)#getting wcs from header




'''for i, (x, y) in enumerate(zip(XX, YY)):#simultanously going through ra and dec coords

    center = SkyCoord(ra=x, dec=y, unit=u.deg) #cutout center coordinates
    size = u.Quantity((101,101),u.pixel) #cutout size
    cutout = Cutout2D(image, center, size, wcs=w)
    filtername = 'F356W_psfmatch'
    plt.imshow(cutout.data, origin = 'lower', ) #displaying image cutout
    plt.savefig('cutout'+str(i)+filtername+'.png') #saving plot as cutout name and filter
    overdensity = is_overdensity[i] # what is the value of our overdensity counter at this position?
    if overdensity == 1:
        filename = 'galaxy' + str(i) + '_' + filtername + '_protocluster.fits'
    elif overdensity == 2:
        filename = 'galaxy' + str(i) + '_' + filtername + '_overdensity1.fits' 
    elif overdensity == 3:
        filename = 'galaxy' + str(i) + '_' + filtername + '_overdensity2.fits'
    elif overdensity == 0 :
        filename = 'galaxy' + str(i) + '_' + filtername + '_field.fits
    

    newheader = cutout.wcs.to_header() #creating new header
    fnew = fits.PrimaryHDU(data=cutout.data, header=newheader) # creating new fits file
    fnew.writeto(filename)''' #writing cutout to fits file 


    
#save each galaxy as fits file 
#globe filter       

#zip for loop
#FW200
#Z03 redshift to one decimal point
center = SkyCoord(ra =XX[65],dec =YY[65], unit = u.deg) #center position coordinate
size = u.Quantity((101,101),u.pixel) #cutout size 
cutout = Cutout2D(image,center,size, wcs=w ) #Test cutout
plt.imshow(cutout.data,origin = 'lower')
plt.savefig('CutoutTest.png')
plt.show()
newheader = cutout.wcs.to_header() #creating new header
fnew = fits.PrimaryHDU(data=cutout.data, header=newheader) # creating new fits file
fnew.writeto('galaxy65_F356W_psfmatch_overdensity2.fits')


    
    