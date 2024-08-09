#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:18:23 2024

@author: lorrainemarcelin
"""
'''creating a table of the parameters to constantly be edited'''
import numpy as np
from glob import glob
from astropy.table import Table
from pathlib import Path
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Column

#creating list of redshift
hdulist = fits.open('jw_o3_candidates_VI_v20240115v2.fits')
center_coords= Table(hdulist[1].data) #putting data to astropy table
good_coords = center_coords[np.where(center_coords['Score_VI']==3)]

XX=(good_coords["ra"])# getting list of X and Y coordinates from table
YY =(good_coords["dec"])


Z03_2D = [''] * 124

is_overdensity = np.zeros(len(good_coords)) # initialize an array of zeros; if it’s in an overdensity then we change it to a nonzero value
protocluster = good_coords[np.where(good_coords['ZO3_2D']>=6.5)]
    
is_overdensity[np.where(good_coords['ZO3_2D']>=6.5)] = 1
overdensity1 =  good_coords[np.where((good_coords['ZO3_2D']>=5.35)&(good_coords['ZO3_2D']<=5.45))]
is_overdensity[np.where((good_coords['ZO3_2D']>=5.35)&(good_coords['ZO3_2D']<=5.45))] = 2
overdensity2 =  good_coords[np.where((good_coords['ZO3_2D']>=6.15)&(good_coords['ZO3_2D']<=6.3))]
is_overdensity[np.where((good_coords['ZO3_2D']>=6.15)&(good_coords['ZO3_2D']<=6.3))] = 3
field = good_coords[np.where(is_overdensity==0)] 


for i,j in enumerate(good_coords['ZO3_2D']):
    print(i)
    if j >=6.5:
        Z03_2D[i]= 'z >=6.5'
    else:
        Z03_2D[i]=' 5.3 < z < 6.5'
        


       



params_table = Table(names=['xc_centroid', 
                            'yc_centroid', 
                            'ellipticity_centroid',
                            'elongation_centroid',
                            'rhalf_circ','rhalf_ellip','Gini','M20','Concentration',
                            'Asymmetry','sersic_n','sersic_chi2_dof',
                            'flag','flag_sersic'], 
                     dtype=('float', 
                            'float', 
                            'float',
                            'float','float','float',
                            'float',
                            'float',
                            'float',
                            'float',
                            'float','float',
                            'float','float'))
line_numbers = (1,2,3,4,13,14,17,18,22,23,28,33,37,38) #indices of desired parameters

num = np.arange(124)

for i in num:
    
     file_path = glob('galaxy'+str(i)+'_*/galaxy*_parameters.txt')[0]
     file = open(file_path,'r') # reading in the text file
     read = file.readlines() #line variable
     for line in line_numbers: # getting numerical value of each parameter
        modified = read[line].split('=')
        read[line]= modified[1] 
        
     parameter_row = [read[1], read[2], read[3],read[4], read[13],read[14],read[17],read[18],read[22],
                     read[23],read[28],read[33],read[37],read[38]]
     params_table.add_row(parameter_row)
     print(file_path)
        
# adding reshift column
params_table.add_column(Z03_2D,name = 'ZO3_2D')
ascii.write(params_table, 'Parameters_Table.dat', overwrite=True)


