#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 11:59:14 2024

@author: lorrainemarcelin
"""

"""individual edits of the statmorph fitting"""

import numpy as np
import matplotlib.pyplot as plt
import statmorph
from astropy.io import fits
from statmorph.utils.image_diagnostics import make_figure
from glob import glob
import shutil

num = "20"
cutout_segmap_f356w = glob('galaxy'+num+'_F356W*Segmap*fits')[0]
cutout_image_f356w = glob('F356W_psfmatch/galaxy'+num+'_F356W_psfmatch_*fits')[0]
seg_image = fits.open(cutout_image_f356w)[0].data
segmap_cutout = fits.open(cutout_segmap_f356w)[0].data
 
psf = fits.open('J0305M3150_F356W_psf.fits')[0].data
gain = 1
source_morphs = statmorph.source_morphology(seg_image[30:70,30:70], segmap_cutout[30:70,30:70], gain=gain, psf=psf)
 # ^^^ list of objects corresponding to each labeled source in the image ^^^
morph = source_morphs[0] #First and only labeled source
fig = make_figure(morph)
fig.savefig('Statmorph_'+num+'.png', dpi=150)
plt.close(fig)
print('flag =', morph.flag)
print('flag_sersic =', morph.flag_sersic)
print ('morph.sersic_chi2_dof =', morph.sersic_chi2_dof)
#creating txt file to save paramaters
new_destination = glob('galaxy'+num+'_*/')[0]
file_path = "galaxy"+num+"_parameters.txt"
'''with open(file_path, 'w') as f:
    # Write content to the file
    print('BASIC MEASUREMENTS (NON-PARAMETRIC)',file=f)
    print('xc_centroid =',morph.xc_centroid,file=f)
    print('yc_centroid =', morph.yc_centroid,file=f)
    print('ellipticity_centroid =', morph.ellipticity_centroid,file=f)
    print('elongation_centroid =', morph.elongation_centroid,file=f)
    print('orientation_centroid =', morph.orientation_centroid,file=f)
    print('xc_asymmetry =', morph.xc_asymmetry,file=f)
    print('yc_asymmetry =', morph.yc_asymmetry,file=f)
    print('ellipticity_asymmetry =', morph.ellipticity_asymmetry,file=f)
    print('elongation_asymmetry =', morph.elongation_asymmetry,file=f)
    print('orientation_asymmetry =', morph.orientation_asymmetry,file=f)
    print('rpetro_circ =', morph.rpetro_circ,file=f)
    print('rpetro_ellip =', morph.rpetro_ellip,file=f)
    print('rhalf_circ =', morph.rhalf_circ,file=f)
    print('rhalf_ellip =', morph.rhalf_ellip,file=f)
    print('r20 =', morph.r20,file=f)
    print('r80 =', morph.r80,file=f)
    print('Gini =', morph.gini,file=f)
    print('M20 =', morph.m20,file=f)
    print('F(G, M20) =', morph.gini_m20_bulge,file=f)
    print('S(G, M20) =', morph.gini_m20_merger,file=f)
    print('sn_per_pixel =', morph.sn_per_pixel,file=f)
    print('C =', morph.concentration,file=f)
    print('A =', morph.asymmetry,file=f)
    print('S =', morph.smoothness,file=f)
    print()
    print('SERSIC MODEL',file=f)
    print('sersic_amplitude =', morph.sersic_amplitude,file=f)
    print('sersic_rhalf =', morph.sersic_rhalf,file=f)
    print('sersic_n =', morph.sersic_n,file=f)
    print('sersic_xc =', morph.sersic_xc,file=f)
    print('sersic_yc =', morph.sersic_yc,file=f)
    print('sersic_ellip =', morph.sersic_ellip,file=f)
    print('sersic_theta =', morph.sersic_theta,file=f)
    print('sersic_chi2_dof =', morph.sersic_chi2_dof,file=f)
    print()
    print('OTHER')
    print('sky_mean =', morph.sky_mean,file=f)
    print('sky_median =', morph.sky_median,file=f)
    print('sky_sigma =', morph.sky_sigma,file=f)
    print('flag =', morph.flag,file=f)
    print('flag_sersic =', morph.flag_sersic,file=f)

 
     
print(f"File '{file_path}' created successfully.")
#moving segmapcutout , images to galaxy folder


shutil.move('Statmorph_'+num+'.png',new_destination)
shutil.move( cutout_segmap_f356w,new_destination)
shutil.move(file_path,new_destination)
'''










