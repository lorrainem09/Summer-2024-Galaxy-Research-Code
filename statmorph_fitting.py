#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 12:42:40 2024

@author: lorrainemarcelin

"""
import numpy as np
import matplotlib.pyplot as plt
import statmorph
from astropy.io import fits
from statmorph.utils.image_diagnostics import make_figure

from glob import glob


index = np.arange(1,124)
badindex = []


for i in index:
        print(i)
        try:
            cutout_segmap_f356w = glob('galaxy'+str(i)+'_F356W*Segmap*fits')[0]
        except IndexError:
            badindex.append(i)
            continue
        try:
            cutout_image_f356w = glob('F356W_psfmatch/galaxy'+str(i)+'_F356W*fits')[0]
        except IndexError:
            continue
        print(cutout_image_f356w)
        print(cutout_segmap_f356w)
        
    
        seg_image = fits.open(cutout_image_f356w)[0].data
        segmap_cutout = fits.open(cutout_segmap_f356w)[0].data
        
        psf = fits.open('J0305M3150_F356W_psf.fits')[0].data
        gain = 1
        source_morphs = statmorph.source_morphology(seg_image[35:65,35:65], segmap_cutout[35:65,35:65], gain=gain, psf=psf)
        # ^^^ list of objects corresponding to each labeled source in the image ^^^
        morph = source_morphs[0] #First and only labeled source
        fig = make_figure(morph)
        fig.savefig('Statmorph_'+str(i)+'.png', dpi=150)
        plt.close(fig)
        
        
 
