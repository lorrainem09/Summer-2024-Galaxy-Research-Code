#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 12:37:07 2024

@author: lorrainemarcelin
"""

import shutil
from glob import glob
import os

cutout_images_f356w = glob('galaxy*F356W*fits')

for c in cutout_images_f356w:
	idx = c.split('galaxy')[1].split('_')[0]
	overdensity = c.split('_')[-1].split('.fits')[0]
	dirname = 'galaxy' + idx + '_' + overdensity
	os.makedirs(dirname)
	shutil.copy2(c, dirname)

