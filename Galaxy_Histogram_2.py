#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:26:58 2024

@author: lorrainemarcelin
"""


#Creating historgrams of the half radius of the galaxies
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from astropy.io import fits
import seaborn as sns
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

data = ascii.read("Parameters_Table.dat")  
data_set = pd.read_csv('Parameters_Table.csv',delimiter=' ')
hdulist = fits.open('jw_o3_candidates_VI_v20240115v2.fits')
center_coords= Table(hdulist[1].data) #putting data to astropy table
good_coords = center_coords[np.where(center_coords['Score_VI']==3)]

XX=(good_coords["ra"])# getting list of X and Y coordinates from table
YY =(good_coords["dec"])
all_galaxy_shift = good_coords['ZO3_2D']



#red shift and clusters of galaxies
is_overdensity = np.zeros(len(good_coords)) # initialize an array of zeros; if it’s in an overdensity then we change it to a nonzero value
protocluster = good_coords[np.where(good_coords['ZO3_2D']>=6.5)]
protocluster_coords = protocluster
protocluster = protocluster['ZO3_2D']
is_overdensity[np.where(good_coords['ZO3_2D']>=6.5)] = 1
overdensity1 =  good_coords[np.where((good_coords['ZO3_2D']>=5.35)&(good_coords['ZO3_2D']<=5.45))]
overdensity1 = overdensity1['ZO3_2D']
is_overdensity[np.where((good_coords['ZO3_2D']>=5.35)&(good_coords['ZO3_2D']<=5.45))] = 2
overdensity2 =  good_coords[np.where((good_coords['ZO3_2D']>=6.15)&(good_coords['ZO3_2D']<=6.3))]
overdensity2 = overdensity2['ZO3_2D']
is_overdensity[np.where((good_coords['ZO3_2D']>=6.15)&(good_coords['ZO3_2D']<=6.3))] = 3
field = good_coords[np.where(is_overdensity==0)] 
field = field['ZO3_2D']

#indices of each overdensity (in order of data taken)
protocluster_index = []
overdensity1_index= [] 
overdensity2_index = []
field_index= []
for i,j in enumerate(all_galaxy_shift):
    if 6.5 <= j:
        protocluster_index.append(i)
    elif 5.35 <= j <= 5.45 :
        overdensity1_index.append(i)
    elif  6.15 <= j <= 6.3 :
        overdensity2_index.append(i)
    else: 
        field_index.append(i)


###############################################################################################

#Calculating Merger Candidates
cosmo = FlatLambdaCDM(H0=70 * u.km/u.s/u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3) 


hdulist = fits.open('jw_o3_candidates_VI_v20240115v2.fits')
center_coords= Table(hdulist[1].data) #putting data to astropy table
good_gals = center_coords[np.where(center_coords['Score_VI']==3)]
good_coords = SkyCoord(ra=good_gals['ra'], dec=good_gals['dec'], unit=u.degree)

XX=(good_gals["ra"])# getting list of X and Y coordinates from table
YY =(good_gals["dec"])


#Getting redshift of galaxies

SED_properties = Table.read('OIII_SED_properties_022624_delayedtau.dat', format='ascii')
redshift = all_galaxy_shift
galaxy_pair=[]
galaxy_pair_redshift=[]
galaxy_pair_zmed=[]
restframe_velocity=[]
nearest_gal_distance=[]

#Calculating separation distance between closest galaxy
galaxy=0 #keep track of first galaxy
flag_array = np.zeros(len(good_coords)) # keep track of galaxies that have a match pair
threshold = 30/3600 # 30 arcsec threshold
flag_value = 0

#getting the XX and YY coords of each sub group to iterate within each one


for gal, (RA1,DEC1) in enumerate(zip(XX,YY)):

    c1 = SkyCoord(RA1*u.deg, DEC1*u.deg, frame='icrs')
   

    
    if flag_array[gal] != 0:
        print("galaxy %i has already been tagged, moving on."%gal)
        continue
     
    angular_separation = c1.separation(good_coords).degree
    
    match = np.where((angular_separation>0)&(angular_separation<=threshold))[0] #find the closest match that is not itself and is within the distance threshold
    
    if len(match) == 0: # no matches, just continue
        continue
    elif len(match) > 0:
        print(match)
        # if there's more than 2 galaxies in the merger,'match' array might have length>1, so just take the minimum
        #print(gal,'match =',match)
        
                  
        #Of those galaxy matched redshifts determine the one with lower degree of separation if len(array) > 1
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
            galaxy_pair.append(int(gal))
            galaxy_pair.append(int(pair))
            galaxy_pair_redshift.append(float(z1))
            galaxy_pair_redshift.append(float(z2[0]))
            galaxy_pair_zmed.append(zmed)
            
###############################################################################################
         
     
for gal, (RA1,DEC1) in enumerate(zip(XX,YY)): #same as previous loop but without thresh hold minimum
    #finding distance to nearest neighbors for ALL galaxies

     c1 = SkyCoord(RA1*u.deg, DEC1*u.deg, frame='icrs')
     
     
      
     angular_separation = c1.separation(good_coords).degree
     match = np.where((angular_separation>0))[0] #find the closest match that is not itself and is within the distance threshold
     
     if len(match) == 0: # no matches, just continue
         continue
     elif len(match) > 0:
         # if there's more than 2 galaxies in the merger,'match' array might have length>1, so just take the minimum
         #print(gal,'match =',match)
     
         
         #Matched galaxies in the same subgroup are assigned to a separate list (added)
         z1 = redshift[gal] #redshift of gal 
         match_redshift = [] #array of match redshifts
         gal_redshift=[]
         if z1 in protocluster: #determining subgroup of gal redshift
             gal_redshift = 'protocluster'
         if z1 in field:
             gal_redshift='field'
         if z1 in overdensity1:
             gal_redshift='overdensity1'
         else:
             gal_redshift='overdensity2'
         #for loop iterating through matches to determine match redshift
         for i in match:
             if redshift[i] in protocluster:
                 match_redshift.append('protocluster')
             if redshift[i] in field:
                 match_redshift.append('field')
             if redshift[i] in overdensity1:
                 match_redshift.append('overdensity1')
             else:
                 match_redshift.append('overdensity2')
         z_match = []   
         for i in match: #comparing match redshifts to gal redshift, removing them as match if not in the same protocluster
             
             if i == gal_redshift:
                 z_match.append(match[i])
     
                   
         #Of those galaxy matched redshifts determine the one with lower degree of separation if len(array) > 1
         match_separation = []
         match_index=0
         for index in match: #creating array of the separation values of the matched galaxies
             match_separation.append(float(angular_separation[index]))
             match_index += 1
         match_separation = np.array(match_separation)
         match_gal = np.where(angular_separation==np.min(match_separation)) #closest galaxy number
         
         pair = match[np.where(match==match_gal[0])]#finding minimum matched angular separation
         
         nearest_gal_distance.append(np.min(match_separation))
         #galaxy_pair_redshift.append(gal_redshift)
         #full_galaxy_pair.append([gal,pair])
         #z2 = redshift[pair]
         
      
         

#Calculating merger pair fraction

merger_fraction=(len(galaxy_pair))/len(all_galaxy_shift)


protocluster_mergers = 0
overdensity1_mergers= 0 
overdensity2_mergers = 0
field_mergers = 0
for i in galaxy_pair_redshift:
    if 6.5 <= i :
        protocluster_mergers += 1
    elif 5.35 <= i <= 5.45 :
        overdensity1_mergers += 1
    elif  6.15 <= i <= 6.3 :
        overdensity2_mergers += 1
    else: 
        field_mergers += 1
protocluster_merger_fraction = protocluster_mergers/(len(galaxy_pair))
''' 
overdensity1_merger_fraction = overdensity1_mergers/(len(galaxy_pair))
overdensity2_merger_fraction = overdensity2_mergers/(len(galaxy_pair)))
field_merger_fraction = field_mergers/(len(galaxy_pair))) 
'''
#calculating Poisson error on percentage
merger_frac_error= np.sqrt(len(galaxy_pair))/len(all_galaxy_shift)


##############################################################################
half_radius = data_set['rhalf_circ']
rhalf=[]
for radius in half_radius:
    arc_radius= radius * 0.031 #converting from pixels to arcseconds
    rhalf.append(arc_radius)
    
protocluster_rhalf =[]
field_rhalf=[]
overdensity1_rhalf= []
overdensity2_rhalf=[]
galaxy_pair_rhalf=[]
for i in protocluster_index:
    protocluster_rhalf.append(rhalf[i])
for i in field_index:
    field_rhalf.append(rhalf[i])
for i in overdensity1_index:
    overdensity1_rhalf.append(rhalf[i])
for i in overdensity2_index:
    overdensity2_rhalf.append(rhalf[i])
    
for i in galaxy_pair:
    if isinstance(i, np.ndarray):  # Check if 'i' is a numpy array
       i = i[0] #take first element of numpy array
       
    galaxy_pair_rhalf.append(rhalf[i])
     
#plotting half radius of all galaxies and overdensities
'''


fig = sns.histplot(data=rhalf,color ='purple')
fig.set(xlabel='Half Radius', ylabel='Number of Galaxies')
plt.suptitle('Galaxy Half Radius')
plt.xlim(0,15)
plt.tick_params(axis='both',direction='in')
#plt.savefig('rhalf_hist.png')
plt.show()
'''


#Plots for rhalf



'''
    
plt.hist(protocluster_rhalf,edgecolor='black',color = 'plum')
plt.xlabel('Half Radius')
plt.suptitle('Protocluster Half Radius')
plt.ylabel('Number of Galaxies')
plt.tick_params(axis='both',direction='in')
#plt.savefig('protocluster_rhalf_hist.png')
plt.show()


plt.hist(field_rhalf,edgecolor='black',color = 'violet')
plt.xlabel('Half Radius')
plt.ylabel('Number of Galaxies')
plt.suptitle('Field Half Radius')
plt.tick_params(axis='both',direction='in')
#plt.savefig('field_rhalf_hist.png')
plt.show()


plt.hist(overdensity1_rhalf,edgecolor='black',color = 'purple')
plt.xlabel('Half Radius')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 1 Half Radius')
plt.tick_params(axis='both',direction='in')
#plt.savefig('Overdensity1_rhalf_hist.png')
plt.show()


plt.hist(overdensity2_rhalf,edgecolor='black',color = 'indigo')
plt.xlabel('Half Radius')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 2 Half Radius')
plt.tick_params(axis='both',direction='in')
#plt.savefig('Overdensity2_rhalf_hist.png')
plt.show()

plt.hist(rhalf, histtype = 'step', stacked='true',color='pink' )
plt.hist(protocluster_rhalf,label='protocluster',color ='plum')
plt.hist(field_rhalf,label = 'field',color='violet')
plt.hist(overdensity1_rhalf,label= 'overdensity 1',color='purple')
plt.hist(overdensity2_rhalf,label='overdensity 2',color='indigo')
plt.xlabel('Half Radius (arcsec)')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity Half Radius')
plt.tick_params(axis='both',direction='in')
plt.legend(loc='upper right')
plt.savefig('Layered_rhalf_hist.png')
'''
###############################################################################

#plotting distance to nearest neighbor vs half radius
distance_in_arcsec = []

#converting the angular distance from degrees to arcsec
for degree in nearest_gal_distance:
    degree_arcsec= degree * 3600 
    distance_in_arcsec.append(degree_arcsec)
    
#Getting distance to nearest neighbor in arcsec of each overdensity
protocluster_distance =[]
field_distance=[]
overdensity1_distance= []
overdensity2_distance=[]
galaxy_pair_distance=[]
for i in protocluster_index:
    protocluster_distance.append(distance_in_arcsec[i])
for i in field_index:
    field_distance.append(distance_in_arcsec[i])
for i in overdensity1_index:
    overdensity1_distance.append(distance_in_arcsec[i])
for i in overdensity2_index:
    overdensity2_distance.append(distance_in_arcsec[i])
for i in galaxy_pair:
    if isinstance(i, np.ndarray):  # Check if 'i' is a numpy array
       i = i[0] 
    galaxy_pair_distance.append(distance_in_arcsec[i])
'''
proper_distance_scale = cosmo.angular_diameter_distance(6).to('kpc') * (np.pi / 180)*(1/3600) # central_z would be the median redshift of your sample - z=6 is probably fine
fig = plt.figure(figsize=(4, 3))

ax = plt.gca()

# create a blank top axis in distance
distance_axis = np.linspace(np.min(protocluster_distance), 
                            np.max(protocluster_distance), 101) * proper_distance_scale 

ax2 = ax.twiny()


top_array = np.zeros(len(distance_axis))
ax2.plot(distance_axis, top_array, alpha=0)
ax2.set_xlabel("proper distance (kpc)", fontsize=11)

# create a blank right axis in half-light radius
ax3 = ax.twinx()
radius_axis = np.linspace(np.min(rhalf), np.max(rhalf), 101) * proper_distance_scale
right_array = np.zeros(len(radius_axis))
ax3.plot(right_array,radius_axis, alpha=0)
ax3.set_ylabel("half light radius (kpc)", fontsize=11)


ax.scatter(protocluster_distance,protocluster_rhalf,color='indigo',edgecolor='k',
           marker = '8', label='Protocluster')
ax.scatter(field_distance,field_rhalf,color='violet',edgecolor='k',
           marker='s', label='Field')
ax.scatter(overdensity1_distance,overdensity1_rhalf,color='mediumvioletred',
            edgecolor='k', marker='^', label='Overdensity 1')
ax.scatter(overdensity2_distance,overdensity2_rhalf,color='mediumpurple',
            edgecolor='k', marker='P', label='Overdensity 2')
ax.scatter(galaxy_pair_distance,galaxy_pair_rhalf,color='gold',
           marker = '*',label='Mergers')
ax.set_xlabel('angular distance (arcsec)',fontsize=11)
ax.set_ylabel('half radius (arcsec)',fontsize = 11)
ax.tick_params(axis='both',direction='in')
ax2.tick_params(axis='both',direction='in')
plt.tick_params(axis='both',direction='in')
ax2.axvline(100,linestyle='--',color='gray',alpha=0.5)
ax.legend(loc='upper right',
          fontsize=8,ncol=(2),labelspacing=0)

#plt.savefig('Distance_to_Nearest_Neighbor.png',bbox_inches='tight',pad_inches=0.1,dpi=250)
plt.show()
'''



###############################################################################

#Plotting Redshift of different overdensities

'''
#bin range
bins = np.arange(5.2,7.1,0.1)
xticks= np.arange(5,7,0.5)
plt.hist(all_galaxy_shift,histtype='step',stacked=True,bins=bins,color='pink')
plt.hist(protocluster, label='protocluster',bins=bins,color='plum')
#plt.hist(field, label='field')
plt.hist(overdensity1, label='overdensity1',bins=bins,color='purple')
plt.hist(overdensity2, label='overdensity2',bins=bins,color='indigo')

plt.xticks((5.5,6.0,6.5))
plt.tick_params(axis='both',direction='in')
plt.legend(loc='upper left')
plt.xlabel('Redshift')
plt.ylabel('Number of Galaxies')
plt.suptitle('Galaxy Overdensities')
plt.tick_params(axis='both',direction='in')
#plt.savefig('Galaxy_Overdensity_hist.png')
plt.show()
'''
###############################################################################

#Plotting Sersic n

'''

sersic_n = data['sersic_n']
bins = np.arange(0,5,0.1)
plt.hist(sersic_n,edgecolor='black',bins=bins,color='blueviolet')
plt.xlim(0,5)
plt.xlabel('Sersic n')
plt.ylabel('Number of Galaxies')
plt.tick_params(axis='both',direction='in')
plt.suptitle('Galaxy Sersic n')
plt.savefig('Galaxy_sersicn_hist.png')
plt.show()

protocluster_sersic =[]
field_sersic=[]
overdensity1_sersic= []
overdensity2_sersic=[]
for i in protocluster_index:
    protocluster_sersic.append(sersic_n[i])
for i in field_index:
    field_sersic.append(sersic_n[i])
for i in overdensity1_index:
    overdensity1_sersic.append(sersic_n[i])
for i in overdensity2_index:
    overdensity2_sersic.append(sersic_n[i])

plt.hist(protocluster_sersic,edgecolor='black',bins=bins,color = 'plum')
plt.xlabel('Sersic n')
plt.xlim(0,5)
plt.suptitle('Protocluster Sersic n')
plt.tick_params(axis='both',direction='in')
plt.ylabel('Number of Galaxies')
plt.savefig('Protocluster_sersicn_hist.png')
plt.show()


plt.hist(field_sersic,edgecolor='black',color = 'violet',bins=bins)
plt.xlabel('Sersic n')
plt.ylabel('Number of Galaxies')
plt.suptitle('Field Sersic n')
plt.tick_params(axis='both',direction='in')
plt.savefig('Field_sersicn_hist.png')
plt.show()


plt.hist(overdensity1_sersic,edgecolor='black',color = 'purple',bins=bins)
plt.xlabel('Sersic n')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 1 Sersic n')
plt.tick_params(axis='both',direction='in')
plt.savefig('Overdensity1_sersicn_hist.png')
plt.show()


plt.hist(overdensity2_sersic,edgecolor='black',color = 'indigo',bins=bins)
plt.xlabel('Sersic n')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 2 Sersic n')
plt.tick_params(axis='both',direction='in')
plt.savefig('Overdensity2_sersicn_hist.png')
plt.show()

plt.hist(sersic_n, histtype = 'step', stacked='true',bins=bins,color='pink')
plt.hist(protocluster_sersic,label='protocluster',bins=bins,color='plum')
plt.hist(field_sersic,label = 'field',bins=bins,color='violet')
plt.hist(overdensity1_sersic,label= 'overdensity 1',bins=bins,color ='purple')
plt.hist(overdensity2_sersic,label='overdensity 2',bins=bins,color='indigo')
plt.xlim(0,5)
plt.xlabel('Sersic n')
plt.ylabel('Number of Galaxies')
plt.suptitle('Sersic Index')
plt.legend(loc='upper right')
plt.tick_params(axis='both',direction='in')
plt.savefig('Layered_sersicn_hist.png')
plt.show()
'''

###############################################################################

#Plotting Asymmetry

asymmetry = data['Asymmetry']
protocluster_asymm =[]
field_asymm=[]
overdensity1_asymm= []
overdensity2_asymm=[]

for i in protocluster_index:
    protocluster_asymm.append(asymmetry[i])
for i in field_index:
    field_asymm.append(asymmetry[i])
for i in overdensity1_index:
    overdensity1_asymm.append(asymmetry[i])
for i in overdensity2_index:
    overdensity2_asymm.append(asymmetry[i])

bins = np.arange(-0.4,1.4,0.1)

'''
#distance to nearest neighbor vs. asymmetry
proper_distance_scale = cosmo.angular_diameter_distance(6).to('kpc') * (np.pi / 180)*(1/3600) # central_z would be the median redshift of your sample - z=6 is probably fine
fig = plt.figure(figsize=(4, 3))
ax = plt.gca()
# create a blank top axis in distance
distance_axis = np.linspace(np.min(protocluster_distance), 
                            np.max(protocluster_distance), 101) * proper_distance_scale 

ax2 = ax.twiny()


top_array = np.zeros(len(distance_axis))
ax2.plot(distance_axis, top_array, alpha=0)
ax2.set_xlabel("proper distance (kpc)", fontsize=11)

ax.scatter(protocluster_distance,protocluster_asymm,color='indigo',edgecolor='k',
           marker = '8', label='Protocluster')
ax.scatter(field_distance,field_asymm,color='violet',edgecolor='k',
           marker='s', label='Field')
ax.scatter(overdensity1_distance,overdensity1_asymm,color='mediumvioletred',
            edgecolor='k', marker='^', label='Overdensity 1')
ax.scatter(overdensity2_distance,overdensity2_rhalf,color='mediumpurple',
            edgecolor='k', marker='P', label='Overdensity 2')
#ax.scatter(galaxy_pair_distance,galaxy_pair_asymm,color='gold',
           #marker = '*',label='Mergers')
ax.set_xlabel('angular distance (arcsec)',fontsize=11)
ax.set_ylabel('Asymmetry',fontsize = 11)
ax.tick_params(axis='both',direction='in')
ax2.tick_params(axis='both',direction='in')
plt.tick_params(axis='both',direction='in')
ax.legend(loc='upper center',
          fontsize=8)
#plt.savefig('Distance_to_Nearest_Neighbor_Asymm.png',bbox_inches='')
plt.show()
'''
'''
plt.hist(asymmetry,edgecolor='black',bins=bins,color='blueviolet')
plt.xlabel('Asymmetry')
plt.ylabel('Number of Galaxies')
plt.suptitle('Galaxy Asymmetry')
plt.tick_params(axis='both',direction='in')
plt.show()


plt.hist(protocluster_asymm,edgecolor = 'black',bins=bins,color='plum')
plt.xlabel('Asymmetry')
plt.ylabel('Number of Galaxies')
plt.suptitle('Protocluster Asymmetry')
plt.tick_params(axis='both',direction='in')
plt.savefig('Protocluster_Asymmetry_hist.png')
plt.show()

plt.hist(field_asymm,edgecolor='black',bins=bins, color='violet')
plt.xlabel('Asymmetry')
plt.ylabel('Number of Galaxies')
plt.tick_params(axis='both',direction='in')
plt.suptitle('Field Asymmetry')
plt.savefig('Field_Asymmetry_hist.png')
plt.show()

plt.hist(overdensity1_asymm,edgecolor='black',bins=bins,color='purple')
plt.xlabel('Asymmetry')
plt.ylabel('Number of Galaxies')
plt.tick_params(axis='both',direction='in')
plt.suptitle('Overdensity 1 Asymmetry')
plt.savefig('Overdensity1_Asymmetry_hist.png')
plt.show()

plt.hist(overdensity2_asymm,edgecolor='black',bins=bins,color='indigo')
plt.xlabel('Asymmetry')
plt.ylabel('Number of Galaxies')
plt.tick_params(axis='both',direction='in')
plt.suptitle('Overdensity 2 Asymm')
plt.savefig('Overdensity2_Asymmetry_hist.png')
plt.show()

plt.hist(asymmetry,histtype = 'step',stacked=True,bins=bins, color='pink')
plt.hist(protocluster_asymm,label = 'protocluster',bins=bins,color='plum')
plt.hist(field_asymm,label='field',bins=bins,color='violet')
plt.hist(overdensity1_asymm,label='overdensity 1',bins=bins,color='purple')
plt.hist(overdensity2_asymm,label='overdensity 2',bins=bins,color='indigo')
plt.xlabel('Asymmetry')
plt.ylabel('Number of Galaxies')
plt.suptitle('Galaxy Asymmetry')
plt.tick_params(axis='both',direction='in')
plt.legend(loc='upper right')
plt.savefig('Layered_asymmetry_hist.png')
plt.show()
'''
#############################################################
#Plotting asymmetry v distance to nearest neighbor
'''
plt.scatter(protocluster_distance,protocluster_asymm,color='indigo')
plt.scatter(field_distance,field_asymm,color='violet')
plt.scatter(overdensity1_distance,overdensity1_asymm,color='mediumvioletred')
plt.scatter(overdensity1_distance,overdensity1_asymm,color='mediumpurple')
plt.savefig('Distance_Nearest_Neighbor_Asymm.png')
plt.legend(['protocluster','field','overdensity1','overdensity2'])
plt.xlabel('Angular Distance (arcsec)')
plt.ylabel('Asymmetry')
plt.tick_params(axis='both',direction='in')
plt.show
'''
#############################################################

#Plotting M20
M20 = data['M20']
protocluster_M20 =[]
field_M20=[]
overdensity1_M20= []
overdensity2_M20=[]
galaxy_pair_M20=[]
for i in protocluster_index:
    protocluster_M20.append(M20[i])
for i in field_index:
    field_M20.append(M20[i])
for i in overdensity1_index:
    overdensity1_M20.append(M20[i])
for i in overdensity2_index:
    overdensity2_M20.append(M20[i])
for i in galaxy_pair:
    galaxy_pair_M20.append(M20[i])
'''    
    
plt.hist(M20, edgecolor = 'black',color='blueviolet')
plt.xlabel('M20')
plt.ylabel('Number of Galaxies')
plt.suptitle('Galaxy M20')
plt.tick_params(axis='both',direction='in')
plt.savefig('Galaxy_M20_hist.png')
plt.show()

plt.hist(protocluster_M20, edgecolor = 'black',color='plum')
plt.xlabel('M20')
plt.ylabel('Number of Galaxies')
plt.suptitle('Protocluster M20')
plt.tick_params(axis='both',direction='in')
plt.savefig('Protocluster_M20.png')
plt.show()


plt.hist(field_M20, edgecolor = 'black',color='violet')
plt.xlabel('M20')
plt.ylabel('Number of Galaxies')
plt.suptitle('Field M20')
plt.tick_params(axis='both',direction='in')
plt.savefig('Field_M20_hist.png')
plt.show()


plt.hist(overdensity1_M20, edgecolor = 'black',color='purple')
plt.xlabel('M20')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 1 M20')
plt.tick_params(axis='both',direction='in')
plt.savefig('Overdenisty1_M20_hist.png')
plt.show()


plt.hist(overdensity2_M20, edgecolor = 'black',color='indigo')
plt.xlabel('M20')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 2 M20')
plt.tick_params(axis='both',direction='in')
plt.savefig('Field_M20_hist.png')
plt.show()

plt.hist(M20, fill=False,stacked=True,color='pink')
plt.hist(protocluster_M20, label = 'Protocluster',color='plum')
plt.hist(field_M20,label='Field',color='violet')
plt.hist(overdensity2_M20, label='Overdensity 1',color='purple')
plt.hist(overdensity1_M20,label = 'Overdensity 2',color='indigo')
plt.xlabel('M20')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity M20')
plt.legend(loc='upper right')
plt.tick_params(axis='both',direction='in')
plt.savefig('Layered_M20_hist.png')
plt.show()
'''


##############################################################

#Plotting Gini
Gini = data['Gini']
protocluster_Gini =[]
field_Gini=[]
overdensity1_Gini= []
overdensity2_Gini=[]
galaxy_pair_Gini=[]
for i in protocluster_index:
    protocluster_Gini.append(Gini[i])
for i in field_index:
    field_Gini.append(Gini[i])
for i in overdensity1_index:
    overdensity1_Gini.append(Gini[i])
for i in overdensity2_index:
    overdensity2_Gini.append(Gini[i])
for i in galaxy_pair:
    galaxy_pair_Gini.append(Gini[i])
'''
plt.hist(Gini, edgecolor = 'black',color='blueviolet')
plt.xlabel('Gini')
plt.ylabel('Number of Galaxies')
plt.suptitle('Galaxy Gini')
plt.tick_params(axis='both',direction='in')
plt.savefig('Galaxy_Gini_hist.png')
plt.show()

plt.hist(protocluster_Gini, edgecolor = 'black',color='plum')
plt.xlabel('Gini')
plt.ylabel('Number of Galaxies')
plt.suptitle('Protocluster Gini')
plt.tick_params(axis='both',direction='in')
plt.savefig('Protocluster_Gini.png')
plt.show()


plt.hist(field_Gini, edgecolor = 'black',color='violet')
plt.xlabel('Gini')
plt.ylabel('Number of Galaxies')
plt.suptitle('Field Gini')
plt.tick_params(axis='both',direction='in')
plt.savefig('Field_Gini_hist.png')
plt.show()


plt.hist(overdensity1_Gini, edgecolor = 'black',color='purple')
plt.xlabel('Gini')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 1 Gini')
plt.tick_params(axis='both',direction='in')
plt.savefig('Overdenisty1_Gini_hist.png')
plt.show()


plt.hist(overdensity2_Gini, edgecolor = 'black',color='indigo')
plt.xlabel('Gini')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity 2 Gini')
plt.tick_params(axis='both',direction='in')
plt.savefig('Field_Gini_hist.png')
plt.show()

plt.hist(Gini, fill=False,stacked=True,color='pink')
plt.hist(protocluster_Gini, label = 'Protocluster',color='plum')
plt.hist(field_Gini,label='Field',color='violet')
plt.hist(overdensity2_Gini, label='Overdensity 1',color='purple')
plt.hist(overdensity1_Gini,label = 'Overdensity 2',color='indigo')
plt.xlabel('Gini')
plt.ylabel('Number of Galaxies')
plt.suptitle('Overdensity Gini')
plt.tick_params(axis='both',direction='in')
plt.legend(loc='upper left')
plt.savefig('Layered_Gini_hist.png')
plt.show()
'''


######################################################################

#Plotting Gini v Asymmetry
'''
Gini_Y=[]
for x in asymmetry: #getting merger candidates
    G = -0.4 * x + 0.68
    Gini_Y.append(G)
    
    
plt.scatter(Gini,asymmetry)
#plt.scatter(galaxy_pair_asymm,galaxy_pair_Gini, marker ='*',label='Mergers')
plt.plot(asymmetry,Gini_Y,'--',color='grey',linewidth='0.5')
plt.xlabel('Asymmetry')
plt.ylable('Gini')
plt.tick_params(axis='both',direction='in')
plt.show()
'''

######################################################################
                                    
#Recreation of fig 3 of Constantine+24 Gini v. M20 using Seaborn



#Gini v. M20 Plot
'''
sns.kdeplot(data = data_set, x='M20', y='Gini', fill=True,color='violet').set(title = 'Gini v. M20')
sns.color_palette("flare", as_cmap=True)
sns.scatterplot(data=data_set,x='M20',y='Gini',color='darkviolet')
plt.plot(M20_x,Gini_Y,'--',color='black',linewidth='0.5')
#plt.savefig('Gini_M20_kde.png')
plt.tick_params(axis='both',direction='in')
plt.show()

#sns.kdeplot(data= data_set,x='M20',y='Gini',hue='ZO3_2D',fill=True).set(title = 'Gini v. M20')
#sns.scatterplot(data=data_set,x='M20',y='Gini',hue='ZO3_2D')
#plt.plot(M20,Gini_Y,'--',color='gray',linewidth='0.5')
#plt.show()

'''



######################################################################

#Calculating a plotting stellar mass surface density sigma Mstar


SED_properties = ascii.read('OIII_SED_properties_022624_delayedtau.dat')
SED_ra=SED_properties['ra']
SED_dec=SED_properties['dec']
SED_Mstars = SED_properties['logMstar_delayedtau']
Mstars=np.zeros(124,dtype='float')
Sigma_Mstars=[]
#getting the 
for index,(ra,dec) in enumerate(zip(XX,YY)):
    for gal,(RA,DEC) in enumerate(zip(SED_ra,SED_dec)): #checks if the galaxies are the same through coords
        if ra==RA and dec==DEC:
            Mstars[index]=(SED_Mstars[gal]) #adds the MSTAR to the list at corresponding galaxy index
        
for Mstar,Reff in zip(Mstars,rhalf):
    Sigma_Mstar = Mstar / (2 * np.pi * Reff)
    Sigma_Mstars.append(Sigma_Mstar)


#Calculating Overdensity rhalf
protocluster_rhalf =[]
field_rhalf=[]
overdensity1_rhalf= []
overdensity2_rhalf=[]
for i in protocluster_index:
    protocluster_rhalf.append(rhalf[i])
for i in field_index:
    field_rhalf.append(rhalf[i])
for i in overdensity1_index:
    overdensity1_rhalf.append(rhalf[i])
for i in overdensity2_index:
    overdensity2_rhalf.append(rhalf[i])

    
#Calculating Overdensity Mstar
protocluster_mstars =[]
field_mstars=[]
overdensity1_mstars= []
overdensity2_mstars=[]
galaxy_pair_mstars=[]
for i in protocluster_index:
    protocluster_mstars.append(Mstars[i])
for i in field_index:
    field_mstars.append(Mstars[i])
for i in overdensity1_index:
    overdensity1_mstars.append(Mstars[i])
for i in overdensity2_index:
    overdensity2_mstars.append(Mstars[i])

# Calculating Overdensity sigma Mstars
protocluster_sig_mstars =[]

for Mstar,Reff in zip(protocluster_mstars,protocluster_rhalf):
    Sigma_Mstar = Mstar / (2 * np.pi * Reff)
    protocluster_sig_mstars.append(Sigma_Mstar)

field_sig_mstars=[]
for Mstar,Reff in zip(field_mstars,field_rhalf):
    Sigma_Mstar = Mstar / (2 * np.pi * Reff)
    field_sig_mstars.append(Sigma_Mstar)
    
overdensity1_sig_mstars= []
for Mstar,Reff in zip(overdensity1_mstars,overdensity1_rhalf):
    Sigma_Mstar = Mstar / (2 * np.pi * Reff)
    overdensity1_sig_mstars.append(Sigma_Mstar)
overdensity2_sig_mstars=[]
for Mstar,Reff in zip(overdensity2_mstars,overdensity2_rhalf):
    
    
#plotting M star and Sig Mstar as histogram
    Sigma_Mstar = Mstar / (2 * np.pi * Reff)
    overdensity2_sig_mstars.append(Sigma_Mstar)

###############################################################################
#Getting stellar mass density of each galaxy and its merger pair

galaxy_Mstar=[]
pair_Mstar=[]
merger_Mstar_ratio=[]
protocluster_merger_ratio=[]
protocluster_merger_redshift=[]
field_merger_ratio=[]
field_merger_redshift=[]
overdensity1_merger_ratio=[]
overdensity1_merger_redshift=[]
overdensity2_merger_ratio=[]
overdensity2_merger_redshift=[]


for i,j in enumerate(galaxy_pair):
    if isinstance(i,np.ndarray):
        i= i[0]
    if i % 2==0:
        print(i)
        galaxy_Mstar.append(10**(Mstars[i]))
       
    else:
        pair_Mstar.append(10**(Mstars[i]))
        
for i,j in zip(galaxy_Mstar,pair_Mstar):#array of stellar mass density ration Mgal/Mpair
    merger_Mstar_ratio.append(float(i)/float(j))
  
    

for i,j in zip(galaxy_pair_zmed,merger_Mstar_ratio):
    if 6.5 < i :
        protocluster_merger_ratio.append(j)
        protocluster_merger_redshift.append(i)
    elif 5.35 < i < 5.45:
        overdensity1_merger_ratio.append(j)
        overdensity1_merger_redshift.append(i)
    elif  6.15 < i < 6.3:
        overdensity2_merger_ratio.append(j)
        overdensity2_merger_redshift.append(i)
    else:
        field_merger_ratio.append(j)
        field_merger_redshift.append(i)
      

#Plot of Merger Stellar Mass Ratio vs. Pair Redshift (Separate Y axis)
'''
fig = plt.figure()
ax = plt.subplot(111)
#ax = plt.gca()
ax2 = ax.twinx() 
ax2.set_ylabel('Number of Galaxies',fontsize=15)
bins = np.arange(5.2,7.1,0.1)
ax2.hist(all_galaxy_shift,histtype='step',bins=bins,fill=False,color='indigo',)
ax2.hist(protocluster, label='protocluster',bins=bins,color='plum',alpha=0.5)
ax2.hist(overdensity1, label='overdensity1',bins=bins,color='mediumslateblue',alpha=0.5)
ax2.hist(overdensity2, label='overdensity2',bins=bins,color='mediumvioletred',alpha=0.5)


ax.scatter(protocluster_merger_redshift,protocluster_merger_ratio,color='plum',
            edgecolor='black')
ax.scatter(field_merger_redshift,field_merger_ratio,color='indigo',
            edgecolor = 'black')
ax.scatter(overdensity1_merger_redshift,overdensity1_merger_ratio,color='mediumslateblue',
            edgecolor='black')
ax.scatter(overdensity2_merger_redshift,overdensity2_merger_ratio,color='mediumvioletred',alpha=0.7,
            edgecolor='black')
#plt.plot(asymmetry,Gini_Y,'--',color='grey',linewidth='0.5')


#ax.scatter(galaxy_pair_zmed,merger_Mstar_ratio,color='indigo')
ax.set_xlabel('Redshift',fontsize=15)
ax.set_ylabel('Merger $M_{star}$ Ratio',fontsize=15)
ax.set_yscale('log')
ax.legend(['protocluster','field','overdensity 1','overdensity2'],
        bbox_to_anchor=(0.01,0.7)  ,loc='lower left',fontsize=12)
minor_merger_ymin=[0.33]*124
X = np.linspace(5.25,7.1,num=124)
#ax.plot(X,minor_merger_ymin,'--',color='gray',linewidth='0.5')
ax2.tick_params(axis='both',direction='in',labelsize=14)
ax.tick_params(axis='both',direction='in',which='both',labelsize=14)
ax.axhline(0.33,color='lightgray',linestyle='--')

# add text to the plot
ax.text(6,4,'Minor Mergers', fontsize=14,color='gray',fontweight=550)
ax.text(5.65,0.08,'Major Mergers',fontsize=14,color='gray',fontweight=550)
#plt.savefig('Merger_Mstar_Ratio.png')
plt.show()


'''


#Histogram of Stellar Mass Density

'''

plt.hist(Sigma_Mstars,stacked = True, fill = False, edgecolor='pink')
plt.hist(protocluster_sig_mstars,label='protocluster',color='plum')
plt.hist(field_sig_mstars,label='field',color='violet')
plt.hist(overdensity1_sig_mstars,label='overdensity 1',color='purple')
plt.hist(overdensity2_sig_mstars, label='overdensity 2',color='indigo')

plt.xlabel('$\Sigma_{M_{star}}$')
plt.ylabel ('Number of Galaxies')

plt.legend(loc='upper right')
plt.suptitle('Galaxy Stellar Mass Density')
plt.tick_params(axis='both',direction='in')
plt.savefig('Galaxy_Sigma_Mstar.png')
plt.show()




plt.hist(protocluster_sig_mstars,edgecolor = 'black',color='plum')
plt.xlabel('$\Sigma_{M_star}$')
plt.ylabel ('Number of Galaxies')
plt.suptitle('Protocluster Stellar Mass Density')
plt.tick_params(axis='both',direction='in')
plt.show()
plt.savefig('Protocluster Sigma Mstar')

plt.hist(field_sig_mstars,edgecolor = 'black',color='violet')
plt.xlabel('$\Sigma_{M_star}$')
plt.ylabel ('Number of Galaxies')
plt.suptitle('Field Stellar Mass Density')
plt.tick_params(axis='both',direction='in')
plt.show()
#plt.savefig('Field Sigma Mstar')

plt.hist(overdensity1_sig_mstars,edgecolor = 'black',color='purple')
plt.xlabel('$\Sigma_{M_star}$')
plt.ylabel ('Number of Galaxies')
plt.suptitle('Overdensity 1 Stellar Mass Density')
plt.tick_params(axis='both',direction='in')
plt.show()
plt.savefig('Overdensity 1 Sigma Mstar')

plt.hist(overdensity2_sig_mstars,edgecolor = 'black',color='indigo')
plt.xlabel('$\Sigma_{M_star}$')
plt.ylabel ('Number of Galaxies')
plt.suptitle('Overdensity 2 Stellar Mass Density')
plt.tick_params(axis='both',direction='in')
plt.show()
plt.savefig('Overdensity 2 Sigma Mstar')
'''




###############################################################################

#Plotting M2O v Gini Scatter plot with concentration colorbar



concentration = data['Concentration'] 
protocluster_concentration =[]
field_concentration=[]
overdensity1_concentration= []
overdensity2_concentration=[]
for i in protocluster_index:
    protocluster_concentration.append(concentration[i])
for i in field_index:
    field_concentration.append(concentration[i])
for i in overdensity1_index:
    overdensity1_concentration.append(concentration[i])
for i in overdensity2_index:
    overdensity2_concentration.append(concentration[i])



#Scatter plot of all Gini v M20 separated by protocluster with concentration map

Gini_Y=[]
M20_x= np.linspace(-2.25,-0.4,101)
for x in M20_x: #getting merger candidates
    G = -0.14 * x + 0.33  #Lotz, J. M., Davis, M., Faber, S. M., et al. 2008, ApJ, 672, 177
    Gini_Y.append(G)


#Calculating Gini M20 merger pair fraction
Gini_M20_mergers=[] #Gini of mergers to find merger fraction
Gini_M20_merger_redshift=[]
Gini_M20_merger_frac= 0
protocluster_Gini_M20=[]
field_Gini_M20=[]
overdensity1_Gini_M20=[]
overdensity2_Gini_M20=[]

for index,(i,j) in enumerate(zip(Gini,M20)):
    if i > -0.14 * j + 0.33:
        Gini_M20_merger_frac += 1
        Gini_M20_mergers.append(index)
        
        Gini_M20_merger_redshift.append(all_galaxy_shift[index]) #storing index of merger
'''
for i in Gini_M20_Merger_Index:
    if i in protocluster_index:
        protocluster_Gini_M20.append(i)
    elif i in field_index:
        field_Gini_M20.append(i)
    elif i in overdensity1_index:
        overdensity1_Gini_M20.append(i)
    else:
        overdensity2_Gini_M20.append(i)
    
'''
Gini_M20_merger_frac = Gini_M20_merger_frac/len(all_galaxy_shift)
Gini_M20_merger_frac_error = np.sqrt(len(Gini_M20_mergers))/len(all_galaxy_shift)

# Finding overlapping galaxies in both GiniM20 and close pair rate criteria
total_mergers=[]
for i in Gini_M20_mergers:
    if i in galaxy_pair:
        total_mergers.append(i)
        

'''
#protocluster_Gini_M20_frac=len(protocluster_Gini_M20)/len(Gini_mergers)
#protocluster_Gini_M20_error= np.sqrt(len(protocluster_Gini))/len(Gini_mergers)
#field_Gini_M20_frac=len(field_Gini_M20)/len(Gini_mergers)
#overdensity1_Gini_M20_frac=len(overdensity1_Gini_M20)/len(Gini_mergers)
#overdensity2_Gini_M20_frac=len(overdensity2_Gini_M20)/len(Gini_mergers)







fig, ax = plt.subplots()
sca=ax.scatter(field_M20, field_Gini,c=field_concentration, marker = '8', cmap='RdPu',edgecolor='black')
sca=ax.scatter(protocluster_M20, protocluster_Gini,c = protocluster_concentration,marker = 's', cmap='RdPu',edgecolor='black')
sca = ax.scatter(overdensity1_M20,overdensity1_Gini, c = overdensity1_concentration, marker ='^',cmap='RdPu',edgecolor='black')
sca = ax.scatter(overdensity2_M20,overdensity2_Gini, c = overdensity2_concentration, marker ='P',cmap='RdPu',edgecolor='black')

plt.colorbar(sca).set_label(label='Concentration', size=14)
plt.plot(M20_x,Gini_Y,'--',color='black',linewidth='0.5')
plt.xlabel('M20',fontsize=14)
plt.ylabel('Gini',fontsize=14)
#plt.suptitle('Gini v M20')
plt.legend(['field','protocluster','overdensity 1','overdensity 2'])
plt.tick_params(axis='both',direction='in',labelsize=13)

# add text to a plot
ax.text(-1.48, 0.35,'Gini > -0.14 * M20 + 0.33', fontsize=11,color='gray',fontweight=550)
ax.text(-2.0,0.63,'Mergers',fontsize=10,color='gray',fontweight=550)


#plt.savefig('Gini_M20_Concentration.png',bbox_inches='tight',dpi=200,pad_inches=0.1)


plt.show()
'''



#plotting Gini and M20 as a function of Sigma Mstar
'''
#plt.scatter(Sigma_Mstars,Gini)
plt.scatter(protocluster_sig_mstars,protocluster_Gini,label='protocluster',color='plum')
plt.scatter(field_sig_mstars,field_Gini,label='field',color='violet')
plt.scatter(overdensity1_sig_mstars,overdensity1_Gini,label='overdensity 1',color='purple')
plt.scatter(overdensity2_sig_mstars,overdensity2_Gini,label='overdensity 2',color='indigo')

plt.xlabel('$\Sigma_{M_{star}}$')
plt.ylabel('Gini')
plt.suptitle('Gini v Mstar')
plt.tick_params(axis='both',direction='in')
plt.legend(loc='upper right')
plt.savefig('Gini_Mstar.png')
plt.show()

M20 = data['M20']
#plt.scatter(Sigma_Mstars,)
plt.scatter(protocluster_sig_mstars,protocluster_M20,label='protocluster',color='plum')
plt.scatter(field_sig_mstars,field_M20,label='field',color='violet')
plt.scatter(overdensity1_sig_mstars,overdensity1_M20,label='overdensity 1',color='purple')
plt.scatter(overdensity2_sig_mstars,overdensity2_M20,label='overdensity 2',color='mediumslateblue')
plt.xlabel('$\Sigma_{M_{star}}$')
plt.ylabel('M20')
plt.suptitle('M20 v MStar')
plt.tick_params(axis='both',direction='in')
plt.legend(loc='upper right')
plt.savefig('Mstar_M20.png')
plt.show()
'''

#Plotting Gini  as a funciton of half radius
'''
plt.scatter(protocluster_rhalf,protocluster_Gini,color ='indigo',marker ='s')
plt.scatter(field_rhalf,field_Gini,color='violet',marker='8')
plt.scatter(overdensity1_rhalf,overdensity1_Gini,color='mediumvioletred',marker='^')
plt.scatter(overdensity2_rhalf,overdensity2_Gini,color='mediumpurple',marker='P')
plt.xlabel('Half Radius')
plt.ylabel('Gini')
plt.tick_params(axis='both',direction='in')
plt.suptitle('Gini v Half Radius')
plt.legend(['Protocluster','Field','Overdensity 1','Overdensity 2'])
plt.savefig('Gini v. Half Radius.png')
plt.show()
'''

#Plotting  M20 as a funciton of half radius
'''
plt.scatter(protocluster_rhalf,protocluster_M20,color ='indigo',marker ='s')
plt.scatter(field_rhalf,field_M20,color='violet',marker='8')
plt.scatter(overdensity1_rhalf,overdensity1_M20,color='mediumvioletred',marker='^')
plt.scatter(overdensity2_rhalf,overdensity2_M20,color='mediumpurple',marker='P')
plt.xlabel('Half Radius')
plt.ylabel('M20')
plt.tick_params(axis='both',direction='in')
plt.suptitle('M20 v Half Radius')
plt.legend(['Protocluster','Field','Overdensity 1','Overdensity 2'])
plt.savefig('M20 V Half Radius.png')
plt.show()
'''
###############################################################################

#Creating Aplypy figure
'''
import aplpy # you may need to pip install aplpy
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM


cosmo = FlatLambdaCDM(H0=70 * u.km/u.s/u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

sources = Table.read('jw_o3_candidates_VI_v20240115v2.fits', format='fits')
good_sources = sources[np.where(sources['Score_VI']==3)]
good_sources


source_113 = good_sources[113]
PAIR_SOURCE = good_sources[112] #galaxy merging pair 
source_113_coords = SkyCoord(ra=source_113['ra'], dec=source_113['dec'], unit='deg')
PAIR_COORDS = SkyCoord(ra=PAIR_SOURCE['ra'], dec=PAIR_SOURCE['dec'], unit='deg')
print(PAIR_COORDS.ra,PAIR_COORDS.dec)
print(source_113_coords.ra,source_113_coords.dec)
proper_distance = cosmo.angular_diameter_distance(5.3).to('kpc') * np.pi/180 * 1/3600

#fig = plt.figure(figsize=(5,5))

obj1 = 'galaxy113_F356W_psfmatch_overdensity2.fits'
obj2 = 'galaxy112_F356W_psfmatch_overdensity2.fits'
wcs = WCS(fits.open(obj1)[0].header)

# initiate the fitsfigure 
JWSTobj = aplpy.FITSFigure(obj1)
#JWSTobj2 = aplpy.FITSFigure(obj2)

#
mean_ra = (source_113_coords.ra.degree+PAIR_COORDS.ra.degree)/2
mean_dec=(source_113_coords.dec.degree+PAIR_COORDS.dec.degree)/2
# this does allow you to zoom in or out. the arguments are center_ra, center_dec, and size
# in degrees - 1.5/3600 = 1.5 arcsec
#JWSTobj.recenter(source_113_coords.ra.degree, source_113_coords.dec.degree, 2/3600)
JWSTobj.recenter(mean_ra, mean_dec, 1.25/3600)

#JWSTobj.recenter(mean_ra, mean_dec, 1.5/3600)

# you can do grayscale or colorscale, see documentation for changing the colormap & scaling
#JWSTobj.show_grayscale(invert=True)
JWSTobj.show_colorscale(cmap='RdPu')

# some aesthetic stuff
ra_ax = JWSTobj.ax.coords[0]
dec_ax = JWSTobj.ax.coords[1]
ra_ax.set_major_formatter('hh:mm:ss')
dec_ax.set_major_formatter('dd:mm:ss')
JWSTobj.ticks.set_xspacing(0.2/3600)
JWSTobj.ticks.set_yspacing(0.2/3600)
JWSTobj.ticks.set_length(10)
ra_ax.set_ticks(direction='in', color='black')
dec_ax.set_ticks(direction='in', color='black')

JWSTobj.axis_labels.hide_x()
JWSTobj.axis_labels.hide_y()

JWSTobj.tick_labels.hide_x()
JWSTobj.tick_labels.hide_y()

# label the galaxy centers
JWSTobj.show_markers(source_113_coords.ra, source_113_coords.dec,
                     edgecolor='darkblue', facecolor='none', marker='^', linewidth=1.5, s=200)
#  coordinates of the pair
JWSTobj.show_markers(PAIR_COORDS.ra, PAIR_COORDS.dec,
                     edgecolor='cyan', facecolor='none', marker='*', linewidth=1.5, s=200)

# add scalebar
JWSTobj.add_scalebar(1/3600)
JWSTobj.scalebar.set_label("1'' (6 kpc)")

plt.plot(10000, 10000, '^', color='darkblue', fillstyle='none', label='[OIII] emitter')
plt.plot(10000, 10000, '*', color='cyan', fillstyle='none', label='Merging Pair')
plt.legend(loc=3, ncol=1, fontsize=13)



#plt.savefig("aplpy_figure2.png", dpi=200,  bbox_inches='tight', 
            #pad_inches=0.05)
'''
###############################################################
#Getting  star formation rates

full_SFR_list = SED_properties['SFR_delayedtau']
SFR = np.zeros((124),dtype='float')
#calculate distance bewteen ra and dec compare that instead
#getting correct SFR based on galaxy coords
for index,(ra,dec) in enumerate(zip(XX,YY)):
    for gal,(RA,DEC) in enumerate(zip(SED_ra,SED_dec)): #checks if the galaxies are the same through coords
        if SkyCoord(ra=ra, dec=dec, unit=u.degree)==SkyCoord(ra=RA, dec=DEC, unit=u.degree):
            SFR[index]=(full_SFR_list[gal]) #adds the SFR to the list at corresponding galaxy index
        
#getting median SFR rates for all mergers and non merger
mergers_SFR=[]
nonmergers_SFR=[]
for i,j in enumerate(SFR):
    if i in Gini_M20_mergers: #finding SFR for gini m20 merger criteria
        mergers_SFR.append(j)
    else:
        nonmergers_SFR.append(j)
median_merger_SFR= np.median(mergers_SFR)      
median_nonmerger_SFR= np.median(nonmergers_SFR) 
           

 
            
            
            