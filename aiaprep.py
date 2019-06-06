#Package Imports-----------------------------------------------------------------------------------------------------------------------------------------
import pandas as pd #databases and Series
import numpy as np 

import matplotlib.pyplot as plt
import matplotlib.colors as colors

#AIA PACKAGES
import sunpy.map #this is the sub-module that holds \'Map\'
from sunpy.map import Map #Map is a spatially-aware 2D array
from sunpy.map.sources.sdo import AIAMap, HMIMap

from sunpy.net import Fido, attrs as a #fido is an objext in sunpy.net - searching & fetching data

import drms #python package that can be used to access HMI/AIA/MDI data

#ASTRONOMY PACKAGES
import astropy.units as u #handles defining/conversion of units
                          #e.g. lambda=171*u.angstrom
                          #lambda.value and lambda.unit can be printed
from astropy.coordinates import SkyCoord
import astropy.wcs
from astropy.wcs import WCS

#MISCELLANEOUS
import datetime
from datetime import timedelta #to support the subtraction and addition of time

import time
from inspect import signature 
import copy
import scipy.constants
import scipy.interpolate as interpolate

import sunpy.physics.differential_rotation
from sunpy.physics.differential_rotation import diffrot_map
#----------------------------------------------------------------------------------------------------------------------------------------------------------
def aia_prep(aiamap,mpo_3Hours=False,mpo_1Week=False,cutout=False): #aiamap is the only positional parameter passed into the code
                                                                    #Additional keyword parameters can be added to perform curated tasks for the user
    """
    This code performs image registration for one or more inputted AIAMaps, returning the updated AIAMaps as output
    
    Calling options (Copy-Paste into command line): 
    1. aia_prep(aiamap, mpo_3Hours=True) #Full-disk image(s) provided, use latest pointing data updated every 3 days
    2. aia_prep(aiamap, mpo_1Week=True) #Full-disk image(s) provided, use pointing data updated weekly
    3. aia_prep(aiamap)                 #Full-disk image(s) provided, 3 day updated pointing data will be used
    4. aia_prep(aiamap, mpo_3Hours=True, cutout=True)#cutout provided, use latest pointing data updated every 3 days
    5. aia_prep(aiamap, mpo_1Week=True, cutout=True)#cutout provided, use pointing data updated weekly
    6. aia_prep(aiamap, cutout=True) #cutout provided, 3 day updated pointing data will be used
    
    Inputs:
    aiamap: Map object containing the 2D spatially-aware data array + metadata (index)
            Mandatory, if not passed the function will raise an error
    
    mpo_1Week: Would you like to use the 1 Week Master Pointing List? If yes, mark as True
    mpo_3Hours: Would you like to use the 3 Day Master Pointing List? If yes, mark as True
    
    Output: Level 1.5 Map/List of Level 1.5 Maps
    """    
    #Step 1: Check to see if the passed file/files are compatible with the tasks aiaprep performs
    sig=signature(aia_prep) #contains the parameters and arguments passed into the function aia_prep
    params=sig.parameters  #writes out the parwillameters and their associated values here    

    if(len(params) < 1):    #making sure something is passed in with the function
        raise ValueError("Not enough parameters passed into the function, kindly provide atleast one AIAMap")
    if not isinstance(aiamap[0], (AIAMap, HMIMap)):
        raise ValueError("Input must be an AIAMap or HMIMap.")
    aiamap_copy=copy.deepcopy(aiamap)
    
    #Step 2: Iterating through each AIAMap & calling associated functions for aia_prep
    n_img=np.size(aiamap_copy) #aiamap is a list
    print('Number of aiamaps are {0}'.format(n_img))
    maps_updated=[] #OUTPUT
    for i in np.arange(0,n_img):
        
        #Step 3: Calling the master pointing list
        updated_header_jsoc=sdo_master_pointing(aiamap_copy[i],mpo_1Week, mpo_3Hours, cutout) #calling function to get pointing info
        print('-----------------------------------------')
        print('The information that needs to be updated in the meta data for aiamap # {0} are listed below:'.format(i+1))
        print(updated_header_jsoc) #this will return only one set of pointing pars (but in list form), since one aiamap was sent
        print('-----------------------------------------')
        
        #Step 4: Updating the header information using the latest pointing data acquired from JSOC
        map_updated1=header_update(aiamap_copy[i],updated_header_jsoc,mpo_1Week,mpo_3Hours,cutout) #calling function that will update the header file
        
        #Step 5: Create a reference header file for registration
        header_reference=copy.deepcopy(map_updated1[0].meta) #creating a copy of the header information
        ref={'crpix1':2048.5,'crpix2':2048.5,'naxis1':4096,'naxis2':4096,'crota2':0,'crval1':0,'crval2':0,'ctype1':'HPLN-TAN','ctype2':'HPLT-TAN','xcen':0.0,'ycen':0.0,'cdelt1':0.6,'cdelt2':0.6}
        header_reference['crpix1'] = 2048.5
        header_reference['crpix2'] = 2048.5
        header_reference['naxis1'] = 4096
        header_reference['naxis2'] = 4096
        header_reference['crval1'] = 0.0
        header_reference['crval2'] = 0.0
        header_reference['xcen']   = 0.0
        header_reference['ycen']   = 0.0
        header_reference['crota2'] = 0.0
        
        #Step 6: Image Registration
        """ This whole section is a mess and needs to be moved to drot_map & updated.
        smap= map_updated1[0] # Sunpy Map
        #dmap = drot_map(smap, header_reference) # rotate the map forward 4 days as seen from Earth

        dmap.peek() # Show the map
    
        deg_to_rad=scipy.constants.pi/180.0 #degrees to radians
        img = map_updated1[0].data #data from sunpy map
        x = (np.arange(1,header_reference['naxis1']+1) - header_reference['crpix1'])*header_reference['cdelt1'] + header_reference['crval1']
        y = (np.arange(1,header_reference['naxis2']+1) - header_reference['crpix2'])*header_reference['cdelt2'] + header_reference['crval2']
        z = np.zeros(shape=(len(x),len(y)),dtype=int)
        func = interpolation.interp(x[:,0],y[0,:],img)
        y, x = np.meshgrid(y,x)
        
        xr = np.cos(CROTA2*DTOR)*x - np.sin(CROTA2*DTOR)*y #rotated axes
        yr =-np.sin(CROTA2*DTOR)*x + np.cos(CROTA2*DTOR)*y #rotated axes
        image_updated=func.ev(xr,yr).reshape(header_reference['naxis1'],header_reference['naxis2'],order='fortran')
                                          #ev evaluates the spline at points, returns interpolated values
        
        #Step 4: Scaling the map
        if (map_updated1[0].scale[0] / 0.6).round() != 1.0 * u.arcsec and map_updated1[0].data.shape != (4096, 4096): #use [0] for map_updated1 since it is returned as a list
            scale = (map_updated1[0].scale[0] / 0.6).round() * 0.6 * u.arcsec
        else:
            scale = 0.6 * u.arcsec 
        scale_factor = map_updated1[0].scale[0] / scale
        
        #Step 5: Rotating the map and finding new range
        tempmap = map_updated1[0].rotate(recenter=True, scale=scale_factor.value,missing=map_updated1[0].min())
        center = np.floor(tempmap.meta['crpix1'])
        range_side = (center + np.array([-1, 1]) * tempmap.data.shape[0] / 2) * u.pix
        
        #Step 6: Creating the final version of the aia level 1.5 map
        newmap = tempmap.submap(u.Quantity([range_side[0], range_side[0]]),u.Quantity([range_side[1], range_side[1]]))
        newmap.meta['r_sun'] = newmap.meta['rsun_obs'] / newmap.meta['cdelt1']
        newmap.meta['bitpix'] = -64 #bits/pixel, negative for floating points
        """
        
        maps_updated.append(dmap) #appending to list. On to the next map!
                                
    return maps_updated
    
#----------------------------------------------------------------------------------------------------------------------
def header_update(aiamap,updated_header_jsoc,mpo_1Week=False,mpo_3Hours=False,cutout=False): 
    
    """
    The header information of the AIAMap are updated using the latest JSOC Pointing Information, concluding the 
    first step of AIAMap update from Level 1 to 1.5.
    Inputs of 1. aiamap(s) 2. the updated header information (record number, sun center coordinates, 
    instrument rotation, image scale in arcsec/pixel of CCD) from sdo_master_pointing.py, and 3. whether image
    is full disk or cutout.
    
    Calling options (Copy-Paste into command line): 
    1. header_update(aiamap,updated_header_jsoc) #Full-disk image(s) provided along with jsoc pointing data
    2. header_update(aiamap,updated_header_jsoc,cutout=False) #Full-disk image(s) provided along with jsoc pointing data
    3. header_update(aiamap,updated_header_jsoc,cutout=True) #cutout image(s) provided along with jsoc pointing data    
    """
    print('Header Update Begins Now')
    #--------------------------------------------------------------------------------------------------------------------    
    #Step 1: Ensuring the AIAMap(s) input is a list, when more than one AIAMap is passed it will be a list, but 
             #if only 1 AIAMap is passed, it will have to be converted to a list of 1. 

    aiamap_list=[]
    aiamap_copy=copy.deepcopy(aiamap)
    if isinstance(aiamap_copy,__builtins__.list)!=True: #in order to loop through one or more AIAMaps, they have to be in list form
        aiamap_list.append(aiamap_copy) #if it's not a list, make it a list (only happens when it's just 1 AIAMap)
    if isinstance(aiamap_copy,__builtins__.list)==True:
        aiamap_list=aiamap_copy 
    #--------------------------------------------------------------------------------------------------------------------    
    #Step 2: Checking Inputs
    sig=signature(header_update) #contains the parameters and arguments passed into the function header_update
    params=sig.parameters  #writes out the parameters and their associated values here    
    #print(params) #print(len(params))
    if(len(params) < 2):    #making sure aiamap and jsoc pointing data are passed
        raise ValueError("Not enough parameters passed into the function, kindly provide atleast one AIAMap and the JSOC Pointing Info file")
    if not isinstance(aiamap_list[0], (AIAMap, HMIMap)):
        raise ValueError("Input must be an AIAMap or HMIMap.")
    #The code will automatically raise an error if the Pointing data isn't produced form sdo_master_pointing    
    
    if(len(aiamap_list)!=len(updated_header_jsoc)):
        raise ValueError("The number of AIAMaps and corresponding JSOC Master Pointing lists donot match")
    
    n_img=len(aiamap_list) #number of AIAMaps
    maps_updatedheader=[] #OUTPUT
    for f in np.arange(0,n_img):       
        #--------------------------------------------------------------------------------------------------------------------    
        #Step 3: Using AIA Map Information used to construct custom JSOC Keyword names        
        wavelnth=aiamap_list[f].wavelength #astropy.units.quantity.Quantity, value+units
        wavelnth_num=wavelnth.value #Just the numerical portion of the wavelength
        x0_pointing='A_'+str(int(wavelnth_num))+'_X0'
        y0_pointing='A_'+str(int(wavelnth_num))+'_Y0'
        instrot_pointing='A_'+str(int(wavelnth_num))+'_INSTROT'
        imscale_pointing='A_'+str(int(wavelnth_num))+'_IMSCALE'
        recnum='*recnum*'
    
        #Step 4: updating the header data of the aiamap passed to this code
        aiamap_list[f].meta['lvl_num']=1.5
        aiamap_list[f].meta['inst_rot']=float(updated_header_jsoc[f][instrot_pointing].values)
        aiamap_list[f].meta['imscl_mp']=float(updated_header_jsoc[f][imscale_pointing].values)
        aiamap_list[f].meta['x0_mp']=float(updated_header_jsoc[f][x0_pointing].values)
        aiamap_list[f].meta['y0_mp']=float(updated_header_jsoc[f][y0_pointing].values)
    
        #Step 5: Updating the MPO Record Number based on the MPO requests by the user
        if mpo_1Week == True:
            aiamap_list[f].meta['mpo_rec']='MPO_1Week_#'+str(updated_header_jsoc[f][recnum].values)
        if mpo_3Hours == True:
            aiamap_list[f].meta['mpo_rec']='MPO_3Hours_#'+str(updated_header_jsoc[f][recnum].values)
        if mpo_1Week == False and mpo_3Hours == False:
            aiamap_list[f].meta['mpo_rec']='MPO_3Hours_#'+str(updated_header_jsoc[f][recnum].values)        
    
        #Step 6: updating certain values based on updated JSOC values
        aiamap_list[f].meta['crota2'] = float(aiamap_list[f].meta['sat_rot'] + updated_header_jsoc[f][instrot_pointing].values) #sat_rot + inst_rot, rotation for array axes to get to image axes (deg)
        aiamap_list[f].meta['crpix1'] = float(updated_header_jsoc[f][x0_pointing].values + 1.) #center of image
        aiamap_list[f].meta['crpix2'] = float(updated_header_jsoc[f][y0_pointing].values + 1.) 
        aiamap_list[f].meta['cdelt1'] = float(updated_header_jsoc[f][imscale_pointing].values) #plate scale
        aiamap_list[f].meta['cdelt2'] = float(updated_header_jsoc[f][imscale_pointing].values)
        #add crval1 and 2 in iteration 2
        
        #Step 7: updating xcen and ycen. Equations sourced from comp_fits_cen [from nomenclature of old solar soft]
        aiamap_list[f].meta['xcen'] = float(aiamap_list[f].meta['crval1'] + aiamap_list[f].meta['cdelt1']*((float(aiamap_list[f].meta['naxis1'])+1.0)/2.0 - aiamap_list[f].meta['crpix1']))
        aiamap_list[f].meta['ycen'] = float(aiamap_list[f].meta['crval2'] + aiamap_list[f].meta['cdelt2']*((float(aiamap_list[f].meta['naxis2'])+1.0)/2.0 - aiamap_list[f].meta['crpix2']))   
    
    #Step 8: Packaging and returning the AIAMaps  
    for g in np.arange(0,n_img):
        map_head=sunpy.map.Map((aiamap_list[g].data, aiamap_list[g].meta))
        maps_updatedheader.append(map_head)
    
    return maps_updatedheader #returning the updated aiamap(s)
#----------------------------------------------------------------------------------------------------------------------
def sdo_master_pointing(aiamap,mpo_1Week=False,mpo_3Hours=False,cutout=False): #or set it to None for default
    """
    This function accesses the master pointing list to update the header of the level 1 data. The following parameters
    are returned: record number, sun center coordinates, instrument rotation, image scale in arcsec/pixel of CCD.
    This is a generated structure contains reference values of WCS keywords, that are to be used as a 
    reference for aligning an SDO image, either full frame or partial frame
    
    Calling options (Copy-Paste into command line): 
    1. sdo_master_pointing(aiamap,mpo_3Hours=True,cutout=False) #Full-disk image(s) provided, 
                                                                   use latest pointing data updated every 3 hourss
    2. sdo_master_pointing(aiamap,mpo_1Week=True,cutout=False) #Full-disk image(s) provided, 
                                                                    use 1 Week pointing data
    3. sdo_master_pointing(aiamap) #Full-disk image(s) provided, use latest pointing data updated every 3 hours
    
    4. sdo_master_pointing(aiamap, mpo_3Hours=True, cutout=True)#cutout provided, use latest pointing data updated every 3 hours
    5. sdo_master_pointing(aiamap, mpo_1Week=True, cutout=True)#cutout provided, use pointing data updated weekly
    6. sdo_master_pointing(aiamap, cutout=True) #cutout provided, 3 hours updated pointing data will be used
    """
    from datetime import datetime, time, date
    print('--------------------MASTER POINTING LIST AT WORK--------------------------')
    print('Accessing Master Pointing Lists')
    jsoc_query_result=[] #contains jsoc master pointing data (returned)
    #--------------------------------------------------------------------------------------------------------------------    
    #Step 1: Ensuring the AIAMap(s) input is a list, when more than one AIAMap is passed it will be a list, but 
             #if only 1 AIAMap is passed, it will have to be converted to a list of 1. 
    aiamap_list=[]
    aiamap_copy=copy.deepcopy(aiamap)
    if isinstance(aiamap_copy,__builtins__.list)!=True: #in order to loop through one or more AIAMaps, they have to be in list form
        aiamap_list.append(aiamap_copy) #if it's not a list, make it a list (only happens when it's just 1 AIAMap)
    if isinstance(aiamap_copy,__builtins__.list)==True:
        aiamap_list=aiamap_copy 
    #--------------------------------------------------------------------------------------------------------------------    
    #Step 2: Checking Inputs
    sig=signature(sdo_master_pointing) #contains the parameters and arguments passed into the function sdo_master_pointing
    params=sig.parameters  #writes out the parameters and their associated values here    
    #print(params) #print(len(params))
    if(len(params) < 1):    #making sure aiamap and jsoc pointing data are passed
        raise ValueError("Not enough parameters passed into the function, kindly provide atleast one AIAMap and the JSOC Pointing Info file")
    if not isinstance(aiamap_list[0], (AIAMap, HMIMap)):
        raise ValueError("Input must be an AIAMap or HMIMap.")
    #-------------------------------------------------------------------------------------------------------------------            
    #Step 3: Looping through each inputed AIAMap
    n_img=np.size(aiamap_list) #List of 1 or more AIAMaps
    print('Number of AIAmaps are {0}'.format(n_img))
    for f in np.arange(0,n_img): #-------------------For each AIAMap-----------------------
        print('For AIAMap # {0}'.format(f+1))    
                                
        #Step 4: Using AIA Map Information used to construct custom JSOC Keyword names & building the query
        wavelnth=aiamap_list[f].wavelength #astropy.units.quantity.Quantity, value+units
        wavelnth_num=wavelnth.value #Just the numerical portion of the wavelength
        x0_pointing='A_'+str(int(wavelnth_num))+'_X0'
        y0_pointing='A_'+str(int(wavelnth_num))+'_Y0'
        instrot_pointing='A_'+str(int(wavelnth_num))+'_INSTROT'
        imscale_pointing='A_'+str(int(wavelnth_num))+'_IMSCALE'
    
    
        query_key='T_START , '+x0_pointing+' , '+y0_pointing+' , '+instrot_pointing+' , '+imscale_pointing+' , '+'*recnum*'
        t_stamp=aiamap_list[f].date #this is the timestamp of the fits file
        print('The query key sent to JSOC is : '+query_key)
        print('The time of this particular AIAMap is %s'%(t_stamp.time()))
        
        #Step 5: Constructing the date time string for 3hour master pointing 
        date=aiamap_list[f].date
        yr=date.year
        day=date.strftime('%d')
        month=date.strftime('%m')
        date_string='{0}-{1}-{2}/1d'.format(yr,month,day)
        #--------------------------------------------------------------------------------------------------------------------
        #Step 5: Downloading the data from JSOC if mpo_3Hours is set to true:
        if mpo_3Hours==False and mpo_1Week==False: #If user doesn't specify preference, set mpo_3Hours as true
            mpo_3Hours=True 
        
        if mpo_3Hours==True:
            
            #Step 6: Pinging the client
            c=drms.Client() #all series are available by calling the Client.series() method
            k=c.query('aia.master_pointing3h['+date_string+']',key=query_key) 
            print(k)
            #---------------------------------------------------------------------------------------------------------------------
            #Step 7: Find the pointing data at the closest timestamp to the AIAMap
            jsoc_datetime_str=k['T_START'] #the datatime is an object, this becomes a Panda Series

            #Step 7.1: Converting the date-time, which is currently in string format to a date-time array:
            i=0
            datetime_dict=[] #empty dictionary which will contain the datetime arrays from JSOC 

            while(i<jsoc_datetime_str.size): #creating date-time array from JSOC data
                first_str=jsoc_datetime_str[i] 
                date_str=first_str[0:10] #just the date
                time_str=first_str[11:19] #just the time
                year_str=int(first_str[0:4]) #forcing it since I am unable to use strptime due to leap sec incapabilities
                month_str=int(first_str[5:7])
                day_str=int(first_str[8:10])
                hour_str=int(first_str[11:13])
                min_str=int(first_str[14:16])
                sec_str=int(first_str[17:19])
    
                if (sec_str > 59): #leap second limitation
                    sec_str=59 
                
                datetime_arr=datetime(year_str,month_str,day_str,hour_str,min_str,sec_str) 
                datetime_dict.append(datetime_arr)
                i+=1

            datetime_pd=pd.Series(datetime_dict) #constructing Panda series

            #Step 7.2: finding the difference between the datetime of the AIAMap with the datetime from JSOC
            delta_list=[] #list for appending
            for date in datetime_pd:
                delta=date-t_stamp
                delta1=(abs(delta.total_seconds())) #absolute value of the time difference in seconds
                delta_list.append(delta1)
    
            timediff=pd.Series(delta_list)

            #Step 7.3: Which timedifference is minimum?
            index_time_min=timediff.index[timediff==min(timediff)].tolist() #at this index in all the panda Series & DataFrames, the time is closest to the AIAMap time
            b=k.iloc[index_time_min]
            jsoc_query_result.append(b)
        #--------------------------------------------------------------------------------------------------------------------
        #Step 8: Downloading the data from JSOC if mpo_1Week is set to true:
        if mpo_1Week==True:
            c=drms.Client() #all series are available by calling the Client.series() method
            #print(c.series(r'aia')) #this gives all series that begin with 'aia'
            k=c.query('sdo.master_pointing',key=query_key) 
            #---------------------------------------------------------------------------------------------------------------------
            #Find the pointing data at the closest timestamp to the AIAMap
            jsoc_datetime_str=k['T_START'] #the datatime is an object, this becomes a Panda Series
            
            #Converting the date-time, which is currently in string format to a date-time array:
            i=0
            datetime_dict=[] #empty dictionary which will contain the datetime arrays from JSOC 

            while(i<jsoc_datetime_str.size): #creating date-time array from JSOC data
                first_str=jsoc_datetime_str[i] 
                date_str=first_str[0:10] #just the date
                time_str=first_str[11:19] #just the time
                year_str=int(first_str[0:4]) #forcing it since I am unable to use strptime due lack of leap sec capability
                month_str=int(first_str[5:7])
                day_str=int(first_str[8:10])
                hour_str=int(first_str[11:13])
                min_str=int(first_str[14:16])
                sec_str=int(first_str[17:19])
    
                if (sec_str > 59):
                    sec_str=59 #forcing the leap second (60) to work
        
                datetime_arr=datetime(year_str,month_str,day_str,hour_str,min_str,sec_str) 
                datetime_dict.append(datetime_arr)
                i+=1

            datetime_pd=pd.Series(datetime_dict)

            #finding the difference between the datetime of the AIAMap with the datetime from JSOC
            delta_list=[] #list for appending
            for date in datetime_pd:
                delta=date-t_stamp
                delta1=(abs(delta.total_seconds())) #absolute value of the time difference in seconds
                delta_list.append(delta1)
    
            timediff=pd.Series(delta_list)
            #Which timedifference is minimum?
            index_time_min=timediff.index[timediff==min(timediff)].tolist() #at this index in all the panda Series & DataFrames, the time is closest to the AIAMap time
            b=k.iloc[index_time_min]
            jsoc_query_result.append(b)
    print('--------------------MASTER POINTING LIST ACCESS WORK COMPLETED--------------------------')
    return jsoc_query_result #returns list of jsoc pointing data for each AIAMap
    #----------------------------------------------------------------------------------------------------------

def drot_map(aiamap,reference_map):
    """
    - This function performs complete image registration, given inputs of (1) the AIAMap and (2) a reference Map
    - Additional functionalities, such as those provided in drot_map.pro need to be discussed
    """
    pass