#!/usr/bin/env python3


"""
Created on Fri Jul  1 13:31:47 2022

@author: mtcontre
"""

import os
#from osgeo import ogr
import geopandas
import numpy as np
import sys
from osgeo import gdal

script_path = '/workspace/chg/'
os.chdir(script_path)
import OM2D_functions

#---------------------------------------
# 1. Set file names and script options
#---------------------------------------
file_path = '/pontus/mcontre3/NOAA-NWM/NWMCoupling/Model_120m_62_8/'
mesh = 'fort'

# Saving files options:
    # If empty (""), doesn't save the file
    # Otherwise, provide the name
fort425_ascii = "/pontus/mcontre3/NOAA-NWM/NWMCoupling/Pp_fort_files/Sandy_storm2.425" # Fort428 and 27 for precipitation
fort425_tiff = ""#"PrecipRate" # Save tiff files for precipitation

# Download time series options
    # If init or end are empty (""), it doesn't download the data nor save fort430
init = "20121025 00:00:00"
end = "20121116 00:00:00"
#end  = "20110820 00:20:00"

# Select repository to extract pp data
pp_repository_version = 1 # 1) Historical data, 2) Operational model
dt_pp = 300# 3600 #300

# Select domain to clip the files
pp_dom = 'bbox' # 'global', 'domain', 'bbox'
#bbox = [[-82.8042,  -61.0238],   [15.6856,   44.5221]] # Provide polygon coordinates of corners of the box
# Irene
# bbox = [[-82.8042, -61.0238], [15.6856, 44.5221]]
# Sandy
bbox =[[-84.8221, -66.4547],[8.8742, 43.1192]]
# Irma
# bbox=[[-87.8094, -59.5365], [12.0610, 35.5033]]
# Harvey
# bbox = [[-98.5150, -59.9330], [9.9770, 31.4200]]

os.chdir(file_path)

#---------------------------------------
# 1. Define domain covered y precipitation
#---------------------------------------   

from shapely.geometry import box

if pp_dom=="domain":
    clip_opt = True
    p,t,tree_p,bd = OM2D_functions.readMesh(mesh)
    min_lat = min(p.geometry.x)
    min_lon = min(p.geometry.y)
    max_lat = max(p.geometry.x)
    max_lon = max(p.geometry.y)
    bbox = box(min_lat,min_lon,max_lat, max_lon)
    bbox_geo = geopandas.GeoDataFrame({'geometry': bbox}, index=[0], crs='epsg:4326')
    coords = OM2D_functions.getFeatures(bbox_geo)
elif pp_dom=="global":
    clip_opt = False
elif pp_dom=="bbox":
    clip_opt = True
    if len(bbox)==2:
        min_lat = min(bbox[:][0])
        min_lon = min(bbox[:][1])
        max_lat = max(bbox[:][0])
        max_lon = max(bbox[:][1])
        bbox = box(min_lat,min_lon,max_lat, max_lon)
        bbox_geo = geopandas.GeoDataFrame({'geometry': bbox}, index=[0], crs='epsg:4326')
        coords = OM2D_functions.getFeatures(bbox_geo)
    else:
        sys.exit('Length of bbox must be 2x2')
else:
    clip_opt = False
    sys.exit('Unrecognized option to define precipitation domain')

os.chdir(script_path)

#---------------------------------------
# 2. Download data and save files
#---------------------------------------     
    
from datetime import datetime
#import netCDF4 as nc
import requests

init_date = datetime.strptime(init, '%Y%m%d %H:%M:%S')
end_date = datetime.strptime(end, '%Y%m%d %H:%M:%S')
df = (end_date-init_date).total_seconds()
ndt = int(np.round(df/dt_pp) + 1)
    

if pp_repository_version ==1:
    url= 'https://mtarchive.geol.iastate.edu/'  #Historical data
    response = requests.get(url+init_date.strftime('%Y')+'/'+init_date.strftime('%m')+'/'+init_date.strftime('%d')+'/mrms')
    if response.status_code != 200:
        print('Historical MRMS is not available \n')
        print('Switching to NSSL (nmq folder)')
        pp_repository_version = 3
elif pp_repository_version ==2:
    url='https://mrms.ncep.noaa.gov/data/2D/PrecipRate/'   # Operational data
else:
    sys.exit('Unrecognized repository option')
    


import datetime
import wget
import gzip
import shutil
#from osgeo import gdal

if fort425_tiff:
    if not os.path.isdir('fort425_tif'):
        print('Creating directory '+ 'fort425_tif')
        os.mkdir('fort425_tif', 0o666)

t_name = [init_date + datetime.timedelta(seconds=i*dt_pp) for i in range(ndt)]
if fort425_ascii:
    f = open(fort425_ascii,'w')
    f.write('Precipitation in Oceanweather Format                 '\
               + init_date.strftime("%Y%m%d%H%M")+'   '\
               + end_date.strftime("%Y%m%d%H%M")+'\n')
        
                
if init !="" and end!="":
    
    fix_id = []
    path_all = []
    print('Checking if files exist...')
    for i in range(ndt):
        #print('Checking file '+str(i)+' of '+str(ndt)+' files \n' )
        # Checking if data exists
        yy = t_name[i].strftime('%Y')
        MM = t_name[i].strftime('%m')
        dd = t_name[i].strftime('%d')
        if pp_repository_version ==1:
            filename =  '/PrecipRate_00.00_'+t_name[i].strftime('%Y%m%d-%H%M%S')+'.grib2.gz'
            path = url + yy + '/' + MM + '/' + dd +'/mrms/reanalysis/PrecipRate/'+ filename 
            
        elif pp_repository_version==2:
            filename =  '/MRMS_PrecipRate_00.00_'+t_name[i].strftime('%Y%m%d-%H%M%S')+'.grib2.gz'
            path = url + filename 
        elif pp_repository_version==3:
            filename =  t_name[i].strftime('%Y%m%d-%H%M%S')+'.nc'
            path = url + yy + '/' + MM + '/' + dd +'/nmq/tile2/data/QPESUMS/grid/q2rad_hsr_nc/short_qpe/'+ filename 
        path_all.append(path)    
        request = requests.get(path)
        if request.status_code!=200:
            fix_id.append(i)
            
    for i in fix_id:
        try:
           if i>0:
               path_or = path_all[i]
               path_all[i] = path_all[i-1]
               print('WARNING: Using '+path_all[i]+' instead of '+path_or)
           else:
               path_or = path_all[i]
               path_all[i] = path_all[i+1]   
               print('WARNING: Using '+path_all[i]+' instead of '+path_or)
        except:
           print('ERROR: File not found. Path: '+ path)
                
   
    for i in range(ndt):
        
        print('Processing file '+str(i)+' of '+str(ndt)+' files \n' )
        # Download data
        try:
            aux = wget.download(path_all[i])
        except:
            aux = wget.download(path_all[i])
            pass
        if pp_repository_version ==1:
           
            with gzip.open(aux, 'rb') as f:
                with open(aux[0:-3],'wb') as f_out:
                    shutil.copyfileobj(f, f_out)
            os.remove(aux)
            aux = aux[0:-3]
            print('\n')
        elif pp_repository_version==2:
            with gzip.open(aux, 'rb') as f:
                with open(aux[0:-3],'wb') as f_out:
                    shutil.copyfileobj(f, f_out)
            os.remove(aux)
            aux = aux[0:-3]
            print('\n')
            
        elif pp_repository_version==3:
            print('\n')
            # Read data
            aux2 = gdal.Open(aux)
            dssubset = aux2.GetSubDatasets()[0][0]
            aux = aux[0:-3]+'.tif'
            OM2D_functions.convert_netcdf_to_geotiff(dssubset,aux, band_number=1, verbose=False)
            os.remove(aux[0:-4]+'.nc')
            
        # Read data
        dsorig = gdal.Open(aux)
        os.remove(aux)
                        
        # Reproject to EPSG:4326 and clip is user-defined
        dst_crs = 'EPSG:4326'
        if clip_opt==True:
            if fort425_tiff:
                fname = './fort425_tif/'+fort425_tiff+t_name[i].strftime("%Y%m%d%H%M")+'.tif'
                ds = gdal.Warp(fname,dsorig,dstSRS = dst_crs,outputBounds=(min_lat,min_lon,max_lat,max_lon))
            else:
                ds = gdal.Warp("dem_clipped_storm.tif",dsorig,dstSRS = dst_crs,outputBounds=(min_lat,min_lon,max_lat,max_lon))
        else:
            if fort425_tiff:
                fname = './fort425_tif/'+fort425_tiff+t_name[i].strftime("%Y%m%d%H%M")+'.tif'
                ds = gdal.Warp(fname,dsorig,dstSRS = dst_crs)
            else:
                ds = gdal.Warp("dem_clipped_storm.tif",dsorig,dstSRS = dst_crs)
        dsorig = None    
        
        ds        
                
        # Extract data to be written
        gt = ds.GetGeoTransform() 
        proj = ds.GetProjection()
        band = ds.GetRasterBand(1)
        pp = band.ReadAsArray()
        pp[pp==-999]=0
        pp = np.flipud(pp)
            
        # Write head of each time step   
        if fort425_ascii:          
            f = open(fort425_ascii,'a')
            f.write('iLat='+'{: >4}'.format(pp.shape[0])+'iLong='+'{: >4}'.format(pp.shape[1])+\
                    'DX='+'{: >.4f}'.format(gt[1])+'DY='+'{: >.4f}'.format(-1*gt[5])+\
                        'SWLat='+'{: >8}'.format(gt[3]+gt[5]*pp.shape[0])+'SWLon='+'{: >8}'.format(gt[0])+\
                     'DT='+t_name[i].strftime("%Y%m%d%H%M")+'\n')  
            pp_array = pp.flatten()    
                 
            # Write matrix for each time step 
            for j in range(int(np.ceil(len(pp_array)/8))):
                aux = pp_array[8*j:8*j+8]
                pp_str = ''
                for k in range(min(8,len(aux))):
                    aux2 = '{:.5f}'.format(aux[k])
                    pp_str = pp_str + '   '+aux2[0:7]
                f.write(pp_str+'\n')
             
        # Closing files and deleting extra files if needed    
        ds = None     
        if clip_opt==False:   
            os.remove("dem_clipped_storm.tif")
             

if fort425_ascii:                
    f.close()
    with open(fort425_ascii, 'rb') as f_in:
        with gzip.open(fort425_ascii+'.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
                     
                        
                    
    
              
                
               
    
       # ds = cfgrib.open_dataset(aux[0:-3])
       
       # ds = cfgrib.open_dataset('tif_clipped_wgs84.tiff')
       # # Read precipitation
       # ds_pp = ds.get("paramId_0")
       # df = ds.to_dataframe()
       # # Transform coordinates from 0-360 to -180 to 180
       # latitudes = df.index.get_level_values("latitude")
       # longitudes = df.index.get_level_values("longitude")
       # if max(longitudes)>180:
       #     map_function = lambda lon: (lon - 360) if (lon > 180) else lon
       #     remapped_longitudes = longitudes.map(map_function)
       #     df["longitude"] = remapped_longitudes
       #     df["latitude"] = latitudes
       
    
    
           
           
           
    
    
    
    
