#!/usr/bin/env python

from matplotlib import use
use('agg')
import pyart
import sys, os, urllib, urllib2
from matplotlib import pyplot as plt
from matplotlib import colors, animation
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic, pyproj
import netCDF4
import copy
from scipy import ndimage
from JSAnimation.IPython_display import display_animation
import shutil

def fetch_gini(dir_string = 'http://motherlode.ucar.edu:8080/thredds/dodsC/satellite/IR/EAST-CONUS_4km/current/',
               pattern_match = 'EAST-CONUS_4km_IR_', hind_time_step = 0):
    #dir_string = 'http://motherlode.ucar.edu:8080/thredds/dodsC/satellite/3.9/EAST-CONUS_4km/current/'
    #pattern_match = 'EAST-CONUS_4km_3.9_'
    dirlisting = urllib2.urlopen(dir_string)
    all_lines = dirlisting.read().splitlines()  #EAST-CONUS_4km_3.9_20140123_2245.gini
    has_gini = []
    for line in all_lines:
        if '.gini' in line:
            #print line
            happy_place = line.find(pattern_match)
            end_place = line.find('.gini')
            my_upper_name = line[happy_place:end_place]+'.gini'
            url = dir_string + my_upper_name
            has_gini.append(url)
    return has_gini[hind_time_step]

def gini_grid(open_dap_dataset):
    xg, yg = np.meshgrid(open_dap_dataset.variables['x'][:]*1000.0, open_dap_dataset.variables['y'][:]*1000.0)
    vals = np.uint8(open_dap_dataset.variables['IR'][:])
    pnyc = pyproj.Proj(proj = 'lcc', 
                       lat_1 = open_dap_dataset.variables['LambertConformal'].latitude_of_projection_origin,
                       lat_2 = open_dap_dataset.variables['LambertConformal'].latitude_of_projection_origin,
                       lat_0 = open_dap_dataset.variables['LambertConformal'].latitude_of_projection_origin,
                       lon_0 = open_dap_dataset.variables['LambertConformal'].longitude_of_central_meridian )
    lon, lat = pnyc(xg, yg, inverse = True)
    mng = pyart.testing.make_empty_grid([1, xg.shape[0], xg.shape[1]], 
                                        ( (0,0), (yg.min(),yg.max()),(xg.min(), xg.max()) ))
    mng.axes['lat']['data'] = open_dap_dataset.variables['LambertConformal'].latitude_of_projection_origin
    mng.axes['lon']['data'] = open_dap_dataset.variables['LambertConformal'].longitude_of_central_meridian
    mng.axes['alt']['data'] = 0.0
    ir_fld = { 'data' : vals,
               'units' : 'Counts',
               'standard_name' : 'Sensor Counts', #non CF
               'long_name' : 'Number_of_counts_in_channel',
               'valid_max' : 256,
               'valid_min' : 0,
               '_FillValue' : 9999}
    lat_fld = {'data' : lat,
               'units' : 'degree_north',
               'standard_name' : 'latitude',
               'long_name' : 'latitude_of_grid_point',
               'valid_max' : 90,
               'valid_min' : -90}
    
    lon_fld = {'data' : lon,
               'units' : 'degree_east',
               'standard_name' : 'longitude',
               'long_name' : 'longitude_of_grid_point',
               'valid_max' : 180,
               'valid_min' : -180}
    mng.fields.update({'IR' : ir_fld, 'lat': lat_fld, 'lon': lon_fld})
    return mng

def plot_gini_grid(gini_grid, box = [-110, -70, 20, 52], resolution = 'l',
                   parallels = np.linspace(10,50, 9),
                   meridians = np.linspace(-110, -80,7),
                   vmin = None, vmax = None,
                   fld = 'IR', title = None):
    m = Basemap(llcrnrlon = box[0] ,llcrnrlat = box[2] , urcrnrlon = box[1],
                   urcrnrlat = box[3] , projection = 'mill', area_thresh =1000 ,
                   resolution='l')
    
    x, y = m(gini_grid.fields['lon']['data'], gini_grid.fields['lat']['data'])
    # create figure.
    m.drawparallels(np.linspace(10,50, 9) ,labels=[1,1,0,0])
    m.drawmeridians(np.linspace(-110, -80,7),labels=[0,0,0,1]) 
    pc = m.pcolormesh(x, y , gini_grid.fields[fld]['data'][0,:], cmap=plt.get_cmap('gray'),
                      vmin = vmin, vmax = vmax)
    m.drawcoastlines(linewidth=1.25)
    m.drawstates()
    plt.title(title)
    plt.colorbar(mappable=pc)

def fetch_heal_plot(n, level = 100, box = [-90, -85, 40, 45]):
    open_dap_11microns = fetch_gini(hind_time_step = n)
    data_11=netCDF4.Dataset(open_dap_11microns)
    grid_11 = gini_grid(data_11)
    new_ir = copy.deepcopy(grid_11.fields['IR'])
    img = new_ir['data'][0,:]
    img[np.where(img > img.max()-level)] = ndimage.median_filter(img, 2)[np.where(img > img.max()-level)]
    new_ir['data'][0, :] = img
    grid_11.fields.update({'IR_filt' : new_ir})
    f = plt.figure(figsize = [15,11])
    plot_gini_grid(grid_11, box = box, vmin = 105, vmax = 210,
                   resolution = 'h', fld = 'IR_filt', 
                   title ='11 Micron channel counts ' + data_11.time_coverage_start) 

if __name__=='__main__':
    nf = int(sys.argv[1])
    odir = sys.argv[2]
    b = 0
    for i in range(nf)[::-1]:
        filename = odir+ '/IR_image_%(d)02d.png' %{'d':b}
        b+=1
        print filename
        fetch_heal_plot(i, level = 200, box = [-100, -85, 35, 45])
        plt.savefig(filename)
    pattern='IR_image_'
    files=os.listdir(odir)
    good_files=[]
    for fl in files:
        if pattern in fl:
            good_files.append(fl)
    good_files.sort()
    i=0
    for fl in good_files:
        nn="/tbmedia_%03d.png" %i
        print odir+fl, odir+nn
        shutil.copyfile(odir+fl, odir+nn)
        i+=1
    print sys.argv
    ffstr = "ffmpeg -y -r 5 -i {0}/tbmedia_\%03d.png  {0}/{1}".format(\
            odir,sys.argv[3])
    os.system(ffstr)