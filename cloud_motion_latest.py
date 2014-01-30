#!/Users/scollis/anaconda/bin/python
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
import cv, cv2

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
    mng.axes['x_disp']['data'] = open_dap_dataset.variables['x'][:]*1000.0
    mng.axes['y_disp']['data'] = open_dap_dataset.variables['y'][:]*1000.0
    mng.axes['z_disp']['data'] = np.array([0.0])
    
    ir_fld = { 'data' : vals,
               'units' : 'Counts',
               'standard_name' : 'Sensor Counts', #non CF
               'long_name' : 'Number_of_counts_in_channel',
               'valid_max' : 256,
               'valid_min' : 0,
               '_FillValue' : 9999}
    lat_fld = {'data' : np.expand_dims(lat,0),
               'units' : 'degree_north',
               'standard_name' : 'latitude',
               'long_name' : 'latitude_of_grid_point',
               'valid_max' : 90,
               'valid_min' : -90}
    
    lon_fld = {'data' : np.expand_dims(lon,0),
               'units' : 'degree_east',
               'standard_name' : 'longitude',
               'long_name' : 'longitude_of_grid_point',
               'valid_max' : 180,
               'valid_min' : -180}
    x_fld = {'data' : np.expand_dims(xg,0),
             'units' : 'meters',
             'standard_name' : 'x',
             'long_name' : 'displacement along x axis'}
    y_fld = {'data' : np.expand_dims(yg,0),
             'units' : 'meters',
             'standard_name' : 'y',
             'long_name' : 'displacement along y axis'}
    mng.fields.update({'IR' : ir_fld, 'latg': lat_fld, 'long': lon_fld, 'xg':x_fld, 'yg':y_fld})
    mng.axes['time']['units'] = 'seconds since '+open_dap_dataset.time_coverage_start
    return mng

def plot_gini_grid(gini_grid, box = [-110, -70, 20, 52], resolution = 'l',
                   parallels = None,
                   meridians = None,
                   vmin = None, vmax = None,
                   fld = 'IR', title = None):
    m = Basemap(llcrnrlon = box[0] ,llcrnrlat = box[2] , urcrnrlon = box[1],
                   urcrnrlat = box[3] , projection = 'mill', area_thresh =1000 ,
                   resolution = resolution)
    
    x, y = m(gini_grid.fields['long']['data'][0], gini_grid.fields['latg']['data'][0])
    # create figure.
    if parallels == None:
        parallels = np.linspace(10,50, 9)
    if meridians == None:
        meridians = np.linspace(-110, -80,7)
    m.drawparallels(parallels, labels=[1,1,0,0])
    m.drawmeridians(meridians,labels=[0,0,0,1]) 
    pc = m.pcolormesh(x, y , gini_grid.fields[fld]['data'][0,:], cmap=plt.get_cmap('gray'),
                      vmin = vmin, vmax = vmax)
    plt.title(title)
    m.drawcoastlines(linewidth=1.25)
    m.drawstates()
    plt.colorbar(mappable=pc, label = 'Counts')

def plot_gini_grid_vectors(gini_grid, box = [-110, -70, 20, 52], resolution = 'l',
                   parallels = None,
                   meridians = None,
                   vmin = None, vmax = None,
                   fld = 'IR', title = None,
                   degrade = 5, u='u', v='v', scale = 200):
    m = Basemap(llcrnrlon = box[0] ,llcrnrlat = box[2] , urcrnrlon = box[1],
                   urcrnrlat = box[3] , projection = 'mill', area_thresh =1000 ,
                   resolution = resolution)
    
    x, y = m(gini_grid.fields['long']['data'][0], gini_grid.fields['latg']['data'][0])
    # create figure.
    if parallels == None:
        parallels = np.linspace(10,50, 9)
    if meridians == None:
        meridians = np.linspace(-110, -80,7)
    m.drawparallels(parallels, labels=[1,1,0,0])
    m.drawmeridians(meridians,labels=[0,0,0,1]) 
    pc = m.pcolormesh(x, y , gini_grid.fields[fld]['data'][0,:], cmap=plt.get_cmap('gray'),
                      vmin = vmin, vmax = vmax)
    qq = m.quiver(x[::degrade,::degrade], y[::degrade,::degrade], 
                  gini_grid.fields[u]['data'][0,::degrade,::degrade], 
                  gini_grid.fields[v]['data'][0,::degrade,::degrade], scale=scale)

    plt.title(title)
    m.drawcoastlines(linewidth=1.25)
    m.drawstates()
    plt.colorbar(mappable=pc, label = 'Counts')

def plot_gini_grid_barbs(gini_grid, box = [-110, -70, 20, 52], resolution = 'l',
                   parallels = None,
                   meridians = None,
                   vmin = None, vmax = None,
                   fld = 'IR', title = None,
                   degrade = 5, u='u', v='v'):
    m = Basemap(llcrnrlon = box[0] ,llcrnrlat = box[2] , urcrnrlon = box[1],
                   urcrnrlat = box[3] , projection = 'mill', area_thresh =1000 ,
                   resolution = resolution)
    
    x, y = m(gini_grid.fields['long']['data'][0], gini_grid.fields['latg']['data'][0])
    # create figure.
    if parallels == None:
        parallels = np.linspace(10,50, 9)
    if meridians == None:
        meridians = np.linspace(-110, -80,7)
    m.drawparallels(parallels, labels=[1,1,0,0])
    m.drawmeridians(meridians,labels=[0,0,0,1]) 
    pc = m.pcolormesh(x, y , gini_grid.fields[fld]['data'][0,:], cmap=plt.get_cmap('gray'),
                      vmin = vmin, vmax = vmax)
    qq = m.barbs(x[::degrade,::degrade], y[::degrade,::degrade], 
                  gini_grid.fields[u]['data'][0,::degrade,::degrade], 
                  gini_grid.fields[v]['data'][0,::degrade,::degrade])

    plt.title(title)
    m.drawcoastlines(linewidth=1.25)
    m.drawstates()
    plt.colorbar(mappable=pc, label = 'Counts')


def cv2array(im):
    depth2dtype = {
    cv.IPL_DEPTH_8U: 'uint8',
    cv.IPL_DEPTH_8S: 'int8',
    cv.IPL_DEPTH_16U: 'uint16',
    cv.IPL_DEPTH_16S: 'int16',
    cv.IPL_DEPTH_32S: 'int32',
    cv.IPL_DEPTH_32F: 'float32',
    cv.IPL_DEPTH_64F: 'float64',
    }

    arrdtype=im.depth
    a = np.fromstring(
    im.tostring(),
    dtype=depth2dtype[im.depth],
    count=im.width*im.height*im.nChannels)
    a.shape = (im.height,im.width,im.nChannels)
    return a

def array2cv(a):
    dtype2depth = {
    'uint8': cv.IPL_DEPTH_8U,
    'int8': cv.IPL_DEPTH_8S,
    'uint16': cv.IPL_DEPTH_16U,
    'int16': cv.IPL_DEPTH_16S,
    'int32': cv.IPL_DEPTH_32S,
    'float32': cv.IPL_DEPTH_32F,
    'float64': cv.IPL_DEPTH_64F,
    }
    try:
        nChannels = a.shape[2]
    except:
        nChannels = 1
        cv_im = cv.CreateImageHeader((a.shape[1],a.shape[0]),
        dtype2depth[str(a.dtype)], nChannels)
        cv.SetData(cv_im, a.tostring(),a.dtype.itemsize*nChannels*a.shape[1])
    return cv_im

def get_optic_flow_fb(im0, im1, winSize = 5, n_iter = 40, levels = 1):
    
    #im0 = (im0).astype('uint8')
    #im1 = (im1).astype('uint8')
    flow = cv2.calcOpticalFlowFarneback(im0, im1, False, levels, winSize, n_iter, 7, 1.5, 0)
    return flow[:,:,0], flow[:,:,1]

def get_optic_flow(im0, im1, winSize = (5,5)):
    #LK method
    
    im0 = (im0).astype('uint8')
    im1 = (im1).astype('uint8')
    
    im0_cv = array2cv(im0)
    im1_cv = array2cv(im1)

    velx = cv.CreateImage((im0_cv.width, im0_cv.height), cv.IPL_DEPTH_32F,1)
    vely = cv.CreateImage((im0_cv.width, im0_cv.height), cv.IPL_DEPTH_32F,1)


    cv.CalcOpticalFlowLK(im0_cv, im1_cv, winSize, velx, vely)
    velx_np = cv2array(velx)
    vely_np = cv2array(vely)
    

    return velx_np, vely_np

def doflow_lk(first_frame, second_frame, winSize = (5,5), filter_len = 10, sig_min = 150 ):
    im0 = copy.deepcopy(first_frame.fields['IR_filt']['data'])
    im0[np.where(im0 < sig_min)] = sig_min
    im1 = copy.deepcopy(second_frame.fields['IR_filt']['data'])
    im1[np.where(im1 < sig_min)] = sig_min
    sim0 = (im0 - im0.min())*(im0.max()/(im0.max()-im0.min()))
    sim1 = (im1 - im1.min())*(im1.max()/(im1.max()-im1.min()))
    
    u, v = get_optic_flow(sim0[0], 
                          sim1[0],
                          winSize = winSize)
    t1 = netCDF4.num2date(second_frame.axes['time']['data'][0], units = second_frame.axes['time']['units'])
    t0 = netCDF4.num2date(first_frame.axes['time']['data'][0], units = first_frame.axes['time']['units'])
    dt = (t1-t0).seconds
    dx = np.expand_dims(np.gradient(second_frame.fields['xg']['data'][0])[1], 0)#lazy.. fix this.. 
    dy = np.expand_dims(np.gradient(second_frame.fields['yg']['data'][0])[0], 0)
    u_fld = {'data' : dt * ndimage.median_filter(u.reshape([1,u.shape[0], u.shape[1]]),filter_len)/dx,
                                        'units' :'pixels',
                                        'standard_name' : 'disp',
                                        'long name' : 'todo'}

    v_fld = {'data' : dt * ndimage.median_filter( v.reshape([1,v.shape[0], v.shape[1]]),filter_len)/dy,
                                        'units' :'pixels',
                                        'standard_name' : 'disp',
                                        'long name' : 'todo'}

    return u_fld, v_fld

def doflow_fb(first_frame, second_frame, winSize = (5,5), filter_len = 10, sig_min = 150, n_iter = 40, levels = 1):
    im0 = copy.deepcopy(first_frame.fields['IR_filt']['data'])
    im0[np.where(im0 < sig_min)] = sig_min
    im1 = copy.deepcopy(second_frame.fields['IR_filt']['data'])
    im1[np.where(im1 < sig_min)] = sig_min
    sim0 = (im0 - im0.min())*(im0.max()/(im0.max()-im0.min()))
    sim1 = (im1 - im1.min())*(im1.max()/(im1.max()-im1.min()))
    
    u, v = get_optic_flow_fb(sim0[0], 
                          sim1[0],
                          winSize = winSize[0], n_iter=n_iter, levels=levels)
    t1 = netCDF4.num2date(second_frame.axes['time']['data'][0], units = second_frame.axes['time']['units'])
    t0 = netCDF4.num2date(first_frame.axes['time']['data'][0], units = first_frame.axes['time']['units'])
    dt = (t1-t0).seconds
    dx = np.expand_dims(np.gradient(second_frame.fields['xg']['data'][0])[1], 0)
    dy = np.expand_dims(np.gradient(second_frame.fields['yg']['data'][0])[0], 0)
    u_fld = {'data' : dt * ndimage.median_filter(u.reshape([1,u.shape[0], u.shape[1]]),filter_len)/dx,
                                        'units' :'pixels',
                                        'standard_name' : 'disp',
                                        'long name' : 'todo'}

    v_fld = {'data' : dt * ndimage.median_filter( v.reshape([1,v.shape[0], v.shape[1]]),filter_len)/dy,
                                        'units' :'pixels',
                                        'standard_name' : 'disp',
                                        'long name' : 'todo'}

    return u_fld, v_fld

def cof(f,s,l):
    first_frame_dap = netCDF4.Dataset(fetch_gini(hind_time_step = f))
    first_frame_grid = gini_grid(first_frame_dap)
    second_frame_dap = netCDF4.Dataset(fetch_gini(hind_time_step = s))
    second_frame_grid = gini_grid(second_frame_dap)
    first_frame_dap.close()
    second_frame_dap.close()
    new_ir = copy.deepcopy(first_frame_grid.fields['IR'])
    img = new_ir['data'][0,:]
    level = 200
    img[np.where(img > img.max()-level)] = ndimage.median_filter(img, 2)[np.where(img > img.max()-level)]
    new_ir['data'][0, :] = img
    first_frame_grid.fields.update({'IR_filt' : new_ir})
    new_ir = copy.deepcopy(second_frame_grid.fields['IR'])
    img = new_ir['data'][0,:]
    level =200
    img[np.where(img > img.max()-level)] = ndimage.median_filter(img, 2)[np.where(img > img.max()-level)]
    new_ir['data'][0, :] = img
    second_frame_grid.fields.update({'IR_filt' : new_ir})
    umot, vmot = doflow_lk(first_frame_grid, second_frame_grid, winSize = (5,5), filter_len = 5)
    umot_fb, vmot_fb = doflow_fb(first_frame_grid, 
                                 second_frame_grid, winSize = (20,20), filter_len = 10,
                                 levels = l)
    second_frame_grid.fields.update({'u' : umot, 'v' : vmot})
    second_frame_grid.fields.update({'u_fb' : umot_fb, 'v_fb' : vmot_fb})
    return second_frame_grid

if __name__=='__main__':
    my_grid = cof(1,0,3)
    dateobj = netCDF4.num2date(my_grid.axes['time']['data'][0], units = my_grid.axes['time']['units'])
    datestr = '{:%Y%m%d%H%M%S}'.format(dateobj)
    dodir = sys.argv[1]
    iodir = sys.argv[2]
    image_file = iodir + '/' + 'flow_' + datestr + '.png'
    data_file = dodir + '/' + 'flow_' + datestr + '.nc' 
    print(data_file)
    pyart.io.write_grid(data_file, my_grid)     
    f = plt.figure(figsize = [15,11])
    plot_gini_grid_vectors(my_grid, box = [-110, -85, 30, 45], resolution = 'h', fld = 'IR_filt',
            meridians =  np.linspace(-110, -80,25), parallels = np.linspace(10,50, 33),
            title = my_grid.axes['time']['units'], degrade = 10, 
            u='u_fb', v='v_fb', scale = 100)
    plt.savefig(image_file)
    