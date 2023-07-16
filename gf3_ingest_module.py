"""This module only concerns the reading part of the data.
"""
import os
import h5py
import re
import xml.dom.minidom
import xml.etree.ElementTree as ET
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from getAOI import orbitFit
import json 

# return a list of path(s) from the top directory
def get_gf3_path(top_dir):
    file_path = []       
    meta_paths = []      
    res_paths = []
    raw_paths = []
    cor_slaves_paths = []
    for top_dir, dirs, files in os.walk(top_dir):                 
        for dir in dirs:                                   
            file_path.append(os.path.join(top_dir, dir))   
    # print(file_path)
    process_paths = file_path[0]           
    # slc_paths = file_path[1]               # 得到 process_paths，slc_paths
    process_sub_paths = file_path[2]       # 得到 process/S01B01 文件夹的路径

    # for file in os.listdir(slc_paths):    
    #     if not os.path.isfile(os.path.join(slc_paths, file)):
    #         slaves_paths = os.path.join(slc_paths, file)       
    #         for file in os.listdir(slaves_paths):
    #             if file.endswith('meta.xml'):                  
    #                 meta_paths.append(os.path.join(slaves_paths, file))     


    for file in os.listdir(process_sub_paths):                
        if not os.path.isfile(os.path.join(process_sub_paths, file)):
            cor_slaves_paths = os.path.join(process_sub_paths, file)       
            print(cor_slaves_paths)
            for file in os.listdir(cor_slaves_paths):
                if file.endswith(('crop.raw', 'rsmp.raw')):                  
                    raw_paths.append(os.path.join(cor_slaves_paths, file))
                elif file.endswith('meta.xml'):
                    meta_paths.append(os.path.join(cor_slaves_paths, file))
                else:
                    continue

    return meta_paths, raw_paths

# return the meta data (in some sort of data structure, i.e., dictionary)
def read_metadata(matadata_path):
    DOMTree = xml.dom.minidom.parse(matadata_path)
    root = DOMTree.documentElement
    root2 = ET.parse(matadata_path).getroot()
    orbit = root.getElementsByTagName("orbitID")[0].childNodes[0].nodeValue
    direction = root.getElementsByTagName("Direction")[0].childNodes[0].nodeValue
    satellite = root.getElementsByTagName("satellite")[0].childNodes[0].nodeValue

    # sensorinfo
    sensor = root.getElementsByTagName("sensor")
    imagingMode = sensor[0].getElementsByTagName("imagingMode")[0].childNodes[0].nodeValue
    lamda = sensor[0].getElementsByTagName("lamda")[0].childNodes[0].nodeValue   
    RCF = sensor[0].getElementsByTagName("RadarCenterFrequency")[0].childNodes[0].nodeValue

    # waveinfo
    waveParams = sensor[0].getElementsByTagName("waveParams")
    wave = waveParams[0].getElementsByTagName("wave")
    waveCode = wave[0].getElementsByTagName("waveCode")[0].childNodes[0].nodeValue
    CLA = wave[0].getElementsByTagName("centerLookAngle")[0].childNodes[0].nodeValue
    prf = wave[0].getElementsByTagName("prf")[0].childNodes[0].nodeValue
    proBandwidth = wave[0].getElementsByTagName("proBandwidth")[0].childNodes[0].nodeValue
    sampleRate = wave[0].getElementsByTagName("sampleRate")[0].childNodes[0].nodeValue

    # productinfo
    productinfo = root.getElementsByTagName("productinfo")
    NoRes = productinfo[0].getElementsByTagName("NominalResolution")[0].childNodes[0].nodeValue
    productType = productinfo[0].getElementsByTagName("productType")[0].childNodes[0].nodeValue
    productFormat = productinfo[0].getElementsByTagName("productFormat")[0].childNodes[0].nodeValue

    # imageinfo
    imageinfo = root.getElementsByTagName("imageinfo")
    imagingTime = imageinfo[0].getElementsByTagName("imagingTime")
    startTime = imagingTime[0].getElementsByTagName("start")[0].childNodes[0].nodeValue
    endTime = imagingTime[0].getElementsByTagName("end")[0].childNodes[0].nodeValue
    acqDate = datetime.strptime(startTime, "%Y-%m-%d %H:%M:%S.%f")
    nearRange = imageinfo[0].getElementsByTagName("nearRange")[0].childNodes[0].nodeValue
    refRange = imageinfo[0].getElementsByTagName("refRange")[0].childNodes[0].nodeValue
    Fs = imageinfo[0].getElementsByTagName("eqvFs")[0].childNodes[0].nodeValue
    PRF = imageinfo[0].getElementsByTagName("eqvPRF")[0].childNodes[0].nodeValue
    center = imageinfo[0].getElementsByTagName("center")
    cen_lat = center[0].getElementsByTagName("latitude")[0].childNodes[0].nodeValue
    cen_lon = center[0].getElementsByTagName("longitude")[0].childNodes[0].nodeValue
    width = imageinfo[0].getElementsByTagName("width")[0].childNodes[0].nodeValue
    height = imageinfo[0].getElementsByTagName("height")[0].childNodes[0].nodeValue
    widthspace = imageinfo[0].getElementsByTagName("widthspace")[0].childNodes[0].nodeValue
    heightspace = imageinfo[0].getElementsByTagName("heightspace")[0].childNodes[0].nodeValue
    imagebit = imageinfo[0].getElementsByTagName("imagebit")[0].childNodes[0].nodeValue
    # print("PRF", PRF)
    # orbitinfo
    orbitinfo = root.getElementsByTagName("GPS")
    x_position = []
    y_position = []
    z_position = []
    x_velocity = []
    y_velocity = []
    z_velocity = []
    time = []
    # for gps_param in root2.findall('.//GPSParam'):
    #     timestamp = gps_param.find("TimeStamp").text
    #     time_mjd = convert_to_mjd(timestamp)
    #     x_pos = gps_param.find("xPosition").text
    #     y_pos = gps_param.find("yPosition").text
    #     z_pos = gps_param.find("zPosition").text
    #     x_vel = gps_param.find("xVelocity").text
    #     y_vel = gps_param.find("yVelocity").text
    #     z_vel = gps_param.find("zVelocity").text
    #     x_position.append(float(x_pos))
    #     y_position.append(float(y_pos))
    #     z_position.append(float(z_pos))
    #     x_velocity.append(float(x_vel))
    #     y_velocity.append(float(y_vel))
    #     z_velocity.append(float(z_vel))

    GPSParams = orbitinfo[0].getElementsByTagName("GPSParam")
    for gps_param in GPSParams:
        timestamp = gps_param.getElementsByTagName('TimeStamp')[0].firstChild.data
        time_mjd = convert_to_mjd(timestamp)
        # get position information
        x_pos = gps_param.getElementsByTagName('xPosition')[0].firstChild.data
        y_pos = gps_param.getElementsByTagName('yPosition')[0].firstChild.data
        z_pos = gps_param.getElementsByTagName('zPosition')[0].firstChild.data
        # get velocity information
        x_vel = gps_param.getElementsByTagName('xVelocity')[0].firstChild.data
        y_vel = gps_param.getElementsByTagName('yVelocity')[0].firstChild.data
        z_vel = gps_param.getElementsByTagName('zVelocity')[0].firstChild.data

        x_position.append(float(x_pos))
        y_position.append(float(y_pos))
        z_position.append(float(z_pos))
        x_velocity.append(float(x_vel))
        y_velocity.append(float(y_vel))
        z_velocity.append(float(z_vel)) 
        time.append(float(time_mjd))
    # print("xposition: ", x_position)       
        # orbit dictionary
    orbit = {}
    orbit = {
            'time': time,
            'x_pos': x_position,
            'y_pos': y_position,
            'z_pos': z_position,
            'x_vel': x_velocity,
            'y_vel': y_velocity,
            'z_vel': z_velocity
        }
    # print("orbit: ", orbit)
    # import pdb; pdb.set_trace()
    orbit_Fit = orbitFit(orbit)
    orbitFit_str = str(orbit_Fit)
    # print("orbitFit_str: ", orbitFit_str)
    # processinfo
    processinfo = root.getElementsByTagName("processinfo")
    DEM = processinfo[0].getElementsByTagName("DEM")[0].childNodes[0].nodeValue
    
    metadata = {
        "orbit": direction,
        "acqDate": startTime,
        "wavelength": lamda,
        "imagingMode": imagingMode,
        "centerLookAngle": CLA,
        "centerLon": cen_lon,
        "centerLat": cen_lat,
        "centerH": DEM,
        "nAzimuth": height,
        "nRange": width,
        "sampleRate": sampleRate,
        "rangeSpacing": widthspace,
        "azimuthSpacing": heightspace,
        "azimuthResolution": Fs,
        "rangeResolution": PRF,
        "orbitFit": orbitFit_str
    }
    return metadata

# return the slc data (in some sort of data structure, i.e., array)
def read_gf3_slc(slc_path):
    with open(slc_path, "rb") as f:        
        slc = np.fromfile(f, dtype=np.dtype([('real', '>i2'), ('imag', '>i2')]), count=-1)
        slc = slc['real'] + 1j * slc['imag']
        print(slc.dtype)
        slc = np.asarray(slc, dtype=np.complex64)
        print('slc:', slc, slc.dtype)
        shape = (26665, 28238)
        slc = np.reshape(slc, shape)
        slc = slc[6000:6500, 6000:6500]
        print(slc.shape, slc.dtype)
    return slc

def calculate_IFG(SLC_data, master_name):
    cal_IFG_data = {}
    master_data = SLC_data[master_name]
    for dset_name, dset_data in SLC_data.items():
        cal_IFG_data[dset_name] = np.multiply(master_data, np.conj(dset_data))
    print("IFG :", cal_IFG_data)
    return cal_IFG_data

def saveHDF(data, HDFpath, group_name):
    with h5py.File(HDFpath, "a") as f:
        if group_name not in f.keys():
            g = f.create_group(group_name)
        else:
            g = f[group_name]
        for dset_name, dset_data in data.items():
            if dset_name not in g.keys():
                g.create_dataset(dset_name, data = dset_data)
            else:
                g[dset_name][:] = dset_data
        f.flush()
        f.close()


# function for converting date to mjd
def convert_to_mjd(date_str):
    date_format = "%Y-%m-%d %H:%M:%S.%f"
    reference_date = datetime(1858, 11, 17, 12, 0, 0)
    given_date = datetime.strptime(date_str, date_format)
    # calculate julian date
    julian_date = given_date - reference_date
    julian_day = julian_date.days + julian_date.seconds / (24 * 3600)
    # convert to modified julian date
    mjd = julian_day + 2400000.5
    return mjd

# function for extracting date from string
def get_date(path, pattern):
    match = re.search(pattern, path)
    if match:      
       date = match.group(1)
       print("date: ", date)
    else:
       print("Not Found!")
    return date

def plot_slc_data(slc_data, date):
    plt.figure(dpi=200)
    plt.imshow(np.abs(slc_data), cmap='gray')
    # print("amplitude:", np.abs(slc_data))
    plt.savefig('/data/tests/hongcan/GF3/inner_mongolia/' + date +'_amplitude.png')

# function to convert complex_short data to complex data
def complex_short2complex(data):
        data = np.copy(data).view(np.float16).astype('float32').view(np.complex64)
        return data

def metadata2HDF(metadata, HDF5_path, date):
    with h5py.File(HDF5_path, 'a') as f:
        if date + "_meta" not in f.keys():
            g = f.create_group(date, metadata, + "_meta")
        else:
            g = f[date + "_meta"]
        for key, value in metadata[-1].items():
            print(key, type(value))
            if key == 'acqDate':               
                value = str(value)
            g.attrs[key] = value    
    f.close()

def slc2HDF(slc_data, HDF5_path, date):
    with h5py.File(HDF5_path, 'a') as f:
        if "SLC" not in f.keys():
            g = f.create_group("SLC")
        else:
            g = f["SLC"]
        if date not in g.keys():
            g.create_dataset(date, data=slc_data)
        else:
            g[date][...] = slc_data        
    f.close()

def readHDF(hdf5_file_path, group_name):
    dset = {}
    with h5py.File(hdf5_file_path, "r") as f:
        for key in f[group_name].keys():
            dset[key] = f[group_name][key][:]
    return dset

def metadata2JSON(metadata):
    with open('/data/tests/hongcan/GF3/inner_mongolia/metadata.json', 'w') as f:
        json.dump(metadata, f, indent=4)


# if __name__ == "__main__":
#     gf3_meta = []
#     gf3_slc = []
#     orbit = []
#     date = None
#     path = "/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii"

#     # f = h5py.File('/data/tests/hongcan/GF3/datong' + '/' + "gf3.hdf5", 'w')
#     # f.close()
#     HDF5_path = "/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii/gf3.hdf5"

#     meta_paths, slc_paths = get_gf3_path(path)
#     print("meta_path: ", meta_paths)
#     print("slc_path: ", slc_paths)

#     for meta_path in (meta_paths):
#         pattern1 = r"_(\d{8})_"
#         date = get_date(meta_path, pattern1)
#         gf3_meta.append(read_metadata(meta_path))  
#         # orbit.append(read_metadata(meta_path)[1])      
#         # metadata2HDF(gf3_meta, HDF5_path, date)
#         # metadata2HDF(orbit, HDF5_path, date)
#         # print(gf3_meta)
#         # print("orbit: ", orbit)
#         metadata2JSON(gf3_meta)
#     for slc_path in (slc_paths):
#         pattern2 = r"/(\d{8})/"
#         date = get_date(slc_path, pattern2)
#         gf3_slc.append(read_gf3_slc(slc_path))
#         slc2HDF(gf3_slc[-1], HDF5_path, date)
#         plot_slc_data(gf3_slc[-1], date)
