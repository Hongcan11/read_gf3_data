"""This module only concerns the reading part of the data.
"""
import os
import h5py
import re
from xml.dom.minidom import parse
import xml.dom.minidom
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

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

    process_paths = file_path[0]           
    slc_paths = file_path[1]               # 得到 process_paths，slc_paths
    process_sub_paths = file_path[2]       # 得到 process/S01B01 文件夹的路径

    for file in os.listdir(slc_paths):    
        if not os.path.isfile(os.path.join(slc_paths, file)):
            slaves_paths = os.path.join(slc_paths, file)       
            for file in os.listdir(slaves_paths):
                if file.endswith('meta.xml'):                  
                    meta_paths.append(os.path.join(slaves_paths, file))     

    for file in os.listdir(process_sub_paths):                
        if not os.path.isfile(os.path.join(process_sub_paths, file)):
            cor_slaves_paths = os.path.join(process_sub_paths, file)       
            print(cor_slaves_paths)
            for file in os.listdir(cor_slaves_paths):
                if file.endswith(('crop.raw', 'rsmp.raw')):                  
                    raw_paths.append(os.path.join(cor_slaves_paths, file))
                else:
                    continue

    return meta_paths, raw_paths

# return the meta data (in some sort of data structure, i.e., dictionary)
def read_metadata(matadata_path):
    DOMTree = xml.dom.minidom.parse(matadata_path)
    root = DOMTree.documentElement

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

    # processinfo
    processinfo = root.getElementsByTagName("processinfo")
    DEM = processinfo[0].getElementsByTagName("DEM")[0].childNodes[0].nodeValue
    
    metadata = {
        "orbit": direction,
        "acqDate": acqDate,
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
    }
    return metadata

# return the slc data (in some sort of data structure, i.e., array)
def read_gf3_slc(slc_path):
    with open(slc_path, "rb") as f:        
        slc = np.fromfile(f, dtype=np.dtype([('real', '>i2'), ('imag', '>i2')]), count=-1)
        slc = slc['real'] + 1j * slc['imag']
        print('slc:', slc)
        shape = (26665, 28238)
        slc = np.reshape(slc, shape)
        print(slc.shape, slc.dtype)
    return slc

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
    print("amplitude:", np.abs(slc_data))
    plt.savefig('/data/tests/hongcan/GF3/map/' + date +'_amplitude.png')

# function to convert complex_short data to complex data
# def complex_short2complex(data):
#         data = np.copy(data).view(np.float16).astype('float32').view(np.complex64)
#         return data

def metadata2HDF(metadata, HDF5_path, date):
    with h5py.File(HDF5_path, 'a') as f:
        if "metadata" not in f.keys():
            g = f.create_group(date + "_meta")
        else:
            g = f[date + "_meta"]
        for key, value in metadata[0].items():
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

if __name__ == "__main__":
    gf3_meta = []
    gf3_slc = []
    date = None
    path = "/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii"

    f = h5py.File(path + '/' + "gf3.hdf5", 'w')
    f.close()
    HDF5_path = path + '/' + "gf3.hdf5"

    meta_paths, slc_paths = get_gf3_path(path)
    print("meta_path: ", meta_paths)
    print("slc_path: ", slc_paths)

    for meta_path in (meta_paths):
        pattern1 = r"_(\d{8})_"
        date = get_date(meta_path, pattern1)
        gf3_meta.append(read_metadata(meta_path))        
        metadata2HDF(gf3_meta, HDF5_path, date)
    
    for slc_path in (slc_paths):
        pattern2 = r"/(\d{8})/"
        date = get_date(slc_path, pattern2)
        gf3_slc.append(read_gf3_slc(slc_path))
        slc2HDF(gf3_slc[-1], HDF5_path, date)
        plot_slc_data(gf3_slc[-1], date)
