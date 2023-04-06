"""This module only concerns the reading part of the data.
"""
import os
from xml.dom.minidom import parse
import xml.dom.minidom
from datetime import datetime
import numpy as np
# return a list of path(s) from the top directory

def get_gf3_path(top_dir):

    file_path = []       
    meta_paths = []      
    res_paths = []
    raw_paths = []
    cor_slaves_paths = []
    for top_dir, dirs, files in os.walk(top_dir):          #顶层目录，文件夹，文件        
        for dir in dirs:                                   #遍历目录下的文件夹
            file_path.append(os.path.join(top_dir, dir))   #获得子目录路径
    # print('file_path', ':', file_path)                     #查看一下子目录
    process_paths = file_path[0]           
    slc_paths = file_path[1]               # 得到 process_paths，slc_paths
    process_sub_paths = file_path[2]       # process_paths 文件夹下的 /S01B01
    # print('slc_path', ':', slc_paths)
    # print('process_sub_paths', ':', process_sub_paths)

    for file in os.listdir(slc_paths):    # 遍历 SLC 文件路径
        if not os.path.isfile(os.path.join(slc_paths, file)):
            slaves_paths = os.path.join(slc_paths, file)       # 得到 SLC 下面两个文件夹的路径
            for file in os.listdir(slaves_paths):
                if file.endswith('meta.xml'):                  # 查找后缀名为 meta.xml 的文件
                    meta_paths.append(os.path.join(slaves_paths, file))     # 添加到 meta_paths 里面
    # print('meta_paths', ':', meta_paths)

    for file in os.listdir(process_sub_paths):                #遍历 /process/S01B01 文件夹
        if not os.path.isfile(os.path.join(process_sub_paths, file)):
            cor_slaves_paths = os.path.join(process_sub_paths, file)
            print('cor_slaves_paths', ':', cor_slaves_paths)
            # for file in os.listdir(cor_slaves_paths):
            #     if file.endswith('slave.res'):                 # 查找后缀名为 slave.res 的文件
            #         res_paths.append(os.path.join(cor_slaves_paths, file))
            #     else:
            #         continue
            for file in os.listdir(cor_slaves_paths):
                if file.endswith('crop.raw'):                  # 查找后缀名为 crop.raw 的文件
                    raw_paths.append(os.path.join(cor_slaves_paths, file))
                else:
                    continue
        else:
            if file.endswith('master.res'):                            
                res_paths.append(os.path.join(process_sub_paths, file))           
            elif file.endswith('master_crop.raw'):
                raw_paths.append(os.path.join(process_sub_paths, file))
            else:
                continue

    # print('res_paths', ':', res_paths)
    # print('raw_paths', ':', raw_paths)
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
    # wavelength = speedOfLight / (float(RCF) * 1e6)

    # waveinfo
    waveParams = sensor[0].getElementsByTagName("waveParams")
    wave = waveParams[0].getElementsByTagName("wave")
    waveCode = wave[0].getElementsByTagName("waveCode")[0].childNodes[0].nodeValue
    CLA = wave[0].getElementsByTagName("centerLookAngle")[0].childNodes[0].nodeValue
    prf = wave[0].getElementsByTagName("prf")[0].childNodes[0].nodeValue
    proBandwidth = wave[0].getElementsByTagName("proBandwidth")[0].childNodes[0].nodeValue
    sampleRate = wave[0].getElementsByTagName("sampleRate")[0].childNodes[0].nodeValue
    sampleDelay = wave[0].getElementsByTagName("sampleDelay")[0].childNodes[0].nodeValue
    bandWidth = wave[0].getElementsByTagName("bandWidth")[0].childNodes[0].nodeValue
    pulseWidth = wave[0].getElementsByTagName("pulseWidth")[0].childNodes[0].nodeValue
    frameLength = wave[0].getElementsByTagName("frameLength")[0].childNodes[0].nodeValue
    groundVelocity = wave[0].getElementsByTagName("groundVelocity")[0].childNodes[0].nodeValue
    averageAltitude = wave[0].getElementsByTagName("averageAltitude")[0].childNodes[0].nodeValue
    
    # platform = root.getElementsByTagName("platform")
    # acqDate = platform[0].getElementsByTagName("CenterTime")[0].childNodes[0].nodeValue

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
    corner = imageinfo[0].getElementsByTagName("corner")
    topLeft = corner[0].getElementsByTagName("topLeft")
    TL_lat = topLeft[0].getElementsByTagName("latitude")[0].childNodes[0].nodeValue
    TL_lon = topLeft[0].getElementsByTagName("longitude")[0].childNodes[0].nodeValue
    topRight = corner[0].getElementsByTagName("topRight")
    TR_lat = topRight[0].getElementsByTagName("latitude")[0].childNodes[0].nodeValue
    TR_lon = topRight[0].getElementsByTagName("longitude")[0].childNodes[0].nodeValue
    bottomLeft = corner[0].getElementsByTagName("bottomLeft")
    BL_lat = bottomLeft[0].getElementsByTagName("latitude")[0].childNodes[0].nodeValue
    BL_lon = bottomLeft[0].getElementsByTagName("longitude")[0].childNodes[0].nodeValue
    bottomRight = corner[0].getElementsByTagName("bottomRight")
    BR_lat = bottomRight[0].getElementsByTagName("latitude")[0].childNodes[0].nodeValue
    BR_lon = bottomRight[0].getElementsByTagName("longitude")[0].childNodes[0].nodeValue
    width = imageinfo[0].getElementsByTagName("width")[0].childNodes[0].nodeValue
    height = imageinfo[0].getElementsByTagName("height")[0].childNodes[0].nodeValue
    widthspace = imageinfo[0].getElementsByTagName("widthspace")[0].childNodes[0].nodeValue
    heightspace = imageinfo[0].getElementsByTagName("heightspace")[0].childNodes[0].nodeValue
    imagebit = imageinfo[0].getElementsByTagName("imagebit")[0].childNodes[0].nodeValue

    # processinfo
    processinfo = root.getElementsByTagName("processinfo")
    incidenceAngleNearRange = processinfo[0].getElementsByTagName("incidenceAngleNearRange")[0].childNodes[0].nodeValue
    incidenceAngleFarRange = processinfo[0].getElementsByTagName("incidenceAngleFarRange")[0].childNodes[0].nodeValue
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

def read_gf3_slc(slc_path):
    with open(slc_path, "rb") as f:        
        slc = np.fromfile(f, dtype = np.float16)
        slc = complex_short2complex(slc)
        shape = (26665, 28238)
        slc = np.reshape(slc, shape)
        print(slc.shape, slc.dtype)
    return slc

def complex_short2complex(data):
        data = np.copy(data).view(np.float16).astype('float32').view(np.complex64)
        return data

if __name__ == "__main__":
    gf3_meta = []
    gf3_slc = []
    meta_paths, slc_paths = get_gf3_path("/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii")
    print("meta_path: ", meta_paths)
    print("slc_path: ", slc_paths)
    for meta_path in (meta_paths):
        print(meta_path)
        gf3_meta.append(read_metadata(meta_path))
        print("gf3_meta: ", gf3_meta)
    for slc_path in (slc_paths):
        print("slc_path: ", slc_path)
        gf3_slc.append(read_gf3_slc(slc_path))
        print("gf3_slc: ", gf3_slc)
