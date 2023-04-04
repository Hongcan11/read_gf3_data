from xml.dom.minidom import parse
import xml.dom.minidom
from datetime import datetime
# speedOfLight = 299792458.0


def read_metadata(matadata_path):
    DOMTree = xml.dom.minidom.parse(matadata_path)
    root = DOMTree.documentElement

    orbit = root.getElementsByTagName("orbitID")[0].childNodes[0].nodeValue
    direction = root.getElementsByTagName("Direction")[0].childNodes[0].nodeValue
    satellite = root.getElementsByTagName("satellite")[0].childNodes[0].nodeValue
    print(orbit, direction, satellite)

    # sensorinfo
    sensor = root.getElementsByTagName("sensor")
    imagingMode = sensor[0].getElementsByTagName("imagingMode")[0].childNodes[0].nodeValue
    lamda = sensor[0].getElementsByTagName("lamda")[0].childNodes[0].nodeValue   
    RCF = sensor[0].getElementsByTagName("RadarCenterFrequency")[0].childNodes[0].nodeValue
    print(sensor, imagingMode, lamda, RCF)
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
    print(waveCode, CLA, prf, proBandwidth, sampleRate, sampleDelay, bandWidth, pulseWidth, frameLength, groundVelocity, averageAltitude)

    # platform = root.getElementsByTagName("platform")
    # acqDate = platform[0].getElementsByTagName("CenterTime")[0].childNodes[0].nodeValue

    # productinfo
    productinfo = root.getElementsByTagName("productinfo")
    NoRes = productinfo[0].getElementsByTagName("NominalResolution")[0].childNodes[0].nodeValue
    productType = productinfo[0].getElementsByTagName("productType")[0].childNodes[0].nodeValue
    productFormat = productinfo[0].getElementsByTagName("productFormat")[0].childNodes[0].nodeValue
    print(NoRes, productType, productFormat)

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
    print(startTime, acqDate)
    
    
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
    print(metadata)
    return metadata
if __name__ == "__main__":
    read_metadata("/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii/slc/GF3C_KSC_FSII_004562_E105.4_N38.2_20230217_L1A_HHHV_L10000115189/GF3C_KSC_FSII_004562_E105.4_N38.2_20230217_L1A_HHHV_L10000115189.meta.xml")


# metadata = {
#         "orbit": geometry, /
#         "acqDate": acqDate, /
#         "azimuth0time": firstLineTime,  
#         "range0time": range0time,    AbsRoot.getAttributeDouble(AbstractMetadata.slant_range_to_first_pixel)/ speedOfLight
#         "PRF": PRF,
#         "RSR": RSR,  /
#         "wavelength": wavelength, /
#         "orbitFit": orbFit,
#         "rangeSpacing": rangeSampling,  height 方位 /
#         "azimuthSpacing": azimuthSampling,  /
#         "centerLon": centerLon,  /
#         "centerLat": centerLat,  /
#         "centerH": centerH,  /
#         "nAzimuth": nAzimuth, /
#         "nRange": nRange, /
#         "swath": swath,
#         "centerAzimuth": centerAzimuth,
#         "beta0": beta0,
#         "azimuthResolution": azRes, /
#         "rangeResolution": rRes, /
#         "nBursts": nBursts,
#         "burstInfo": burstInfo,
#         "steeringRate": steeringRate,
#         "azFmRateArray": azFmRateArray,
#         "dcPolyArray": dcPolyArray,
#         "PRI": pri,
#         "rank": rank,
#         "chirpRate": chirpRate,
#     }