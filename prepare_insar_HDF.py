import os
from selectPSC import selectPSc, plot_psc
import gf3_ingest_module
import json 
import numpy as np
import h5py
from datetime import datetime
import classes_2    

def prepare_insar_HDF(parms, path):
    # prepare dir:
    if not os.path.exists(parms["outDir"]):
        os.makedirs(parms["outDir"])
    # load stack:
    stackHDF = parms["outDir"] + os.sep + "stack_" + parms["stackId"] + ".hdf5"
    # if not os.path.isfile(stackHDF):
    stack, stackHDF = data2HDF(parms, path)
    print("stackHDF: ", stackHDF)
    psc = selectPSc(stackHDF, D_A_thr=parms["D_A_thr"])
    plot_psc(psc, parms["outDir"])
    # save to HDF:
    HDFfile = parms["outDir"] + os.sep + "insar_" + parms["stackId"] + ".hdf5"
    # meta_data = fromJSON(jsonpath)
    print("JSONmetadata: ", stack.metadata)
    saveHDF(psc, stack, HDFfile)
    #  open hdf5
    print("Saved.")
    data = openHDF(HDFfile)
    print("Opened!")

def data2HDF(parms, path):
    metapath = []
    slcpath = []
    metapath, slcpath = gf3_ingest_module.get_gf3_path(path)
    print("path: ", metapath, slcpath)
    # prepare output:
    outHDF = (
            parms["outDir"] + os.sep + "stack_" + parms["stackId"] + ".hdf5"
        )
    # initialise stack:
    stack = classes_2.Stack(
        parms["stackId"],
        parms["stackDir"]
    )
    stack.readData()
    print("metadata: ", stack.metadata)   
    print("Metadata read.")
    for slc_path in (slcpath):
        pattern2 = r"/(\d{8})/"
        slcdata = []
        date = gf3_ingest_module.get_date(slc_path, pattern2)
        slcdata.append(gf3_ingest_module.read_gf3_slc(slc_path))
        gf3_ingest_module.slc2HDF(slcdata[-1], outHDF, date)
    SLCdata = gf3_ingest_module.readHDF(outHDF, 'SLC')
    print("SLCdata: ", SLCdata)
    IFG = gf3_ingest_module.calculate_IFG(SLCdata, parms["master_name"])
    gf3_ingest_module.saveHDF(IFG, outHDF, "IFG")
    print("SLC data read.")     
    return stack, outHDF  

def saveHDF(psc, stack, HDFfile):
    with h5py.File(HDFfile, "w") as f:
        for k, v in psc.items():
            f.create_dataset("psc/" + k, data=v)
    print("acqDate: ", stack.acqDates)
    formatted_dates = []
    for date_string in stack.acqDates:
        datetime_obj = datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S.%f")
        formatted_date = datetime_obj.strftime("%Y%m%d")
        formatted_dates.append(formatted_date)
    print("formatted_date: ", formatted_dates)
    with h5py.File(HDFfile, "a") as f:
        if "stack/slcs" not in f:
            f.create_dataset("stack/slcs", data=np.array(formatted_dates, dtype='S8'))
        else:
            f["stack/slcs"][...] = np.array(formatted_dates, dtype='S8')
        master_acqDate = stack.acqDates[1]
        print("master_acqDate: ", master_acqDate)
        f["stack"].attrs["orbit"] = stack.metadata[master_acqDate]["orbit"]
        f["stack"].attrs["master_date"] = stack.metadata[master_acqDate]["acqDate"]
        f["stack"].attrs["wavelength"] = stack.metadata[master_acqDate]["wavelength"]
        f["stack"].attrs["rangeSpacing"] = stack.metadata[master_acqDate]["rangeSpacing"]
        f["stack"].attrs["azimuthSpacing"] = stack.metadata[master_acqDate]["azimuthSpacing"]
        f["stack"].attrs["centerLat"] = stack.metadata[master_acqDate]["centerLat"]
        f["stack"].attrs["centerLon"] = stack.metadata[master_acqDate]["centerLon"]
        f["stack"].attrs["centerH"] = stack.metadata[master_acqDate]["centerH"]

def openHDF(HDFfile):
            pass

def toJSON(cr, outDir):

    outFile = outDir + "/" + "gf3_metadata.json"
    with open(outFile, "w") as file:
        json.dump(cr, file, indent=2)
    return outFile

def fromJSON(filePath):
    with open(filePath, "r") as file:
        data = json.load(file)
    return data

if __name__ == "__main__":
    parms = {
        "stackDir": '/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii',
        "stackId": "GF3",
        "outDir": "/data/tests/hongcan/GF3/inner_mongolia/test",
        "startDate": "20230119",
        "D_A_thr": 0.25,
        "master_name": '20230217'
    }
    path = "/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii"
    prepare_insar_HDF(parms, path) 
    print("Done.")