"""This module only concerns the reading part of the data.
"""
import os
import glob
# return a list of path(s) from the top directory
def get_gf3_path(top_dir):
    """
    Function for getting meta_paths and slc_paths from top_dir. 

    Parameters
    ----------
    top_dir : string
        Top path of gf3_data.

    Returns
    -------
    meta_paths : string
    slc_paths : string

    """    
    file_path = []
    data_path = []
    for top_dir, dirs, files in os.walk(top_dir):           #顶层目录，文件夹，文件
         for dir in dirs:                                   #遍历目录下的文件夹
             file_path.append(os.path.join(top_dir, dir))   #获得子目录路径
         for file in files:                                 #遍历目录下的文件
             file_path = os.path.join(top_dir, file)        #获得文件的路径
         process_paths = file_path[0]                       #得到 process_paths
         slc_paths = file_path[1]                           #得到 slc_paths 
         for file in os.listdir(process_paths):
             if not os.path.isfile(file):
                 sub_paths.append(os.path.join(process_paths, file))
         for file in os.listdir(slave_paths):  
             if not os.path.isfile(file):
                 slave_paths.append(os.path.join(sub_paths, file))
         for file in os.listdir(slave_paths):
             if os.path.isfile(file):
                 meta_paths = glob.glob('*.xml')            #得到 meta_paths .xml文件的地址
                 data_paths = glob.glob('*.tiff')           #得到 tiff 文件的地址        
    return meta_paths, slc_paths, dates

# return the meta data (in some sort of data structure, i.e., dictionary)
def read_gf3_meta(meta_paths):
    """
    Function for reading gf3_metadata

    Parameters
    ----------
    meta_paths : string

    Returns
    -------
    meta : dictionary
         Return a dictionary storing the metadata
    """    
    meta = {}
    #get required metadata
    get_orbit;
    get_acqDate;
    get_resolutions;
    ...

    # create meta_data dictionary
    meta = {...}
    return meta

# return the SLC data (either in memory or on disk)
def read_gf3_slc(slc_path):
    """
    Function for reading gf3_slc

    Parameters
    ----------
    slc_path : string

    Returns
    -------
    slc : array(range, azimuth)
         Returns an array storing slc_data
    """    
    slc = []     #array, 根据 boundingBox 给定数组的大小
    files = os.listdir(meta_paths)           #得到文件夹下所有文件的名称
    for file in files:                       #遍历文件夹
        if not os.path.isfile(file):         #判断是不是文件
            sub_path = os.path.join(meta_paths, file)
            slc.append = ReadSLC(sub_path)   #读取文件
        
    return slc

if __name__ == "__main__":
    meta_paths, slc_paths, dates = get_gf3_path(top_dir)
    for meta_path, slc_path, date in (meta_paths, slc_paths, dates):
      gf3_meta = read_gf3_meta(meta_path)
      gf3_slc = read_gf3_slc(gf3_meta, slc_path)