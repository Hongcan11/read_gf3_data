"""This module only concerns the reading part of the data.
"""
import os
# return a list of path(s) from the top directory
def get_gf3_path(top_dir):
    file_path = []
    for top_dir, dirs, files in os.walk(top_dir):           #顶层目录，文件夹，文件
         for dir in dirs:                                   #遍历目录下的文件夹
             file_path.append(os.path.join(top_dir, dir))   #获得子目录路径
         for file in files:                                 #遍历目录下的文件
             file_path = os.path.join(top_dir, file)        #获得文件的路径
         meta_paths = file_path[0]                          #得到 meta_paths
         slc_paths = file_path[1]                           #得到 slc_paths
    return meta_paths, slc_paths, dates

# return the meta data (in some sort of data structure, i.e., dictionary)
def read_gf3_meta(meta_paths):
    meta = []
    files = os.listdir(meta_paths)                     #得到文件夹下所有文件的名称
    for file in files:                                 #遍历文件夹
        if os.path.isfile(file):                       #判断是不是文件
           meta.append = ReadData(file)                #读取文件
        sub_path = get_gf3_path(meta_paths)[0]
        read_gf3_meta(sub_path)                        #迭代读取子文件夹下的文件
    return meta

# return the SLC data (either in memory or on disk)
def read_gf3_slc(gf3_meta, slc_path):
    slc = []
    files = os.listdir(meta_paths)           #得到文件夹下所有文件的名称
    for file in files:                       #遍历文件夹
        if os.path.isfile(file):             #判断是不是文件
            slc.append = ReadData(file)      #读取文件
        sub_path = get_gf3_path(slc_paths)[0]
        read_gf3_slc(sub_path)               #迭代读取子文件夹下的文件
    return slc

if __name__ == "__main__":
    meta_paths, slc_paths, dates = get_gf3_path(top_dir)
    for meta_path, slc_path, date in (meta_paths, slc_paths, dates):
      gf3_meta = read_gf3_meta(meta_path)
      gf3_slc = read_gf3_slc(gf3_meta, slc_path)