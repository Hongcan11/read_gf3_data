import os
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
    string
        paths of some data required.
    """  

    file_path = []       
    meta_paths = []      
    res_paths = []
    raw_paths = []
    cor_slaves_paths = []
    for top_dir, dirs, files in os.walk(top_dir):          #顶层目录，文件夹，文件        
        for dir in dirs:                                   #遍历目录下的文件夹
            file_path.append(os.path.join(top_dir, dir))   #获得子目录路径
    print('file_path', ':', file_path)                     #查看一下子目录
    process_paths = file_path[0]           
    slc_paths = file_path[1]               # 得到 process_paths，slc_paths
    process_sub_paths = file_path[2]       # process_paths 文件夹下的 /S01B01
    print('slc_path', ':', slc_paths)
    print('process_sub_paths', ':', process_sub_paths)

    for file in os.listdir(slc_paths):    # 遍历 SLC 文件路径
        if not os.path.isfile(os.path.join(slc_paths, file)):
            slaves_paths = os.path.join(slc_paths, file)       # 得到 SLC 下面两个文件夹的路径
            for file in os.listdir(slaves_paths):
                if file.endswith('meta.xml'):                  # 查找后缀名为 meta.xml 的文件
                    meta_paths.append(os.path.join(slaves_paths, file))     # 添加到 meta_paths 里面
    print('meta_paths', ':', meta_paths)

    for file in os.listdir(process_sub_paths):                #遍历 /process/S01B01 文件夹
        if not os.path.isfile(os.path.join(process_sub_paths, file)):
            cor_slaves_paths = os.path.join(process_sub_paths, file)
            print('cor_slaves_paths', ':', cor_slaves_paths)
            for file in os.listdir(cor_slaves_paths):
                if file.endswith('slave.res'):                 # 查找后缀名为 slave.res 的文件
                    res_paths.append(os.path.join(cor_slaves_paths, file))
                else:
                    continue
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

    print('res_paths', ':', res_paths)
    print('raw_paths', ':', raw_paths)
    return meta_paths, res_paths, raw_paths


if __name__ == "__main__":
    meta_paths= get_gf3_path('/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii')