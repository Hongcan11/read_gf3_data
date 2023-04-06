import numpy as np

def read_gf3_slc(slc_path):
    with open(slc_path, "rb") as f:        
        # slc = np.fromfile(f, dtype = np.complex64)
        slc = np.fromfile(f, dtype = np.float16)
        slc = complex_short2complex(slc)
        shape = (26665, 28238)
        slc = np.reshape(slc, shape)
        print(slc.shape, slc.dtype, slc)
    return slc


def complex_short2complex(data):
        """
        Convert complex short values to regular numpy complex64 data.

        :param np.ndarray data: Input data to be converted
        :return: data
        """
        data = np.copy(data).view(np.float16).astype('float32').view(np.complex64)
        return data

if __name__ == "__main__":
    read_gf3_slc("/data/tests/hongcan/GF3/cn_inner_mongolia_gf3c_dsc_fsii/process/S01B01/master_crop.raw")