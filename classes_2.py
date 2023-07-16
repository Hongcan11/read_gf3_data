from collections import OrderedDict
import numpy as np
from gf3_ingest_module import read_metadata, get_gf3_path

class Stack:
    def __init__(self, id, stackDir):
        self.id = id
        self.stackDir = stackDir
        # self.sensor = sensor
        # self.subswath = subswath
        # self.type = stackType
        self._nSLC = None
        self._masterIdx = None

    def readData(self):
        d = dict()
        gf3_meta = []
        path = get_gf3_path(self.stackDir)[0]
        print("path: ", path)      
        for meta_path in path:
            # pattern1 = r"_(\d{8})_"
            # date = get_date(meta_path, pattern1)
            gf3_meta = read_metadata(meta_path)
            acqStr = gf3_meta["acqDate"] 
            d[acqStr] = gf3_meta.copy()
        print("metadata: ", d)
        self.metadata = OrderedDict(sorted(d.items(), key=lambda t: t[0]))
        self.acqDates = list(self.metadata.keys())
        self.files = sorted(path)
        self.orbit = gf3_meta["orbit"]


    def reduce(self, startDate):
        """
        startDate in 'YYYYMMDD' format
        """
        if hasattr(self, "acqDates"):
            # remove dates before startDate:
            start_idx = next(
                x for x, val in enumerate(self.acqDates) if val > startDate
            )
            if self.masterDate not in self.acqDates[start_idx:]:
                print("Start date specified after master date. Quiting.")
                raise
            else:
                self.acqDates = self.acqDates[start_idx:]
            # adapt files and metadata dicts:
            self.files = self.files[start_idx:]
            tmp = dict()
            for k, v in self.metadata.items():
                if k > startDate:
                    tmp[k] = v
            self.metadata = tmp
        else:
            print("Perform .readData() first!")

    def remove_idx(self, indices):
        if hasattr(self, "acqDates"):
            if self.masterIdx in indices:
                print("Cannot delete master date!")
                raise
            else:
                removeDates = np.array(self.acqDates)[np.array(indices)]
                self.acqDates = [
                    j for i, j in enumerate(self.acqDates) if i not in indices
                ]
                self.files = [
                    j for i, j in enumerate(self.files) if i not in indices
                ]
                for k in removeDates:
                    del self.metadata[k]

    def remove_date(self, acqDates):
        if hasattr(self, "acqDates"):
            remove_idx = [self.acqDates.index(d) for d in acqDates]
            self.remove_idx(remove_idx)

    @property
    def nSLC(self):
        if self._nSLC is None:
            if hasattr(self, "files"):
                self._nSLC = len(self.files)
            else:
                print("Perform .readData() first!")
        return self._nSLC

    @property
    def masterIdx(self):
        # if self._masterIdx is None:
        if hasattr(self, "files"):
            self._masterIdx = np.where(
                np.array(self.acqDates) == self.masterDate
            )[0][0]
        else:
            print("Perform .readData() first!")
        return self._masterIdx
