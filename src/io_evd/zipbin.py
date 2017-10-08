import cPickle
import os
import zlib
import warnings

class Stream:
    cmpLevel = 9
    
    def __init__(self, fileName=None, cmpLevel = 9):
        if fileName == None:
            raise IOError("Valid file name not provided")
        fileName = os.path.abspath(fileName)
        dirName = os.path.dirname(fileName)
        if not os.access(dirName, os.W_OK):
            raise IOError("No write privileges to directory '%s'" % dirName)
        if not os.path.isdir(dirName):
            raise IOError("Directory '%s' does not exist" % (dirName))
        
        if not isinstance(cmpLevel, int):
            raise IOError("cmpLevel must be an integer between 0-9")
        if (cmpLevel < 0) or (cmpLevel > 9):
            raise IOError("cmpLevel must be an integer between 0-9")
        
        self._file = fileName
        self._cmpLevel = cmpLevel
        
    def write(self, obj):
        if os.path.exists(self._file):
            warnings.warn("Output file already exist and will be overwritten.", UserWarning)
        
        datString = cPickle.dumps(obj, protocol=cPickle.HIGHEST_PROTOCOL)
        cmpString = zlib.compress(datString, self._cmpLevel)
        with open(self._file, "wb") as f:
            f.write(cmpString)
            
    def read(self):
        if not os.path.exists(self._file):
            raise IOError("Input file does not exist")
        with open(self._file, "rb") as f:
            cmpString = f.read()
        datString = zlib.decompress(cmpString)
        obj = cPickle.loads(datString)
        return obj
            