import os.path
import csv
import warnings

class Data:
    
    validFileExt = (".TXT", ".CSV", ".TSV")
    validNAs     = ("", "NA", "NAN")
    maxIter = 10e6
    sniffRange = 8192
    
    def __init__(self, inFile=None, outFile=None, inNAsAllowed=False, outNAsOnly=True, firstLine=None):
        self._inValid = False
        self._outValid = False
        self._loaded = False
        self._dialect = csv.excel_tab
        self._dialect_fixed = False
        self._dictData = None
        
        self._inNAsAllowed = inNAsAllowed
        self._outNAsOnly = outNAsOnly
        
        self._firstLine = firstLine
        
        self._inputFields = None
        self._outputFields = None
        self._allFields = None

        
        self._iterI = []
        self._iterInLoop = False
        self._numIters = 0
        self._currIterI = 0
        

        self._inFile = rel_2_abs_path(inFile)
        self._outFile = rel_2_abs_path(outFile)
        
        self._check_input_valid()
        self._check_output_valid()
                 
    def _check_input_valid(self):
        if self._inFile == None:
            self._inValid = False
            return
        if os.path.isfile(self._inFile):
            _, ext = os.path.splitext(self._inFile)
            if ext.upper() in Data.validFileExt:
                self._inValid = True
            else:
                self._inValid = False
                raise IOError("Input file extension '%s' not supported" % ext)
        else:
            self._inValid = False
            raise IOError("Input file '%s' does not exist" % self._inFile)
            
    def _check_output_valid(self):
        self._outValid = False
        if self._outFile == None:
            return
        dirName = os.path.dirname(self._outFile)
        _, ext = os.path.splitext(self._outFile)
        if not (ext.upper() in Data.validFileExt):
            raise IOError("Output file extension '%s' not supported" % ext)
        if not os.access(dirName, os.W_OK):
            raise IOError("No write privileges to directory '%s'" % dirName)
        if os.path.exists(self._outFile):
            warnings.warn("Output file already exist and will be overwritten.", UserWarning)
        self._outValid = True
        
    def _get_NA_Is(self, fields):
        naIs = set()
        if fields == None or not fields:
            return naIs
        for i, item in enumerate(self._dictData):
            for field in fields:
                if not field in item.keys():
                    raise TypeError("Input or output fields not properly instantiated")
                value = item[field]
                if value.upper() in Data.validNAs:
                    naIs.add(i)
                    break
        return naIs
    
    def _filter_NAs(self):
        inNaIs  = self._get_NA_Is(self._inputFields)
        outNaIs = self._get_NA_Is(self._outputFields)
        
        allIs = set(range(len(self._dictData)));
        
        if self._inNAsAllowed and not self._outNAsOnly:
            iterI = allIs
        elif self._inNAsAllowed and self._outNAsOnly:
            iterI = outNaIs
        elif not self._inNAsAllowed and not self._outNAsOnly:
            iterI = allIs - inNaIs
        else:
            iterI = outNaIs - inNaIs  
        self._iterI = list(iterI)
        
    def _check_dictData_is_consistant(self, dictData):
        if not (type(dictData) is list):
            raise TypeError("Data must be presented as a list of dictionaries")
        for i, item in enumerate(dictData):
            if not (type(item) == dict):
                raise TypeError("Each row entry must be in dictionary format")
            if i == 0:
                fields = item.keys()
            for field in fields:
                if not field in item:
                    raise TypeError("Each row entry must contain all fields")
                if not isinstance(item[field], basestring):
                    raise TypeError("Each cell entry must be in string format")
        
    def __iter__(self):
#         if self._inputFields == None or not self._inputFields:
#             raise IOError("Input fields not specified")
        
        self._filter_NAs()
   
        self._numIters = len(self._iterI)
        self._currIterI = -1
        self._iterInLoop = True
        return self
    
    def next(self):
        self._currIterI += 1
        if (self._currIterI == self._numIters) or (self._currIterI == Data.maxIter):
            self._iterInLoop = False
            raise StopIteration
        currData = self._dictData[self._iterI[self._currIterI]]
        inList = [currData[fieldKey] for fieldKey in self._inputFields]
        return inList
    
    def set_dialect(self, dialect):
        self._dialect = dialect
        self._dialect_fixed = True
        
    def set_in_file(self, inFile):
        inFile = rel_2_abs_path(inFile)
        if (self._inValid == True) and (self._inFile != inFile):
            raise IOError("Input file can only be assigned once and is set to: '%s'" % self._inFile)
        
        self._inFile = inFile
        self._check_input_valid()
    
    def set_out_file(self, outFile):
        outFile = rel_2_abs_path(outFile)
        if (self._outValid == True) and (self._outFile != outFile):
            raise IOError("Output file only allowed to be set once and is set to: '%s'" % self._outFile)
        self._outFile = outFile
        self._check_output_valid()
        
    def get_in_file(self):
        return self._inFile
    
    def get_out_file(self):
        return self._outFile
    
    def set_data(self, dictData, **kwargs):
        print kwargs.keys()
        
        if self._iterInLoop:
            raise IOError("Not allowed to change data dictionary while iterating")
        if self._loaded == True:
            raise IOError("Only allowed to load data once. All updates allowed only through 'write_out_fields()'")
        self._check_dictData_is_consistant(dictData)
        self._dictData = dictData[:]
        if "keyList" in kwargs.keys():
            self._allFields = kwargs["keyList"]
        else:
            self._allFields = self._dictData[0].keys()
        self._loaded = True
        
    def get_data(self):
        if not self._loaded:
            self.read()
        return self._dictData

    def get_field(self, fieldName):
        if not self._loaded:
            self.read()
        fieldEntries = list()
        for item in self._dictData:
            if fieldName in item.keys():
                fieldEntries.append(item[fieldName])
        return fieldEntries
        
    def get_field_names(self):
        if not self._loaded:
            self.read()
        fieldNames = self._allFields
        return fieldNames
        
    def read(self):
        if self._iterInLoop:
            raise IOError("Not allowed to read input file while iterating")
        if not self._inValid:
            raise IOError("No valid input file to read")
        if self._loaded:
            raise IOError("Only allowed to load data once. All updates allowed only through 'write_out_fields()'")
            
        dictData = []
        with open(self._inFile, "rU") as f:
            if not self._firstLine == None:
                for ignoreLine in range(self._firstLine):
                    ignoreText = f.readline()
            if not self._dialect_fixed:
                self._dialect = csv.Sniffer().sniff(f.read(Data.sniffRange))
                f.seek(0)
                if not self._firstLine == None:
                    for ignoreLine in range(self._firstLine):
                        ignoreText = f.readline()
                
            dictReader = csv.DictReader(f, dialect=self._dialect)
            self._allFields = dictReader.fieldnames
            for row in dictReader:
                dictData.append(row)
        self._check_dictData_is_consistant(dictData)
        self._dictData = dictData
        self._loaded = True
        
#         self.set_out_fields(self._outputFields)
        
    def write(self):
        if not self._outValid:
            raise IOError("No valid output file to write to")
        if not self._loaded:
            self.read()
        if len(self._dictData) == 0:
            warnings.warn("Nothing to write. Output file is not created.", UserWarning)
            return
        if self._allFields == None:
            self._allFields = self._dictData[0].keys()
        with open(self._outFile, "wb") as f:
            dictWriter = csv.DictWriter(f, dialect=self._dialect, fieldnames=self._allFields)
            dictWriter.writeheader()
            for row in self._dictData:
                dictWriter.writerow(row)
    
    def set_in_fields(self, fields=[]):
        if self._iterInLoop:
            raise IOError("Not allowed to change input fields while iterating")
        fields = typeset2strList(fields)
        if not self._loaded:
            self.read()
        if len(self._dictData) == 0:
            return
        fieldsExist = [0 for field in fields]
        for i in range(len(self._dictData)):
            for j, field in enumerate(fields):
                if field in self._dictData[i].keys():
                    fieldsExist[j] += 1           
        for j, field in enumerate(fields):
            if not fieldsExist[j]:
                raise TypeError("Input field '%s' does not exist" % field)
        self._inputFields = fields      
        
    def allow_input_NAs(self):
        if self._iterInLoop and not self._inNAsAllowed:
            raise IOError("Cannot suddenly allow input NAs while iterating")
        self._inNAsAllowed = True
             
    def ignore_input_NAs(self):
        if self._iterInLoop and self._inNAsAllowed:
            raise IOError("Cannot suddenly ignore input NAs while iterating")
        self._inNAsAllowed = False
            
    def set_out_fields(self, fields=[]):
        if self._iterInLoop:
            raise IOError("Not allowed to change output fields while iterating")
        fields = typeset2strList(fields)
        self._outputFields = fields
        if not self._loaded:
            self.read()
        if fields == None or not fields:
            return
        if len(self._dictData) == 0:
            return   
    
        for i in range(len(self._dictData)):
            for field in fields:
                if not field in self._allFields:
                    if self._allFields == None:
                        self._allFields = [field]
                    else:
                        self._allFields.append(field)
                if not field in self._dictData[i].keys():
                    self._dictData[i][field] = Data.validNAs[0]

    def set_fields(self, fields=[]):
        if self._iterInLoop:
            raise IOError("Not allowed to change output fields while iterating")
        fields = typeset2strList(fields)
        if not self._loaded:
            self.read()
        if fields == None or not fields:
            return
        if len(self._dictData) == 0:
            return   
        for i in range(len(self._dictData)):
            for field in fields:
                if not field in self._allFields:
                    if self._allFields == None:
                        self._allFields = [field]
                    else:
                        self._allFields.append(field)
                if not field in self._dictData[i].keys():
                    self._dictData[i][field] = Data.validNAs[0]
        self._allFields = fields
                    
    def overwrite_defined_outputs(self):
        if self._iterInLoop and self._outNAsOnly:
            raise IOError("Cannot suddenly allow to overwrite defined outputs while iterating")
        self._outNAsOnly = False
             
    def ignore_defined_outputs(self):
        if self._iterInLoop and not self._outNAsOnly:
            raise IOError("Cannot suddenly allow to overwrite defined outputs while iterating")
        self._outNAsOnly = True
                
    def write_out_fields(self, fields=[]):
        if not self._iterInLoop:
            raise IOError("Writing output fields are only permitted while iterating through items")
        if self._outputFields == None or not self._outputFields:
            raise IOError("Output fields not set.")
        fields = typeset2strList(fields)
        if not len(fields) == len(self._outputFields):
            raise TypeError("The number of output fields do not match those set by 'set_out_fields()'")
        
        for i, field in enumerate(self._outputFields):
            self._dictData[self._iterI[self._currIterI]][field] = fields[i]
            
    def append_out_fields(self, fields=[]):
        if self._iterInLoop:
            raise IOError("Appending output fields are not permitted while iterating through items")
        if self._outputFields == None or not self._outputFields:
            raise IOError("Output fields not set.")
        fields = typeset2strList(fields)
        if not len(fields) == len(self._outputFields):
            raise TypeError("The number of output fields do not match those set by 'set_out_fields()'")
        if not self._loaded:
            self.read()
        currDict = dict()
        for field in self._allFields:
            currDict[field] = Data.validNAs[0]  
        for i, field in enumerate(self._outputFields):
            currDict[field] = fields[i]   
        self._dictData.append(currDict)
        
def tsvDict2fieldDict(tsvDict):
    numEntries = len(tsvDict)
    if numEntries < 1:
        raise IOError("tsvDict must contain at least one entry")
    fields = tsvDict[0].keys()
    fieldDict = dict()
    for currField in fields:
        fieldDict[currField] = list()
    for currEntry in range(numEntries):
        for currField in fields:
            fieldDict[currField].append(tsvDict[currEntry][currField])
    return fieldDict

def fieldDict2tsvDict(fieldDict):
    tsvDict = list()
    keys = fieldDict.keys()
    fieldLen = 0
    for key in keys:
        fieldLen = max(fieldLen, len(fieldDict[key]))
    
    for i in range(fieldLen):
        currDict = dict()
        for key in keys:
            try:
                currDict[key] = str(fieldDict[key][i])
            except IndexError:
                currDict[key] = ""
                    
        tsvDict.append(currDict)
    return tsvDict   

def rel_2_abs_path(fileName):
    try:
        absFileName = os.path.abspath(fileName)
    except AttributeError:
        absFileName = None
    return absFileName
                     
def typeset2strList(fields):
    if isinstance(fields, basestring):
            fields = [fields]
    elif not isinstance(fields, list):
        raise TypeError("Output fields are not passed in a list format")
     
    listItemsAreStrings = True
    for value in fields:
        if not isinstance(value, basestring):
            listItemsAreStrings = False
            break
    if not listItemsAreStrings:
        raise TypeError("Not all output fields are in string format")           
    return fields  
    