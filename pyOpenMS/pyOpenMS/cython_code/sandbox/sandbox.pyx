
from cython.operator cimport dereference  as deref
from libcpp.vector cimport vector as cppvector
from libcpp.string cimport string as cppstring
#cimport vector

from pxd.DataValue cimport DataValue as _DataValue
from pxd.String cimport String as _String
from pxd.DoubleList cimport DoubleList as _DoubleList
from pxd.FileHandler cimport FileHandler as _FileHandler
from pxd.MSExperiment cimport MSExperiment as _MSExperiment
from pxd.Peak1D cimport Peak1D as _Peak1D
from pxd.ChromatogramPeak cimport ChromatogramPeak as _ChromatogramPeak
from pxd.Feature cimport Feature as _Feature
from pxd.FeatureMap cimport FeatureMap as _FeatureMap
from pxd.TransformationDescription cimport TransformationDescription as _TransformationDescription
from pxd.MapAlignmentAlgorithmPoseClustering cimport MapAlignmentAlgorithmPoseClustering as _MapAlignmentAlgorithmPoseClustering

cdef cppvector[float] py_list_to_vector_float(li):
    cdef cppvector[float]  res  
    for item in li:
        res.push_back(<float>item)
    return res

cdef cppvector[double] py_list_to_vector_double(li):
    cdef cppvector[double]  res  
    for item in li:
        res.push_back(<double>item)
    return res

cdef cppvector[_String] py_list_to_vector__String(li):
    cdef cppvector[_String]  res  
    for item in li:
        arg0 = py_str_to__String(item)
        res.push_back(arg0)
    return res

cdef vector__String_to_py_list(cppvector[_String] arg):
    return [ <cppstring>i for i in arg ]


#        cdef cppvector[unsigned int] arg0 = py_list_to_vector_unsigned_int(keys)
#        self.inst.getKeys(arg0)
#        return vector_unsigned_int_to_py_list(arg0)

cdef _DataValue py_object_to__DataValue(obj):
    if isinstance(obj, int):
        return _DataValue(<int> obj)
    if isinstance(obj, float):
        return _DataValue(<double> obj)
    if isinstance(obj, str):
        return _DataValue(<cppstring> obj)
    raise Exception("can not convert type %s to DataValue arg" % type(obj))


cdef DataValue_to_py(_DataValue val):
    type = val.valueType()
    if type == DataType.STRING_VALUE:
        return <cppstring> val
    elif type == DataType.INT_VALUE:
        return <int> val
    elif type == DataType.DOUBLE_VALUE:
        return <double> val
    elif type == DataType.EMPTY_VALUE:
        return None
    raise Exception("can not handle type %d" % type)
    

cdef cppvector[unsigned int] py_list_to_vector_unsigned_int(li):
    cdef cppvector[unsigned int]  res  
    for item in li:
        arg0 = <unsigned int>item
        res.push_back(arg0)
    return res
    
cdef vector_unsigned_int_to_py_list(cppvector[unsigned int] arg):
    return [ <int>i for i in arg ]
cdef class FeatureMap:
    cdef _FeatureMap[_Feature] * inst
    def __cinit__(self):
        self.inst = new _FeatureMap[_Feature]()
    def __dealloc__(self):
        del self.inst



cdef class str_conv:

    cdef _String * inst

    def __cinit__(self, str val):
        self.inst = new _String(<cppstring> val)

    cdef _String value(self):
        return deref(self.inst)

    def __dealloc__(self):
        del self.inst

cdef _String py_str_to__String(str val):
    return _String(<cppstring> val)
        

cdef class TransformationDescription:
    cdef _TransformationDescription * inst
    def __cinit__(self):
        self.inst = new _TransformationDescription()
    def __dealloc__(self):
        del self.inst

cdef class MapAlignmentAlgorithmPoseClustering:
    cdef _MapAlignmentAlgorithmPoseClustering * inst
    def __cinit__(self):
        self.inst = new _MapAlignmentAlgorithmPoseClustering()
    def __dealloc__(self):
        del self.inst

    def alignFeatureMaps(self, fms, tds=None):
        if type(fms) != list or any(not isinstance(fm, FeatureMap) for fm in fms):
            raise Exception()
        if tds is None:
            tds = []    
        if type(tds) != list or any(not isinstance(td, TransformationDescription) for td in tds):
            raise Exception()

        #cdef arg0 = py_list_to_vector__Feature_Map(fms)
        #cdef arg1 = py_list_to_vector__TransformationDescription(tds)
        #self.inst.alignFeatureMaps(arg0, arg1)
        
        #return vector__TransformationDescription_to_py(tds)



cdef class DataType:
    STRING_VALUE=0
    INT_VALUE=1
    DOUBLE_VALUE=2
    STRING_LIST=3
    INT_LIST=4
    DOUBLE_LIST=5
    EMPTY_VALUE=6

cdef class DataValue:
    cdef _DataValue * inst

    def __cinit__(self, val):
        try:
            if type(val) == int:
                self.inst = new _DataValue(<int>val)
            elif type(val) == float:
                self.inst = new _DataValue(<double>val)
            elif type(val) == DoubleList:
                self.inst = new _DataValue(deref((<DoubleList>val).inst))
            elif type(val) == list and all(type(v) == float for v in val):
                self.inst = new _DataValue(py_list_to_vector_float(val))
            else:
                raise Exception("input args do not match")
        except:
            raise
            

    def __dealloc__(self):
        del self.inst

    def toInt(self):
        return <int>deref(self.inst)

    def toDouble(self):
        return <double>deref(self.inst)

    def valueType(self):
        return self.inst.valueType()



cdef class DoubleList:

    cdef _DoubleList * inst

    def __cinit__(self, val):
        if type(val) == list and all(type(v) == float for v in val):
            #py_list_to_vector_double(val)
            self.inst = new _DoubleList((py_list_to_vector_double(val)))
        else:
            raise Exception("argument type does not match")

    def __dealloc__(self):
        del self.inst

    def at(self, int arg0):
        return self.inst.at(arg0)


cdef class FileHandler:

    cdef _FileHandler * inst

    def __cinit__(self):
        self.inst = new _FileHandler()

    def __dealloc__(self):
        del self.inst

    def loadExperiment(self, str path, MSExperiment exp):
        self.inst.loadExperiment(<cppstring>path, deref(exp.inst))



cdef class MSExperiment:

    cdef _MSExperiment[_Peak1D, _ChromatogramPeak] * inst

    def __cinit__(self):
        self.inst = new _MSExperiment[_Peak1D, _ChromatogramPeak]()

    def __dealloc__(self):
        del self.inst

    def getMetaValue(self, str name):
        #cdef str_conv arg0x = str_conv(name)
        cdef _String arg0 = py_str_to__String(name)
        cdef _DataValue res = self.inst.getMetaValue(arg0)
        return DataValue_to_py(res)

    def setMetaValue(self, str name, val):
        cdef _String arg0 = py_str_to__String(name)
        cdef _DataValue arg1 = py_object_to__DataValue(val)
        self.inst.setMetaValue(arg0, arg1)

    def getStringKeys(self, keys=None):
        if keys is None:
            keys = []
        if type(keys) != list or any(type(v) != str for v in keys):
            raise Exception("arg keys has wrong type")
        cdef cppvector[_String] arg0 = py_list_to_vector__String(keys)
        self.inst.getKeys(arg0)
        return vector__String_to_py_list(arg0)

    def getIntKeys(self, keys=None):
        if keys is None:
            keys = []

        if type(keys) != list or any(type(v) != int for v in keys):
            raise Exception("arg keys has wrong type")
        
        cdef cppvector[unsigned int] arg0 = py_list_to_vector_unsigned_int(keys)
        self.inst.getKeys(arg0)
        return vector_unsigned_int_to_py_list(arg0)

    def isMetaEmpty(self):
        return <int> self.inst.isMetaEmpty()
        

