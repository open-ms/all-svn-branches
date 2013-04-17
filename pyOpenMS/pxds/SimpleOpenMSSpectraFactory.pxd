from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from ISpectrumAccess cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>" namespace "OpenMS":

    cdef cppclass SimpleOpenMSSpectraFactory(ProgressLogger):

        SimpleOpenMSSpectraFactory() nogil except +

        # shared_ptr[ISpectrumAccess] getSpectrumAccessOpenMSPtr(MSExperiment[Peak1D,ChromatogramPeak] exp) # wrap-ignore
        # OPENSWATHALGO_DLLAPI typedef boost::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;


