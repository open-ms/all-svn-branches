from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>" namespace "OpenMS":

  cdef cppclass SpectrumAccessOpenMS:
        SpectrumAccessOpenMS(MSExperiment[Peak1D, ChromatogramPeak] & ms_experiment)
        SpectrumAccessOpenMS(SpectrumAccessOpenMS)

        shared_ptr[Spectrum] getSpectrumById(int id)  #wrap-ignore
        libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT)
        size_t getNrSpectra()



