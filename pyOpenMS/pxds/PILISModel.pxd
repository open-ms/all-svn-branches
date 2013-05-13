from Types cimport *
from libcpp cimport bool
from Map cimport *
from Types cimport *
from HiddenMarkovModel cimport *
from ProtonDistributionModel cimport *
from TheoreticalSpectrumGenerator cimport *
# from SpectrumAlignment cimport *
from TextFile cimport *
from DefaultParamHandler cimport *
from PILISNeutralLossModel cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PILISModel.h>" namespace "OpenMS":
    
    cdef cppclass PILISModel(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PILISModel() nogil except +
        PILISModel(PILISModel) nogil except +
        void train(MSSpectrum[RichPeak1D] & , AASequence & peptide, UInt charge) nogil except +
        void readFromFile(String & filename) nogil except +
        void writeGraphMLFile(String & filename) nogil except +
        void writeToFile(String & filename) nogil except +
        void init(bool generate_models) nogil except +
        void getSpectrum(MSSpectrum[RichPeak1D] & spec, AASequence & peptide, UInt charge) nogil except +
        void evaluate() nogil except +

