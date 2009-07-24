// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XMASSFILE_H
#define OPENMS_FORMAT_XMASSFILE_H

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>
#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
 	/**
 		@brief File adapter for 'XMass Analysis (fid)' files.
 		test commit
 		MZ calculs are based on article :
      A database application for pre-processing, storage and comparison of mass spectra derived from patients and controls
      Mark K Titulaer, Ivar Siccama, Lennard J Dekker, Angelique LCT van Rijswijk, Ron MA Heeren, Peter A Sillevis Smitt, and Theo M Luider
      BMC Bioinformatics. 2006; 7: 403
      http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=1594579&blobtype=pdf
  	
  	@ingroup FileIO
  */
  
  class OPENMS_DLLAPI XMassFile
  {
// 	  private:
//		  PeakFileOptions options_;

    public:
      /// Default constructor
      XMassFile();
      
      template <class PeakType>
      void load(const String& filename, MSSpectrum<PeakType>& spectrum)
      {					        
        Internal::AcqusHandler acqus(filename.prefix(filename.length()-3) + String("acqus"));

        Internal::FidHandler fid(filename);
        if (!fid)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}      
				
        //  Delete old spectrum
				spectrum.clear();
				
				//temporary variables
				PeakType p;
							
        while( spectrum.size() < acqus.getSize() )
        {
				  //fill peak
				  p.setPosition( (typename PeakType::PositionType) acqus.getPosition(fid.getIndex()) );
				  p.setIntensity( (typename PeakType::IntensityType) fid.getIntensity() );
				  spectrum.push_back(p);
        }
        fid.close();
        
        // import metadata
        spectrum.setType(SpectrumSettings::RAWDATA);
        spectrum.setNativeID(acqus.getParam("ID_raw"));
        
        InstrumentSettings& instrumentSettings = spectrum.getInstrumentSettings();
          instrumentSettings.setScanMode(InstrumentSettings::MASSSPECTRUM);
          instrumentSettings.setZoomScan(false);
          if(acqus.getParam(".IONIZATION MODE") == "LD+")      
            instrumentSettings.setPolarity(IonSource::POSITIVE);
          else if(acqus.getParam(".IONIZATION MODE") == "LD-")      
            instrumentSettings.setPolarity(IonSource::NEGATIVE);
          else 
            instrumentSettings.setPolarity(IonSource::POLNULL);     
            
        // AcquisitionInfo& acquisitionInfo = spectrum.getAcquisitionInfo();
          
        SourceFile& sourceFile = spectrum.getSourceFile();
          sourceFile.setNameOfFile("fid");
          sourceFile.setPathToFile(filename.prefix(filename.length()-3));
          sourceFile.setFileSize(100.0);
          sourceFile.setFileType("Xmass analysis file (fid)");
          // sourceFile.setChecksum(const String &checksum, ChecksumType type);
          sourceFile.setNativeIDType(SourceFile::BRUKER_FID);
        
        // std::vector<Product>& product = spectrum.getProducts();
        // std::vector<Precursor>& precursor = spectrum.getPrecursors();
        // std::vector<PeptideIdentification>& peptideIdentification = spectrum.getPeptideIdentifications();
      }
        
      template <class PeakType>
      void importExperimentalSettings(const String& filename, MSExperiment<PeakType>& exp)
      {		        
        Internal::AcqusHandler acqus(filename.prefix(filename.length()-3) + String("acqus"));
        ExperimentalSettings& experimentalSettings = exp.getExperimentalSettings();
        
        Sample& sample = experimentalSettings.getSample();
          sample.setName("unknow");
          sample.setOrganism("unknow");
          sample.setNumber("unknow");
          sample.setComment("none");
          sample.setState(Sample::SAMPLENULL);          
          sample.setMass(0.0);
          sample.setVolume(0.0);
          sample.setConcentration(0.0);
        
        std::vector<SourceFile>& sourceFileList = experimentalSettings.getSourceFiles();
          sourceFileList.clear();
          for(typename MSExperiment<PeakType>::Iterator it=exp.begin(); it!=exp.end(); ++it)
          {
            sourceFileList.push_back(it->getSourceFile());
          }
                
        std::vector<ContactPerson>& contactPerson = experimentalSettings.getContacts();
          contactPerson.clear();
        
        Instrument& instrument = experimentalSettings.getInstrument();
          instrument.setName(acqus.getParam("SPECTROMETER/DATASYSTEM"));
          instrument.setVendor(acqus.getParam("ORIGIN"));
          instrument.setModel(acqus.getParam("$InstrID"));
          
          std::vector<IonSource>& ionSourceList = instrument.getIonSources();
            ionSourceList.clear();
            ionSourceList.resize(1);
            if(acqus.getParam(".INLET") == "DIRECT")
              ionSourceList[0].setInletType(IonSource::DIRECT);
            else
              ionSourceList[0].setInletType(IonSource::INLETNULL);
            ionSourceList[0].setIonizationMethod(IonSource::MALDI);
            if(acqus.getParam(".IONIZATION MODE") == "LD+")      
              ionSourceList[0].setPolarity(IonSource::POSITIVE);
            else if(acqus.getParam(".IONIZATION MODE") == "LD-")      
              ionSourceList[0].setPolarity(IonSource::NEGATIVE);
            else 
              ionSourceList[0].setPolarity(IonSource::POLNULL);
            ionSourceList[0].setOrder(0);
            
          std::vector<MassAnalyzer>& massAnalyzerList = instrument.getMassAnalyzers();
            massAnalyzerList.clear();
            massAnalyzerList.resize(1);  
            if(acqus.getParam(".SPECTROMETER TYPE") == "TOF")
              massAnalyzerList[0].setType(MassAnalyzer::TOF);
            else
              massAnalyzerList[0].setType(MassAnalyzer::ANALYZERNULL);
                               
        // HPLC& hplc = experimentalSettings.getHPLC();
        
        DateTime date;
          date.set(acqus.getParam("$AQ_DATE").remove('<').remove('>') );
          experimentalSettings.setDateTime(date);
        
        String comment("comment");
          experimentalSettings.setComment(comment);
        
        // std::vector<ProteinIdentification>& proteinIdentification = experimentalSettings.getProteinIdentifications();
        

      }

      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("'Xmass' file not writable."));
      }
      
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_XMASSFILE_H

