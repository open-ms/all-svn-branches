// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Import settings
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

#ifndef OPENMS_FORMAT_SETTINGSPARSER_H
#define OPENMS_FORMAT_SETTINGSPARSER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <fstream>

using std::cout;
using std::endl;

namespace OpenMS
{
 	/**
  */
  
  class OPENMS_DLLAPI SettingsParser
  {
    private:
      int findIndex_(const String& value, const std::string* table, const Size tableSize) const;
      void stringToIntVector_(const String& str, std::vector<Int>& list) const;
      void addMetaInfo_(MetaInfoInterface& metaInfo, const String& name, const String& value) const;
      
      template <class T>
      void resizeList_(std::vector<T>& list, const int size) const
      {
        list.resize(size<0 ? Size(0) : Size(size));
      }
      
      template <class T>
      void setIndex_(Size index, const int value, std::vector<T>& list) const
      {
        index = value<1 ? Size(0) : Size(value-1);
        if(index >= list.size())
          list.resize(index + 1);
      }
      
    public:
      /// Default constructor
      SettingsParser();
      
      template <class PeakType>
      void importSpectrumSettings(const String& line, MSSpectrum<PeakType>& spectrum) const
      {
        //temporary variables
        String command = "";
        String value = "";
        SourceFile::ChecksumType checksum_type_ = SourceFile::UNKNOWN_CHECKSUM;
        String checksum_ = "";
        AASequence aaSequence_;
        Size indexScanWindows = 0;
        Size indexPrecursor = 0;
        Size indexProduct = 0;
        Size indexPeptideIdentification = 0;
        Size indexPeptideHit = 0;

        try
        {            
          {
            std::vector<String> strings;
            if(!line.split('=', strings))
              throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" );
              
            if( strings.size() != 2 )
            {
cout << "BBB: " << strings.size() << endl;
              throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" );
            }
              
            command = strings[0].trim();
            value = strings[1].trim();
          }    
                 
          if(command.hasPrefix("#"))
            return;          
          
          if(command.hasPrefix("SpectrumSettings."))
          {
            if(command.hasPrefix("SpectrumSettings.Type"))
              spectrum.setType((SpectrumSettings::SpectrumType) findIndex_(
                value, SpectrumSettings::NamesOfSpectrumType, SpectrumSettings::SIZE_OF_SPECTRUMTYPE));
            else if(command.hasPrefix("SpectrumSettings.NativeID"))
              spectrum.setNativeID(value);
            else if(command.hasPrefix("SpectrumSettings.Comment"))
              spectrum.setComment(value);
            else if(command.hasPrefix("SpectrumSettings.InstrumentSettings."))
            {
              InstrumentSettings& instrumentSettings = spectrum.getInstrumentSettings();
              if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ScanMode"))
                instrumentSettings.setScanMode((InstrumentSettings::ScanMode) findIndex_(
                  value, InstrumentSettings::NamesOfScanMode, InstrumentSettings::SIZE_OF_SCANMODE));
              else if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ZoomScan"))
                instrumentSettings.setZoomScan(value=="TRUE" ? true : false);
              else if(command.hasPrefix("SpectrumSettings.InstrumentSettings.Polarity"))
                instrumentSettings.setPolarity((IonSource::Polarity) findIndex_(
                  value, IonSource::NamesOfPolarity, IonSource::SIZE_OF_POLARITY));
              else if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ScanWindows."))
              {
                std::vector<ScanWindow>& scanWindowsList = instrumentSettings.getScanWindows();
                
                if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ScanWindows.Size"))
                  resizeList_(scanWindowsList, value.toInt());
                else if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ScanWindows.Index"))
                  setIndex_(indexScanWindows, value.toInt()-1, scanWindowsList);              
                else if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ScanWindows.Begin"))
                  scanWindowsList[indexScanWindows].begin = value.toDouble();
                else if(command.hasPrefix("SpectrumSettings.InstrumentSettings.ScanWindows.End"))
                  scanWindowsList[indexScanWindows].end = value.toDouble();
                else
                  addMetaInfo_(/*(MetaInfoInterface&)*/ scanWindowsList[indexScanWindows], command, value);
              }                    
            }
            else if(command.hasPrefix("SpectrumSettings.AcquisitionInfo.MethodOfCombination"))
            {
                AcquisitionInfo& acquisitionInfo = spectrum.getAcquisitionInfo();
                acquisitionInfo.setMethodOfCombination(value);
            }         
            else if(command.hasPrefix("SpectrumSettings.SourceFile."))
            {
              SourceFile& sourceFile = spectrum.getSourceFile();
              
              if(command.hasPrefix("SpectrumSettings.SourceFile.NameOfFile"))
                sourceFile.setNameOfFile(value);
              else if(command.hasPrefix("SpectrumSettings.SourceFile.PathToFile"))
                sourceFile.setPathToFile(value);
              else if(command.hasPrefix("SpectrumSettings.SourceFile.FileSize"))
                sourceFile.setFileSize(value.toFloat());
              else if(command.hasPrefix("SpectrumSettings.SourceFile.FileType"))
                sourceFile.setFileType(value);
              else if(command.hasPrefix("SpectrumSettings.SourceFile.ChecksumType"))
              {
                checksum_type_ = findIndex_(value, SourceFile::NamesOfChecksumType, SourceFile::SIZE_OF_CHECKSUMTYPE);
                sourceFile.setChecksum(checksum_, checksum_type_);
              }
              else if(command.hasPrefix("SpectrumSettings.SourceFile.Checksum"))
              {
                checksum_ = value;
                sourceFile.setChecksum(checksum_, checksum_type_);
              }
              else if(command.hasPrefix("SpectrumSettings.SourceFile.NativeIDType"))
                sourceFile.setNativeIDType((SourceFile::NativeIDType) findIndex_(
                  value, SourceFile::NamesOfNativeIDType, SourceFile::SIZE_OF_NATIVEIDTYPE));
            }
            else if(command.hasPrefix("SpectrumSettings.Precursors."))
            {
              std::vector<Precursor>& precursorsList = spectrum.getPrecursors();
              
              if(command.hasPrefix("SpectrumSettings.Precursors.Size"))
                resizeList_(precursorsList, value.toInt());
              else if(command.hasPrefix("SpectrumSettings.Precursors.Index"))
                setIndex_(indexPrecursor, value.toInt()-1, precursorsList);
              else /*if(command.hasPrefix("SpectrumSettings.Precursors.ActivationMethod"))
                precursorsList[indexPrecursor].setActivationMethod(findIndex_(value, Precursor::NamesOfActivationMethod, Precursor::SIZE_OF_ACTIVATIONMETHOD));
              else*/ if(command.hasPrefix("SpectrumSettings.Precursors.ActivationEnergy"))
                precursorsList[indexPrecursor].setActivationEnergy(value.toDouble());
              else if(command.hasPrefix("SpectrumSettings.Precursors.IsolationWindowLowerOffset"))
                precursorsList[indexPrecursor].setIsolationWindowLowerOffset(value.toDouble());
              else if(command.hasPrefix("SpectrumSettings.Precursors.IsolationWindowUpperOffset"))
                precursorsList[indexPrecursor].setIsolationWindowUpperOffset(value.toDouble());
              else if(command.hasPrefix("SpectrumSettings.Precursors.Charge"))
                precursorsList[indexPrecursor].setCharge(value.toInt());
              else if(command.hasPrefix("SpectrumSettings.Precursors.PossibleChargeStates"))
              {
                std::vector<Int> intVector;
                stringToIntVector_(value, intVector);
               	precursorsList[indexPrecursor].setPossibleChargeStates(intVector);
              }
            }    
            else if(command.hasPrefix("SpectrumSettings.Products."))
            {
              std::vector<Product>& productsList = spectrum.getProducts();
              
              if(command.hasPrefix("SpectrumSettings.Products.Size"))
                resizeList_(productsList, value.toInt());
              else if(command.hasPrefix("SpectrumSettings.Products.Index"))
                setIndex_(indexProduct, value.toInt()-1, productsList);
              else if(command.hasPrefix("SpectrumSettings.Products.MZ"))
                productsList[indexProduct].setMZ(value.toDouble());
              else if(command.hasPrefix("SpectrumSettings.Products.IsolationWindowLowerOffset"))
                productsList[indexProduct].setIsolationWindowLowerOffset(value.toDouble());
              else if(command.hasPrefix("SpectrumSettings.Products.IsolationWindowUpperOffset"))
                productsList[indexProduct].setIsolationWindowUpperOffset(value.toDouble());
            }
            else if(command.hasPrefix("SpectrumSettings.DataProcessing."))
            {
              std::vector<DataProcessing>& dataProcessingList = spectrum.getDataProcessing();

              if(command.hasPrefix("SpectrumSettings.DataProcessing.Size"))
                resizeList_(dataProcessingList, value.toInt());
              else if(command.hasPrefix("SpectrumSettings.DataProcessing.Index"))
                setIndex_(indexDataProcessing, value.toInt()-1, dataProcessingList);                
              else if(command.hasPrefix("SpectrumSettings.DataProcessing.Software."))
              {
                Software& software = dataProcessingList[indexDataProcessing].getSoftware();

                if(command.hasPrefix("SpectrumSettings.DataProcessing.Software.Name"))
                 software.setName(value);
                else if(command.hasPrefix("SpectrumSettings.DataProcessing.Software.Version"))
                 software.setVersion(value);
              }                  
              else if(command.hasPrefix("SpectrumSettings.DataProcessing.ProcessingActions"))
              {
                std::vector<String> processingActionsStrings;
                std::set<DataProcessing::ProcessingAction>& processingActionsList = dataProcessingList[indexDataProcessing].getProcessingActions();

                value.split(',', processingActionsStrings);
                for(Size index=0; index<processingActionsStrings.size(); ++index)
                  processingActionsList.insert((DataProcessing::ProcessingAction) findIndex_(
                    processingActionsStrings[index], DataProcessing::NamesOfProcessingAction, DataProcessing::SIZE_OF_PROCESSINGACTION));
              }
              else if(command.hasPrefix("SpectrumSettings.DataProcessing.CompletionTime"))
                dataProcessingList[indexDataProcessing].setCompletionTime(DateTime(value));
            }              
            else if(command.prefix(26) == "SpectrumSettings.PeptideIdentification.")
            {
              std::vector<PeptideIdentification>& peptideIdentificationsList = spectrum.getPeptideIdentifications();
              
              if(command.hasPrefix("SpectrumSettings.PeptideIdentification.Size"))
                resizeList_(peptideIdentificationsList, value.toInt());
              else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.Index"))
                setIndex_(indexPeptideIdentification, value.toInt()-1, peptideIdentificationsList);
              else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit."))
              {
                std::vector<PeptideHit>& peptideHitList = peptideIdentificationsList[indexPeptideIdentification].getHits();
                
                if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.Size"))
                  resizeList_(peptideHitList, value.toInt());
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.Index"))
                  setIndex_(indexPeptideHit, value.toInt()-1, peptideHitList);
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.ProteinAccessions"))
                {
                  std::vector<String> accessions;
                  value.split(',', accessions);
                  peptideHitList[indexPeptideHit].setProteinAccessions(accessions);
                }  
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.Score"))
                  peptideHitList[indexPeptideHit].setScore(value.toDouble());
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.Rank"))
                  peptideHitList[indexPeptideHit].setRank((UInt) value.toInt());
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.AASequence.StringSequence"))
                {
                  aaSequence_.setStringSequence(value);
                  peptideHitList[indexPeptideHit].setSequence(aaSequence_);
                }
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.AASequence.Modification"))
                {
                  std::vector<String> modification;
                  value.split(',', modification);
                  if(modification.size() == 2)
                  {
                    aaSequence_.setModification((Size) modification[0], modification[1]);
                    peptideHitList[indexPeptideHit].setSequence(aaSequence_);
                  }
                }
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.AASequence.NTerminalModification"))
                {
                  aaSequence_.setNTerminalModification(value);
                  peptideHitList[indexPeptideHit].setSequence(aaSequence_);
                }
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.AASequence.CTerminalModification"))
                {
                  aaSequence_.setCTerminalModification(value);
                  peptideHitList[indexPeptideHit].setSequence(aaSequence_);
                }
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.Charge"))
                  peptideHitList[indexPeptideHit].setCharge(value.toInt());                                             
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.AABefore"))
                  peptideHitList[indexPeptideHit].setAABefore(value[0]);  
                else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.PeptideHit.AAAfter"))
                  peptideHitList[indexPeptideHit].setAAAfter(value[0]);
              }
              else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.SignificanceThreshold"))
                peptideIdentificationsList[indexPeptideIdentification].setSignificanceThreshold(value.toDouble());
              else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.ScoreType"))
                peptideIdentificationsList[indexPeptideIdentification].setScoreType(value);                  
              else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.HigherScoreBetter"))
                peptideIdentificationsList[indexPeptideIdentification].setHigherScoreBetter(value=="TRUE" ? true : false);                       
              else if(command.hasPrefix("SpectrumSettings.PeptideIdentification.Identifier"))
                peptideIdentificationsList[indexPeptideIdentification].setIdentifier(value);
            }                                                              
          }
        } 
			  catch(const std::exception& e)
			  {
cout << "excep: " << e.what() << endl;					  
				  throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"");       
        }
      }

      template <class PeakType>
      void importSpectrumSettingsFromFile(const String& filename, MSSpectrum<PeakType>& spectrum) const
      {					        
        std::ifstream is(filename.c_str());
        if (!is)
        {
          throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
        }
                
        //temporary variables
        String line;
        
        //read lines
        while (getline(is, line, '\n'))
        {
          importSpectrumSettings(filename, spectrum);
        }
        
        is.close();
      }
              
      template <class PeakType>
      void importExperimentalSettings(const String& line, MSExperiment<PeakType>& exp) const
      {				        
        //temporary variables
        SourceFile::ChecksumType checksum_type_ = SourceFile::UNKNOWN_CHECKSUM;
        String checksum_ = "";
        Size indexProteinIdentification = 0;
        Size indexProteinHit = 0;        
        Size indexTreatment = 0;
        Size indexSourceFiles = 0;
        Size indexContacts = 0;
        Size indexDataProcessing = 0;
        Size indexIonSource = 0;
        Size indexMassAnalyzer = 0;
        Size indexIonDetector = 0;
        SampleTreatment* treatment_ = NULL;
        String command = "";
        String value = "";
        
        try
        {        
          if( line.empty() )
            return;
        
          {
            std::vector<String> strings;
            if(!line.split('=', strings))
              throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" );
              
            if( strings.size() != 2 )
            {
cout << "BBB: " << strings.size() << endl;
              throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" );
            }
              
            command = strings[0].trim();
            value = strings[1].trim();
          }

          if(command.hasPrefix("#"))
            return;             
          
          if(command.hasPrefix("ExperimentalSettings."))
          {
            if(command.hasPrefix("ExperimentalSettings.Sample."))
            {
              Sample& sample = exp.getSample();
              
              if(command.hasPrefix("ExperimentalSettings.Sample.Name"))
                sample.setName(value);
              else if(command.hasPrefix("ExperimentalSettings.Sample.Organism"))
                sample.setOrganism(value);
              else if(command.hasPrefix("ExperimentalSettings.Sample.Number"))
                sample.setNumber(value);
              else if(command.hasPrefix("ExperimentalSettings.Sample.Comment"))
                sample.setComment(value);
              else if(command.hasPrefix("ExperimentalSettings.Sample.State"))
                sample.setState((Sample::SampleState) findIndex_(value, Sample::NamesOfSampleState, Sample::SIZE_OF_SAMPLESTATE ));
              else if(command.hasPrefix("ExperimentalSettings.Sample.Mass"))
                sample.setMass(value.toDouble());
              else if(command.hasPrefix("ExperimentalSettings.Sample.Volume"))
                sample.setVolume(value.toDouble());
              else if(command.hasPrefix("ExperimentalSettings.Sample.Concentration"))
                sample.setConcentration(value.toDouble());
              else if(command.hasPrefix("ExperimentalSettings.Sample.Subsample."))
              {
                // std::vector< Sample > & 	getSubsamples ()
              }
              else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment."))
              {
                if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.Type"))
                {
                  if(value == "Digestion")
                    treatment_ = (SampleTreatment*) new Digestion();
                  else if(value == "Modification")
                    treatment_ = (SampleTreatment*) new Modification();
                  else if(value == "Tagging")
                    treatment_ = (SampleTreatment*) new Tagging();
                  sample.addTreatment(*treatment_);
                }
                else if(treatment_ != NULL)
                {
                  if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.Comment"))
                    treatment_->setComment(value);
                  else if(treatment_->getType() == "Digestion")
                  {
                    Digestion* digestion = (Digestion*) treatment_;
                    
                    if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.Enzyme"))
                      digestion->setEnzyme(value);
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.DigestionTime"))
                      digestion->setDigestionTime(value.toDouble());
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.Temperature"))
                      digestion->setTemperature(value.toDouble());
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.pH"))
                      digestion->setPh(value.toDouble());                                                        
                  }
                  else if(treatment_->getType() == "Modification" || treatment_->getType() == "Tagging")
                  {
                    Modification* modification = (Modification*) treatment_;
                  
                    if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.ReagentName"))
                      modification->setReagentName(value);
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.Mass"))
                      modification->setMass(value.toDouble());
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.SpecificityType"))
                      modification->setSpecificityType((Modification::SpecificityType) findIndex_(
                        value, Modification::NamesOfSpecificityType, Modification::SIZE_OF_SPECIFICITYTYPE));
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.AffectedAminoAcids"))
                      modification->setAffectedAminoAcids(value);                         
                  }
                  else if(treatment_->getType() == "Tagging")
                  {
                    Tagging* tagging = (Tagging*) treatment_;
                  
                    if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.MassShift"))
                      tagging->setMassShift(value.toDouble());
                    else if(command.hasPrefix("ExperimentalSettings.Sample.Treatment.IsotopeVariant"))
                      tagging->setVariant((Tagging::IsotopeVariant) findIndex_(
                        value, Tagging::NamesOfIsotopeVariant, Tagging::SIZE_OF_ISOTOPEVARIANT));                                                      
                  }     
                }                 
              }
            }
            else if(command.hasPrefix("ExperimentalSettings.SourceFiles."))
            {
              std::vector<SourceFile>& sourceFilesList = exp.getSourceFiles();

              if(command.hasPrefix("ExperimentalSettings.SourceFiles.Size"))
                resizeList_(sourceFilesList, value.toInt());
              else if(command.hasPrefix("ExperimentalSettings.SourceFiles.Index"))
                setIndex_(indexSourceFiles, value.toInt()-1, sourceFilesList);
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.NameOfFile"))
                sourceFilesList[indexSourceFiles].setNameOfFile(value);
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.PathToFile"))
                sourceFilesList[indexSourceFiles].setPathToFile(value);
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.FileSize"))
                sourceFilesList[indexSourceFiles].setFileSize(value.toFloat());
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.FileType"))
                sourceFilesList[indexSourceFiles].setFileType(value);
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.ChecksumType"))
              {
                checksum_type_ = findIndex_(value, SourceFile::NamesOfChecksumType, SourceFile::SIZE_OF_CHECKSUMTYPE);
                sourceFilesList[indexSourceFiles].setChecksum(checksum_, checksum_type_);
              }
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.Checksum"))
              {
                checksum_ = value;
                sourceFilesList[indexSourceFiles].setChecksum(checksum_, checksum_type_);
              }
              else if(command.hasPrefix("ExperimentalSettings.SourceFile.NativeIDType"))
                sourceFilesList[indexSourceFiles].setNativeIDType((SourceFile::NativeIDType) findIndex_(
                  value, SourceFile::NamesOfNativeIDType, SourceFile::SIZE_OF_NATIVEIDTYPE));
            } 
            else if(command.hasPrefix("ExperimentalSettings.Contacts."))
            {
              std::vector<ContactPerson>& contactsList = exp.getContacts();

              if(command.hasPrefix("ExperimentalSettings.Contacts.Size"))
                resizeList_(contactsList, value.toInt());
              else if(command.hasPrefix("ExperimentalSettings.Contacts.Index"))
                setIndex_(indexContacts, value.toInt()-1, contactsList);                
              else if(command.hasPrefix("ExperimentalSettings.Contacts.FirstName"))
                contactsList[indexContacts].setFirstName(value);
              else if(command.hasPrefix("ExperimentalSettings.Contacts.LastName"))
                contactsList[indexContacts].setLastName(value);
              else if(command.hasPrefix("ExperimentalSettings.Contacts.Institution"))
                contactsList[indexContacts].setInstitution(value);
              else if(command.hasPrefix("ExperimentalSettings.Contacts.Email"))
                contactsList[indexContacts].setEmail(value);
              else if(command.hasPrefix("ExperimentalSettings.Contacts.URL"))
                contactsList[indexContacts].setURL(value);
              else if(command.hasPrefix("ExperimentalSettings.Contacts.Address"))
                contactsList[indexContacts].setAddress(value);
              else if(command.hasPrefix("ExperimentalSettings.Contacts.ContactInfo"))
                contactsList[indexContacts].setContactInfo(value);
            }    
            else if(command.hasPrefix("ExperimentalSettings.Instrument."))
            {
              Instrument& instrument = exp.getInstrument();
              
              if(command.hasPrefix("ExperimentalSettings.Instrument.Name"))
                instrument.setName(value);
              else if(command.hasPrefix("ExperimentalSettings.Instrument.Vendor"))
                instrument.setVendor(value);
              else if(command.hasPrefix("ExperimentalSettings.Instrument.Model"))
                instrument.setModel(value);
              else if(command.hasPrefix("ExperimentalSettings.Instrument.Customizations"))
                instrument.setCustomizations(value);
              else if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources."))
              {
                std::vector<IonSource>& ionSourceList = instrument.getIonSources();

                if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources.Size"))
                  resizeList_(ionSourceList, value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources.Index"))
                  setIndex_(indexIonSource, value.toInt()-1, ionSourceList);                
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources.InletType"))
                  ionSourceList[indexIonSource].setInletType((IonSource::InletType) findIndex_(
                    value, IonSource::NamesOfInletType, IonSource::SIZE_OF_INLETTYPE));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources.IonizationMethod"))
                  ionSourceList[indexIonSource].setIonizationMethod((IonSource::IonizationMethod) findIndex_(
                    value, IonSource::NamesOfIonizationMethod, IonSource::SIZE_OF_IONIZATIONMETHOD));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources.Polarity"))
                  ionSourceList[indexIonSource].setPolarity((IonSource::Polarity) findIndex_(
                    value, IonSource::NamesOfPolarity, IonSource::SIZE_OF_POLARITY));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonSources.Order"))
                  ionSourceList[indexIonSource].setOrder(value.toInt());               
              }
              else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers."))
              {
                std::vector<MassAnalyzer>& massAnalyzerList = instrument.getMassAnalyzers();

                if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.Size"))
                  resizeList_(massAnalyzerList, value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.Index"))
                  setIndex_(indexMassAnalyzer, value.toInt()-1, massAnalyzerList);                
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.Type"))
                  massAnalyzerList[indexMassAnalyzer].setType((MassAnalyzer::AnalyzerType) findIndex_(
                    value, MassAnalyzer::NamesOfAnalyzerType, MassAnalyzer::SIZE_OF_ANALYZERTYPE));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ResolutionMethod"))
                  massAnalyzerList[indexMassAnalyzer].setResolutionMethod((MassAnalyzer::ResolutionMethod) findIndex_(
                    value, MassAnalyzer::NamesOfResolutionMethod, MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ResolutionType"))
                  massAnalyzerList[indexMassAnalyzer].setResolutionType((MassAnalyzer::ResolutionType) findIndex_(
                    value, MassAnalyzer::NamesOfResolutionType, MassAnalyzer::SIZE_OF_RESOLUTIONTYPE));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ScanDirection"))
                  massAnalyzerList[indexMassAnalyzer].setScanDirection((MassAnalyzer::ScanDirection) findIndex_(
                    value, MassAnalyzer::NamesOfScanDirection, MassAnalyzer::SIZE_OF_SCANDIRECTION));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ScanLaw"))
                  massAnalyzerList[indexMassAnalyzer].setScanLaw((MassAnalyzer::ScanLaw) findIndex_(
                    value, MassAnalyzer::NamesOfScanLaw, MassAnalyzer::SIZE_OF_SCANLAW));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ReflectronState"))
                  massAnalyzerList[indexMassAnalyzer].setReflectronState((MassAnalyzer::ReflectronState) findIndex_(
                    value, MassAnalyzer::NamesOfReflectronState, MassAnalyzer::SIZE_OF_REFLECTRONSTATE));
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.Resolution"))
                  massAnalyzerList[indexMassAnalyzer].setResolution(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.Accuracy"))
                  massAnalyzerList[indexMassAnalyzer].setAccuracy(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ScanRate"))
                  massAnalyzerList[indexMassAnalyzer].setScanRate(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.ScanTime"))
                  massAnalyzerList[indexMassAnalyzer].setScanTime(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.TOFTotalPathLength"))
                  massAnalyzerList[indexMassAnalyzer].setTOFTotalPathLength(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.IsolationWidth"))
                  massAnalyzerList[indexMassAnalyzer].setIsolationWidth(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.FinalMSExponent"))
                  massAnalyzerList[indexMassAnalyzer].setFinalMSExponent(value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.MagneticFieldStrength"))
                  massAnalyzerList[indexMassAnalyzer].setMagneticFieldStrength(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.MassAnalyzers.Order"))
                  massAnalyzerList[indexMassAnalyzer].setOrder(value.toInt());                
              }
              else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetectors."))
              {
                std::vector<IonDetector>& IonDetectorList = instrument.getIonDetectors();

                if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.Size"))
                  resizeList_(IonDetectorList, value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.Index"))
                  setIndex_(indexIonDetector, value.toInt()-1, IonDetectorList);                
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.Type"))
                  IonDetectorList[indexIonDetector].setType((IonDetector::Type) findIndex_(
                    value, IonDetector::NamesOfType, IonDetector::SIZE_OF_TYPE));  
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.AcquisitionMode"))
                  IonDetectorList[indexIonDetector].setAcquisitionMode((IonDetector::AcquisitionMode) findIndex_(
                    value, IonDetector::NamesOfAcquisitionMode, IonDetector::SIZE_OF_ACQUISITIONMODE));  
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.Resolution"))
                  IonDetectorList[indexIonDetector].setResolution(value.toDouble());  
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.ADCSamplingFrequency"))
                  IonDetectorList[indexIonDetector].setADCSamplingFrequency(value.toDouble());  
                else if(command.hasPrefix("ExperimentalSettings.Instrument.IonDetector.Order"))
                  IonDetectorList[indexIonDetector].setOrder(value.toInt());              
              }
              else if(command.hasPrefix("ExperimentalSettings.Instrument.Software."))
              {
                Software& software = instrument.getSoftware();

                if(command.hasPrefix("ExperimentalSettings.Instrument.Software.Name"))
                 software.setName(value);
                else if(command.hasPrefix("ExperimentalSettings.Instrument.Software.Version"))
                 software.setVersion(value);
              }
              else if(command.hasPrefix("ExperimentalSettings.Instrument.IonOptics"))
                instrument.setIonOptics((Instrument::IonOpticsType) findIndex_(
                  value, Instrument::NamesOfIonOpticsType, Instrument::SIZE_OF_IONOPTICSTYPE));
            }  
            else if(command.hasPrefix("ExperimentalSettings.HPLC."))
            {
              HPLC& hplc = exp.getHPLC();
              
              if(command.hasPrefix("ExperimentalSettings.HPLC.Instrument"))
                hplc.setInstrument(value);
              else if(command.hasPrefix("ExperimentalSettings.HPLC.Column"))
                hplc.setColumn(value);
              else if(command.hasPrefix("ExperimentalSettings.HPLC.Temperature"))
                hplc.setTemperature(value.toInt());
              else if(command.hasPrefix("ExperimentalSettings.HPLC.Pressure"))
                hplc.setPressure((UInt) value.toInt());
              else if(command.hasPrefix("ExperimentalSettings.HPLC.Flux"))
                hplc.setFlux((UInt) value.toInt());
              else if(command.hasPrefix("ExperimentalSettings.HPLC.Comment"))
                hplc.setComment(value);
              else if(command.hasPrefix("ExperimentalSettings.HPLC.Gradient."))
              {
                Gradient& gradient = hplc.getGradient();
                
                if(command.hasPrefix("ExperimentalSettings.HPLC.Gradient.Eluent"))
                  gradient.addEluent(value);
                else if(command.hasPrefix("ExperimentalSettings.HPLC.Gradient.Timepoint"))
                  gradient.addTimepoint(value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.HPLC.Gradient.Percentage"))
                {
                  std::vector<String> percentageInfo;
                  value.split(',', percentageInfo);
                  if(percentageInfo.size() == 3)
                    gradient.setPercentage(percentageInfo[0], percentageInfo[1].toInt(), (UInt) percentageInfo[2].toInt());
                }
              }
            }    
            else if(command.hasPrefix("ExperimentalSettings.DateTime"))
              exp.setDateTime(DateTime(value));
            else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification."))
            {
              std::vector<ProteinIdentification>& ProteinIdentificationsList = exp.getProteinIdentifications();
              
              if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.Size"))
                resizeList_(ProteinIdentificationsList, value.toInt());
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.Index"))
                setIndex_(indexProteinIdentification, value.toInt()-1, ProteinIdentificationsList);
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit."))
              {
                std::vector<ProteinHit>& ProteinHitList = ProteinIdentificationsList[indexProteinIdentification].getHits();
                
                if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit.Size"))
                  resizeList_(ProteinHitList, value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit.Index"))
                  setIndex_(indexProteinHit, value.toInt()-1, ProteinHitList);
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit.Score"))
                  ProteinHitList[indexProteinHit].setScore(value.toDouble());
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit.Rank"))
                  ProteinHitList[indexProteinHit].setRank((UInt) value.toInt());
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit.Accession"))
                  ProteinHitList[indexProteinHit].setAccession(value);                  
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ProteinHit.Sequence"))
                  ProteinHitList[indexProteinHit].setSequence(value);
              }
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SignificanceThreshold"))
                ProteinIdentificationsList[indexProteinIdentification].setSignificanceThreshold(value.toDouble());
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.ScoreType"))
                ProteinIdentificationsList[indexProteinIdentification].setScoreType(value);                  
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.HigherScoreBetter"))
                ProteinIdentificationsList[indexProteinIdentification].setHigherScoreBetter(value=="TRUE" ? true : false);                       
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.DateTime"))
                ProteinIdentificationsList[indexProteinIdentification].setDateTime(DateTime(value));
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchEngine"))
                ProteinIdentificationsList[indexProteinIdentification].setSearchEngine(value);
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchEngineVersion"))
                ProteinIdentificationsList[indexProteinIdentification].setSearchEngineVersion(value);
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters."))
              {
                ProteinIdentification::SearchParameters& searchParameters = ProteinIdentificationsList[indexProteinIdentification].getSearchParameters();
                
                if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.Database"))
                  searchParameters.db = value;
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.DatabaseVersion"))
                  searchParameters.db_version = value;
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.Taxonomy"))
                  searchParameters.taxonomy = value;
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.Charges"))
                  searchParameters.charges = value;
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.MassType"))
                  searchParameters.mass_type = (ProteinIdentification::PeakMassType) findIndex_(
                    value, ProteinIdentification::NamesOfPeakMassType, ProteinIdentification::SIZE_OF_PEAKMASSTYPE);
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.FixedModifications"))
                  value.split(',', searchParameters.fixed_modifications);
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.VariableModifications"))
                  value.split(',', searchParameters.variable_modifications);
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.Enzyme"))
                  searchParameters.enzyme = (ProteinIdentification::DigestionEnzyme) findIndex_(
                    value, ProteinIdentification::NamesOfDigestionEnzyme, ProteinIdentification::SIZE_OF_DIGESTIONENZYME);
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.MissedCleavages"))
                  searchParameters.missed_cleavages = (UInt) value.toInt();
                else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.SearchParameters.PeakMassTolerance"))
                  searchParameters.peak_mass_tolerance = value.toDouble();
              } 
              else if(command.hasPrefix("ExperimentalSettings.ProteinIdentification.Identifier"))
                ProteinIdentificationsList[indexProteinIdentification].setIdentifier(value);                                  
            }
          }
        }
			  catch(const std::exception& e)
			  {
cout << "excep: " << e.what() << endl;					  
				  throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" );
			  }        
      }

      template <class PeakType>
      void importExperimentalSettingsFromFile(const String& filename, MSExperiment<PeakType>& exp) const
      {				        
        std::ifstream is(filename.c_str());
        if (!is)
        {
          throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
        }
        
        //temporary variables
        String line;
        
        //read lines
        while (getline(is, line, '\n'))
        {
          importExperimentalSettingsFromFile(filename, exp);
        }
        is.close();
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_SETTINGSPARSER_H

