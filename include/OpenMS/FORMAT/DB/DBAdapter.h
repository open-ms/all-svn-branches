// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DB_DBADAPTER_H
#define OPENMS_FORMAT_DB_DBADAPTER_H

//#include <OpenMS/config.h>

//OpenMS includes
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

//QT includes
#include <qsqlquery.h>
#include <qdatetime.h>
#include <QVariant>

//std and STL includes
#include <string>
#include <map>

namespace OpenMS
{	
	class DBConnection;
	
  /** 
  	@brief A class for accessing and storing data in a SQL database
  	
    It can be used to create objects from the DB or store them in the DB.
    
    @ingroup DatabaseIO
  */
 
  class DBAdapter
  {
    public:
			/// Constructor
      DBAdapter(DBConnection& db_con);

      /// Destructor
      ~DBAdapter();

			/// Stores a MSExperiment
			template <class ExperimentType>
			void storeExperiment(ExperimentType& exp);

			/// Reads a MSExperiment according to a set of options
			template <class ExperimentType>
			void loadExperiment(UID id, ExperimentType& exp, PeakFileOptions options = PeakFileOptions());

			template <class SpectrumType>
			/// Reads a MSSpectrum according to a set of options
			void loadSpectrum(UID id, SpectrumType& spec, PeakFileOptions options = PeakFileOptions());

			/**
				@brief Returns true if the DB is up-to-date (Checks the version in ADMIN_Version table).
				
				@param warning if this is set, a warning is issued to stderr if the db is not up-to-date. 
			*/
			bool checkDBVersion(bool warning);
			
			/// Deletes all tables in the database and creates a new OpenMS database
			void createDB();
			
    private:
    	/// Reference to the DB connection handed over in the constructor
    	DBConnection& db_con_;
			
			/// Not implemented, thus private
			DBAdapter();
			
			/**
				@brief Stores, updates or deletes MetaInfo data
				
				@return the id of the new META_MetaInfo table row
			*/
			UID storeMetaInfo_(const String& parent_table, UID parent_id, const MetaInfoInterface& info);	
			
			/**
				@brief Loads MetaInfo data from database
				
			*/
			void loadMetaInfo_(UID id, MetaInfoInterface& info);
			
			/**
				@brief Conditionally deletes MetaInfo data from database
				
			*/
			void deleteMetaInfo_(const String& parent_table, const String& condition);
			
			/**
				@brief Stores, updates or deletes file information
				
				@return the id of the new META_File table row
			*/
			UID storeFile_(const String& parent_table, UID parent_id, const SourceFile& file);
			
			/**
				@brief Loads file information
				
			*/
			void loadFile_(UID id, SourceFile& file);
  
  			/**
				@brief Stores, updates or deletes sample information
				
				@return the id of the new META_Sample table row
			*/
			UID storeSample_(const Sample& sample, UID exp_id, UID parent_id);
			
			/**
				@brief Loads sample information
				
			*/
			void loadSample_(UID id, Sample& sample);
};


//------------------------------------------- IMPLEMENTATION OF TEMPLATE METHODS ----------------------------------

	template <class ExperimentType>
	void DBAdapter::storeExperiment(ExperimentType& exp)
	{
		std::stringstream query; // query to build
		String end;              // end of the query that is added afer all fields
		String tmp;              // temporary data
		bool new_entry(false);          // stores if the current object is already in the DB
		QSqlQuery result;        // place to store the query results in
		int parent_id(-1);           // stores parent_id of meta information
		UID acquisition_info_id(0); // stores id of acquisition_info
		UID meta_id(0);             // stores MetaInfo id of meta information that was just stored
		UID meta_parent_id(0);      // stores parent ID of MetaInfo that will be stored

		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION --------------------------------------- 
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;
		
		//----------------------------------------------------------------------------------------
		//------------------------------- EXPERIMENTAL SETTINGS ---------------------------------- 
		//----------------------------------------------------------------------------------------
		
		new_entry = (exp.getPersistenceId()==0);
		if (new_entry)
		{
			query << "INSERT INTO META_MSExperiment SET ";
			end = "";
		}
		else
		{
			query << "UPDATE META_MSExperiment SET ";
			end = " WHERE id='" + String(exp.getPersistenceId()) + "'";
		}
		//type
		query << "Type=" << (1u+exp.getType());
		//date
		exp.getDate().get(tmp);
		query << ",Date='" << tmp << "'";
		//description
		query << ",Description='" << exp.getComment() << "'";
		
		query << end;
		db_con_.executeQuery(query.str(),result);
		if (new_entry)
		{
			exp.setPersistenceId(db_con_.getAutoId());
		}
		
		storeMetaInfo_("META_MSExperiment", exp.getPersistenceId(), exp);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- SAMPLE ------------------------------------------ 
		//----------------------------------------------------------------------------------------
		
		query.str("");
		deleteMetaInfo_("META_Sample", "fid_MSExperiment=" + String(exp.getPersistenceId()));
		// this also deletes all references in META_SampleTreatment, META_Digestion and META_Modification by constraint
		query << "DELETE FROM META_Sample WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
		storeSample_(exp.getSample(), exp.getPersistenceId(), 0);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- CONTACTPERSON ----------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const std::vector<ContactPerson>& contacts = exp.getContacts();
		
		query.str("");
		deleteMetaInfo_("META_ContactPerson", "fid_MSExperiment=" + String(exp.getPersistenceId()));
		query << "DELETE FROM META_ContactPerson WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
		db_con_.executeQuery(query.str(), result);
			
		for (std::vector<ContactPerson>::const_iterator contact_it = contacts.begin(); contact_it != contacts.end(); contact_it++)
		{
			query.str("");
			query << "INSERT INTO META_ContactPerson SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "'";
			query << ",PreName='" << contact_it->getFirstName() << "'";
			query << ",LastName='" << contact_it->getLastName() << "'";
			query << ",Affiliation='" << contact_it->getInstitution() << "'";
			query << ",Email='" << contact_it->getEmail() << "'";
			query << ",Comment='" << contact_it->getContactInfo() << "'";
			
			db_con_.executeQuery(query.str(), result);
			parent_id = db_con_.getAutoId();
			
			storeMetaInfo_("META_ContactPerson", parent_id, *contact_it);
		}
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- HPLC -------------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const HPLC& hplc = exp.getHPLC();
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_HPLC SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_HPLC WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
			db_con_.executeQuery(query.str(),result);
			result.first();
			parent_id = result.value(0).asInt();
			
			query.str("");
			query << "UPDATE META_HPLC SET ";
			end = " WHERE fid_MSExperiment='" + String (exp.getPersistenceId()) + "'";
		}
		
		query << "InstrumentName='" << hplc.getInstrument() << "'";
		query << ",ColumnName='" << hplc.getColumn() << "'";
		query << ",Description='" << hplc.getComment() << "'";
		query << ",Flux=" << hplc.getFlux();
		query << ",Pressure=" << hplc.getPressure();
		query << ",Temperature=" << hplc.getTemperature();
		
		query << end;
		db_con_.executeQuery(query.str(),result);
		
		if (new_entry)
		{
			parent_id = db_con_.getAutoId();
		}
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- GRADIENT* --------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const Gradient& gradient = exp.getHPLC().getGradient();
		const std::vector <String>& eluents = gradient.getEluents();
		const std::vector <SignedInt>& time = gradient.getTimepoints();
		const std::vector< std::vector< UnsignedInt > >& percentages = gradient.getPercentages();
		std::stringstream query_eluents, query_time, query_percentages;
		UID eluents_id, time_id;
		
		// this also deletes all references in META_GradientPercentage by constraint
		query.str("");
		query << "DELETE FROM META_GradientEluent WHERE fid_HPLC=" << parent_id;
		db_con_.executeQuery(query.str(),result);
		query.str("");
		query << "DELETE FROM META_GradientTime WHERE fid_HPLC=" << parent_id;
		db_con_.executeQuery(query.str(),result);
		
		if (! eluents.empty())
		{
			query_eluents.str("");
			query_eluents << "INSERT INTO META_GradientEluent (fid_HPLC, Name) VALUES ";
			for (std::vector<String>::const_iterator eluents_it = eluents.begin(); eluents_it != eluents.end(); eluents_it++)
			{
				query_eluents << "(";
				query_eluents << parent_id;
				query_eluents << ",'" << *eluents_it << "'";
				query_eluents << "),";
			}
			db_con_.executeQuery(String(query_eluents.str()).substr(0,-1),result);
			eluents_id = db_con_.getAutoId();
		}
		
		if (! time.empty())
		{
			query_time.str("");
			query_time << "INSERT INTO META_GradientTime (fid_HPLC, Time) VALUES ";
			for (std::vector<SignedInt>::const_iterator time_it = time.begin(); time_it != time.end(); time_it++)
			{
				query_time << "(";
				query_time << parent_id;
				query_time << "," << *time_it;
				query_time << "),";
			}
			db_con_.executeQuery(String(query_time.str()).substr(0,-1),result);
			time_id = db_con_.getAutoId();
		}
		
		if (! percentages.empty() && ! eluents.empty() && ! time.empty())
		{
			query_percentages.str("");
			query_percentages << "INSERT INTO META_GradientPercentage (fid_GradientEluent, fid_GradientTime, Percentage) VALUES ";
			int i = 0;
			// iterate over eluents
			for (std::vector< std::vector< UnsignedInt> >::const_iterator percent_outer_it = percentages.begin(); percent_outer_it != percentages.end(); percent_outer_it++)
			{
				int j = 0;
				// iterate over timepoints
				for (std::vector< UnsignedInt>::const_iterator percent_inner_it = (*percent_outer_it).begin(); percent_inner_it != (*percent_outer_it).end(); percent_inner_it++)
				{
					query_percentages << "(";
					query_percentages << eluents_id + i;
					query_percentages << "," << time_id + j;
					query_percentages << "," << *percent_inner_it;
					query_percentages << "),";
					j++;
				}
				i++;
			}
			db_con_.executeQuery(String(query_percentages.str()).substr(0,-1),result);
		}
			
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- INSTRUMENT -------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const Instrument& instrument = exp.getInstrument();
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_MSInstrument SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_MSInstrument WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
			db_con_.executeQuery(query.str(),result);
			result.first();
			parent_id = result.value(0).asInt();
			
			query.str("");
			query << "UPDATE META_MSInstrument SET ";
			end = " WHERE fid_MSExperiment='" + String (exp.getPersistenceId()) + "'";
		}
		
		query << "Model='" << instrument.getModel() << "'";
		query << ",Vendor='" << instrument.getVendor() << "'";
		query << ",Description='" << instrument.getCustomizations() << "'";
		
		query << end;
		db_con_.executeQuery(query.str(),result);
		
		if (new_entry)
		{
			parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_MSInstrument", parent_id, instrument);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- IONDETECTOR ------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const IonDetector& detector = exp.getInstrument().getIonDetector();
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_IonDetector SET ";
			query << "fid_MSInstrument='" << parent_id << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_IonDetector WHERE fid_MSInstrument='" << parent_id << "'";
			db_con_.executeQuery(query.str(),result);
			result.first();
			meta_parent_id = result.value(0).asInt();
			
			query.str("");
			query << "UPDATE META_IonDetector SET ";
			end = " WHERE fid_MSInstrument='" + String (parent_id) + "'";
		}
		
		query << "AcquisitionMode=" << (1u+detector.getAcquisitionMode());
		query << ",Type=" << (1u+detector.getType());
		query << ",Resolution=" << detector.getResolution();
		query << ",ADCSamplingFrequency=" << detector.getADCSamplingFrequency();
		
		query << end;
		db_con_.executeQuery(query.str(),result);
		
		if (new_entry)
		{
			meta_parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_IonDetector", meta_parent_id, detector);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- IONSOURCE --------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const IonSource& source = exp.getInstrument().getIonSource();
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_IonSource SET ";
			query << "fid_MSInstrument='" << parent_id << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_IonSource WHERE fid_MSInstrument='" << parent_id << "'";
			db_con_.executeQuery(query.str(),result);
			
			query.str("");
			query << "UPDATE META_IonSource SET ";
			end = " WHERE fid_MSInstrument='" + String (parent_id) + "'";
			
			result.first();
			meta_parent_id = result.value(0).asInt();
		}
		
		query << "InletType=" << (1u+source.getInletType());
		query << ",IonizationMethod=" << (1u+source.getIonizationMethod());
		query << ",IonizationMode=" << (1u+source.getPolarity());
		
		query << end;
		db_con_.executeQuery(query.str(),result);
		
		if (new_entry)
		{
			meta_parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_IonSource", meta_parent_id, source);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- MASSANALYZER ------------------------------------ 
		//----------------------------------------------------------------------------------------
		
		const std::vector<MassAnalyzer>& analyzers = exp.getInstrument().getMassAnalyzers();
		query.str("");
		
		deleteMetaInfo_("META_MassAnalyzer", "fid_MSInstrument=" + String(parent_id));
		query << "DELETE FROM META_MassAnalyzer WHERE fid_MSInstrument='" << parent_id << "'";
		db_con_.executeQuery(query.str(),result);
			
		for (std::vector<MassAnalyzer>::const_iterator analyzer_it = analyzers.begin(); analyzer_it != analyzers.end(); analyzer_it++)
		{
			query.str("");
			query << "INSERT INTO META_MassAnalyzer SET ";
			query << "fid_MSInstrument='" << parent_id << "'";
			query << ",Accuracy=" << analyzer_it->getAccuracy();
			query << ",FinalMSExponent=" << analyzer_it->getFinalMSExponent();
			query << ",IsolationWidth=" << analyzer_it->getIsolationWidth();
			query << ",MagneticFieldStrength=" << analyzer_it->getMagneticFieldStrength();
			query << ",ReflectronState=" << (1u+analyzer_it->getReflectronState());
			query << ",Resolution=" << analyzer_it->getResolution();
			query << ",ResolutionMethod=" << (1u+analyzer_it->getResolutionMethod());
			query << ",ResolutionType=" << (1u+analyzer_it->getResolutionType());
			query << ",ScanDirection=" << (1u+analyzer_it->getScanDirection());
			query << ",ScanFunction=" << (1u+analyzer_it->getScanFunction());
			query << ",ScanLaw=" << (1u+analyzer_it->getScanLaw());
			query << ",ScanRate=" << analyzer_it->getScanRate();
			query << ",ScanTime=" << analyzer_it->getScanTime();
			query << ",TandemScanningMethod=" << (1u+analyzer_it->getTandemScanMethod());
			query << ",TOFPathLength=" << analyzer_it->getTOFTotalPathLength();
			query << ",Type=" << (1u+analyzer_it->getType());
			
			db_con_.executeQuery(query.str(),result);
			storeMetaInfo_("META_MassAnalyzer", db_con_.getAutoId(), *analyzer_it);
		}
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- SPECTRUM ---------------------------------------- 
		//----------------------------------------------------------------------------------------
		for (typename ExperimentType::Iterator exp_it = exp.begin(); exp_it != exp.end(); ++exp_it)
		{
			query.str("");
			new_entry = (exp_it->getPersistenceId()==0);
			if (new_entry)
			{
				query << "INSERT INTO DATA_Spectrum SET ";
				end = "";
			}
			else
			{
				query << "UPDATE DATA_Spectrum SET ";
				end = " WHERE id='" + String (exp_it->getPersistenceId()) + "'";
			}
			//FC (MSExperiment)
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "'";
			//type
			query << ",Type=" << (1u+exp_it->getType());
			//RT
			query << ",RetentionTime='" << exp_it->getRetentionTime() << "'";
			//MS-Level
			query << ",MSLevel='" << exp_it->getMSLevel() << "'";
			//Description
			query << ",Description='" << exp_it->getComment() << "'";
			
			//TODO: MassType (average/monoisotopic)
			//TODO: TIC
			
			
			query << end;
			db_con_.executeQuery(query.str(),result);
			if (new_entry)
			{
				exp_it->setPersistenceId(db_con_.getAutoId());
			}
			storeFile_("DATA_Spectrum", exp_it->getPersistenceId(), exp_it->getSourceFile());
			meta_id = storeMetaInfo_("DATA_Spectrum", exp_it->getPersistenceId(), *exp_it);
			
			
			//----------------------------------------------------------------------------------------
			//--------------------------------- METAINFODESCRIPTION ---------------------------------- 
			//----------------------------------------------------------------------------------------
			
			const std::map<String,MetaInfoDescription>& desc = exp_it->getMetaInfoDescriptions();
			
			for (std::map<String,MetaInfoDescription>::const_iterator desc_it = desc.begin(); desc_it != desc.end(); ++desc_it)
			{
				// first check if there is already an entry in META_MetaInfoDescription for this spectrum and this name
				// We cannot simply delete all entries for the spectrum because this might leave unreferenced META_TNVs 
				// and META_MetaInfos in the database.
				query.str("");
				query << "SELECT id FROM META_MetaInfoDescription WHERE fid_Spectrum=";
				query << exp_it->getPersistenceId();
				query << " AND Name='" << desc_it->first << "'";
				
				db_con_.executeQuery(query.str(), result);
				
				query.str("");
				
				if (result.size() > 0)
				{
					parent_id = result.value(0).asInt();
					new_entry = false;
					query << "UPDATE META_MetaInfoDescription SET ";
					end  = " WHERE fid_Spectrum=" + String(exp_it->getPersistenceId());
					end += " AND Name='" + String(desc_it->first) + "'";
				}
				else
				{
					new_entry = true;
					query << "INSERT INTO META_MetaInfoDescription SET ";
					query << "fid_Spectrum='" << exp_it->getPersistenceId() << "', ";
					query << "Name='" << desc_it->first << "',";
					end = "";
				}
				
				query << "Description='" << desc_it->second.getComment() << "'";
				query << end;
				
				db_con_.executeQuery(query.str(), result);
				if (new_entry)
				{
					parent_id = db_con_.getAutoId();
				}
				
				storeFile_("META_MetaInfoDescription", parent_id, desc_it->second.getSourceFile());
				storeMetaInfo_("META_MetaInfoDescription", parent_id, desc_it->second);
			}
			

			//----------------------------------------------------------------------------------------
			//-------------------------------------- PRECURSOR --------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			if (exp_it->getMSLevel()>1)
			{
				query.str("");
				if (new_entry)
				{
					query << "INSERT INTO DATA_Precursor SET fid_Spectrum='" << exp_it->getPersistenceId() << "',";
					end = "";
				}
				else
				{
					// first determine parent_id for meta information.
					// TODO: determination of parent_id can be optimized as soon as precursor becomes a PersistentObject
					query << "SELECT id FROM DATA_Precursor WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
					db_con_.executeQuery(query.str(),result);
					result.first();
					parent_id = result.value(0).asInt();
					
					query.str("");
					query << "UPDATE DATA_Precursor SET ";
					end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
				}
				//Intensity
				query << "Intensity='" << exp_it->getPrecursorPeak().getIntensity() << "'";
				//mz
				query << ",mz='" << exp_it->getPrecursorPeak().getPosition()[0] << "'";
				//charge
				query << ",Charge='" << exp_it->getPrecursorPeak().getCharge() << "'";
				//activation method
				query << ",ActivationMethod=" << (1u+exp_it->getPrecursor().getActivationMethod());		
				//activation energy unit
				query << ",ActivationEnergyUnit=" << (1u+exp_it->getPrecursor().getActivationEnergyUnit());		
				//activation energy
				query << ",ActivationEnergy='" << exp_it->getPrecursor().getActivationEnergy() << "'";
				//window size
				query << ",WindowSize='" << exp_it->getPrecursor().getWindowSize() << "'";
				
				query << end;
				db_con_.executeQuery(query.str(),result);
				if (new_entry) parent_id = db_con_.getAutoId();
				storeMetaInfo_("DATA_Precursor",parent_id, exp_it->getPrecursor());
				//TODO store persistence ID => Precusor class a persistent object
			}
			
			//----------------------------------------------------------------------------------------
			//---------------------------------------- PEAKS ----------------------------------------- 
			//----------------------------------------------------------------------------------------
			query.str("");
			deleteMetaInfo_("DATA_Peak", "fid_Spectrum=" + String(exp_it->getPersistenceId()));
			db_con_.executeQuery("DELETE FROM DATA_Peak WHERE fid_Spectrum=" + String(exp_it->getPersistenceId()), result);
			query << "INSERT INTO DATA_Peak (fid_Spectrum,Intensity,mz) VALUES ";
			tmp = "(" + String(exp_it->getPersistenceId()) + ",'";
			for (typename ExperimentType::SpectrumType::Iterator spec_it = exp_it->begin(); spec_it != exp_it->end(); ++spec_it)
			{
				//Foreign Key (Spectrum)
				query << tmp;
				//Intensity
				query << spec_it->getIntensity() << "','";
				//mz
				query << spec_it->getPosition()[0] << "'),";
			}
			db_con_.executeQuery(String(query.str()).substr(0,-1),result);
			
			// We know that all inserted peaks have IDs beginning from last_insert_id() (= ID of first inserted entry
			// of last insert operation), so we can insert Meta Information without actually fetching the ID
			UID insert_id = db_con_.getAutoId();
			for (typename ExperimentType::SpectrumType::Iterator spec_it = exp_it->begin(); spec_it != exp_it->end(); ++spec_it)
			{
				storeMetaInfo_("DATA_Peak", insert_id, *spec_it);
				insert_id++;
			}
			
			//----------------------------------------------------------------------------------------
			//---------------------------------- INSTRUMENT SETTINGS --------------------------------- 
			//----------------------------------------------------------------------------------------
			
			const InstrumentSettings & settings = exp_it->getInstrumentSettings();
			
			query.str("");
			
			if (new_entry)
			{
				query << "INSERT INTO META_InstrumentSettings SET fid_Spectrum=" << exp_it->getPersistenceId() << ",";
				end = "";
			}
			else
			{
				query << "SELECT id FROM META_InstrumentSettings WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				db_con_.executeQuery(query.str(),result);
				result.first();
				parent_id = result.value(0).asInt();
					
				query.str("");
				query << "UPDATE META_InstrumentSettings SET ";
				end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
			}

			query << "MZRangeBegin=" << settings.getMzRangeStart() << ",";
			query << "MZRangeEnd=" << settings.getMzRangeStop() << ",";
			query << "Polarity=" << (1u+settings.getPolarity()) << ",";
			query << "ScanMode=" << (1u+settings.getScanMode());
			query << end;
			
			db_con_.executeQuery(query.str(), result);
		
			if (new_entry) parent_id = db_con_.getAutoId();
			storeMetaInfo_("META_InstrumentSettings", parent_id, exp_it->getInstrumentSettings());
			
			//----------------------------------------------------------------------------------------
			//--------------------------------- ACQUISITIONINFO -------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			const AcquisitionInfo & info = exp_it->getAcquisitionInfo();
			
			query.str("");
			
			if (new_entry)
			{
				query << "INSERT INTO META_AcquisitionInfo SET fid_Spectrum=" << exp_it->getPersistenceId() << ",";
				end = "";
			}
			else
			{
				query << "SELECT id FROM META_AcquisitionInfo WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				db_con_.executeQuery(query.str(),result);
				result.first();
				acquisition_info_id = result.value(0).asInt();
					
				query.str("");
				query << "UPDATE META_AcquisitionInfo SET ";
				end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
			}

			query << "MethodOfCombination='" << info.getMethodOfCombination() << "'";
			query << end;
			
			db_con_.executeQuery(query.str(), result);
			if (new_entry)
			{
				acquisition_info_id = db_con_.getAutoId();
			}
		
			//----------------------------------------------------------------------------------------
			//----------------------------------- ACQUISITION ---------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			query.str("");
			deleteMetaInfo_("META_Acquisition", "fid_AcquisitionInfo=" + String(parent_id));
			query << "DELETE FROM META_Acquisition WHERE fid_AcquisitionInfo=" << parent_id;
			db_con_.executeQuery(query.str(), result);
				
			for (std::vector<Acquisition>::const_iterator info_it = info.begin(); info_it != info.end(); info_it++)
			{
				query.str("");
				query << "INSERT INTO META_Acquisition SET fid_AcquisitionInfo=" << acquisition_info_id << ",";
				query << "Number=" << info_it->getNumber();
				
				db_con_.executeQuery(query.str(), result);
				parent_id = db_con_.getAutoId();
				
				storeMetaInfo_("META_Acquisition", parent_id, *info_it);
			}
		}
	}

	template <class ExperimentType>
	void DBAdapter::loadExperiment(UID id, ExperimentType& exp, PeakFileOptions options)
	{
		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION --------------------------------------- 
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;
		
		std::stringstream query; // query to build
		String tmp;              // temporary data
		QSqlQuery result;        // place to store the query results in
		UID parent_id;					 // holds ID of parent data set
		
		query << "SELECT Type-1,Date,fid_MetaInfo,Description FROM META_MSExperiment WHERE id='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		//Experiment meta info
		exp.setType((ExperimentalSettings::ExperimentType)(result.value(0).asInt()));
		if (result.value(1).asDate().isValid())
		{
			Date d;
			d.set(result.value(1).asDate().toString(Qt::ISODate).ascii());
			exp.setDate(d);
		}
		exp.setComment(result.value(3).asString().ascii());
		loadMetaInfo_(result.value(2).asInt(),exp);
		
		// Sample
		Sample sample;
		query.str("");
		// finding root of recursive sample tree
		query << "SELECT id FROM META_Sample WHERE fid_MSExperiment='" << id << "' AND fid_Sample IS NULL";
		db_con_.executeQuery(query.str(),result);
		result.first();
		loadSample_ (result.value(0).asInt(), sample);
		exp.setSample(sample);
		
		// ContactPerson
		ContactPerson contact;
		query.str("");
		query << "SELECT PreName,LastName,Affiliation,Email,Comment,fid_MetaInfo FROM META_ContactPerson WHERE fid_MSExperiment='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		result.first();
		while(result.isValid())
		{
			contact.setFirstName(result.value(0).toString().ascii());
			contact.setLastName(result.value(1).toString().ascii());
			contact.setInstitution(result.value(2).toString().ascii());
			contact.setEmail(result.value(3).toString().ascii());
			contact.setContactInfo(result.value(4).toString().ascii());
			loadMetaInfo_(result.value(5).asInt(),contact);
			result.next();
			exp.getContacts().push_back(contact);
		}
		
		// HPLC
		query.str("");
		query << "SELECT id,InstrumentName,ColumnName,Description,Flux,Pressure,Temperature FROM META_HPLC WHERE fid_MSExperiment='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		parent_id=result.value(0).asInt();
		exp.getHPLC().setInstrument(result.value(1).toString().ascii());
		exp.getHPLC().setColumn(result.value(2).toString().ascii());
		exp.getHPLC().setComment(result.value(3).toString().ascii());
		exp.getHPLC().setFlux(result.value(4).asInt());
		exp.getHPLC().setPressure(result.value(5).asInt());
		exp.getHPLC().setTemperature(result.value(6).asInt());
		
		// Gradient*
		// I tried taking the big query apart in order to skip the double join, but this leads to
		// the problem of saving all requested keys in a vector in order to request the percentages (complex).
		// I'll still preserve the code in order to optimize in the future. Maybe I was just being blind. ;-)
		
		String last_name;
		bool timepoints_done = false;
		query.str("");
		/*
		query << "SELECT id,Name FROM META_GradientEluent WHERE fid_HPLC=" << parent_id;
		db_con_.executeQuery(query.str(),result);
		result.first();
		while(result.isValid())
		{
			exp.getHPLC().getGradient().addEluent(result.value(0).toString().ascii());
			result.next();
		}
		
		query.str("");
		query << "SELECT id,Time FROM META_GradientTime WHERE fid_HPLC=" << parent_id;
		db_con_.executeQuery(query.str(),result);
		result.first();
		while(result.isValid())
		{
			exp.getHPLC().getGradient().addTimepoint(result.value(0).asInt());
			result.next();
		}
		*/
		
		query << "SELECT Name,Time,Percentage FROM META_GradientEluent, META_GradientTime, META_GradientPercentage WHERE META_GradientEluent.fid_HPLC=" << parent_id << " AND fid_GradientEluent=META_GradientEluent.id AND fid_GradientTime=META_GradientTime.id";
		db_con_.executeQuery(query.str(),result);		
		result.first();
		
		if (result.isValid())
		{
			last_name = result.value(0).toString().ascii();
			exp.getHPLC().getGradient().addEluent(last_name);
		}
		
		while(result.isValid())
		{
			if (result.value(0).toString().ascii() != last_name)
			{
				exp.getHPLC().getGradient().addEluent(result.value(0).toString().ascii());
				timepoints_done = true;
			}
			
			if (timepoints_done == false)
			{
				exp.getHPLC().getGradient().addTimepoint(result.value(1).asInt());
			}
			
			exp.getHPLC().getGradient().setPercentage(result.value(0).toString().ascii(), result.value(1).asInt(), result.value(2).asInt());
			
			last_name = result.value(0).toString().ascii();
			result.next();
		}
		
		// Instrument
		query.str("");
		query << "SELECT id,Model,Vendor,Description,fid_MetaInfo FROM META_MSInstrument WHERE fid_MSExperiment='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		parent_id = result.value(0).asInt();
		exp.getInstrument().setModel(result.value(1).toString().ascii());
		exp.getInstrument().setVendor(result.value(2).toString().ascii());
		exp.getInstrument().setCustomizations(result.value(3).toString().ascii());
		loadMetaInfo_(result.value(4).asInt(),exp.getInstrument());
		
		// IonDetector
		query.str("");
		query << "SELECT AcquisitionMode-1,Type-1,Resolution,ADCSamplingFrequency,fid_MetaInfo FROM META_IonDetector WHERE fid_MSInstrument='" << parent_id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		exp.getInstrument().getIonDetector().setAcquisitionMode((IonDetector::AcquisitionMode) result.value(0).asInt());
		exp.getInstrument().getIonDetector().setType((IonDetector::Type) result.value(1).asInt());
		exp.getInstrument().getIonDetector().setResolution(result.value(2).toDouble());
		exp.getInstrument().getIonDetector().setADCSamplingFrequency(result.value(3).toDouble());
		loadMetaInfo_(result.value(4).asInt(),exp.getInstrument().getIonDetector());
		
		// IonSource
		query.str("");
		query << "SELECT InletType-1,IonizationMethod-1,IonizationMode-1,fid_MetaInfo FROM META_IonSource WHERE fid_MSInstrument='" << parent_id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		exp.getInstrument().getIonSource().setInletType((IonSource::InletType) result.value(0).asInt());
		exp.getInstrument().getIonSource().setIonizationMethod((IonSource::IonizationMethod) result.value(1).asInt());
		exp.getInstrument().getIonSource().setPolarity((IonSource::Polarity) result.value(2).toDouble());
		loadMetaInfo_(result.value(3).asInt(),exp.getInstrument().getIonSource());
		
		// MassAnalyzers
		MassAnalyzer analyzer;
		std::vector<MassAnalyzer> analyzers;
		query.str("");
		query << "SELECT Accuracy,FinalMSExponent,IsolationWidth,MagneticFieldStrength,ReflectronState-1,Resolution,ResolutionMethod-1,ResolutionType-1,ScanDirection-1,ScanFunction-1,ScanLaw-1,ScanRate,ScanTime,TandemScanningMethod-1,TOFPathLength,Type-1,fid_MetaInfo FROM META_MassAnalyzer WHERE fid_MSInstrument='" << parent_id << "'";
		db_con_.executeQuery(query.str(),result);
		
		result.first();
		while(result.isValid())
		{
			analyzer.setAccuracy(result.value(0).toDouble());
			analyzer.setFinalMSExponent(result.value(1).asInt());
			analyzer.setIsolationWidth(result.value(2).toDouble());
			analyzer.setMagneticFieldStrength(result.value(3).toDouble());
			analyzer.setReflectronState((MassAnalyzer::ReflectronState) result.value(4).asInt());
			analyzer.setResolution(result.value(5).toDouble());
			analyzer.setResolutionMethod((MassAnalyzer::ResolutionMethod) result.value(6).asInt());
			analyzer.setResolutionType((MassAnalyzer::ResolutionType) result.value(7).asInt());
			analyzer.setScanDirection((MassAnalyzer::ScanDirection) result.value(8).asInt());
			analyzer.setScanFunction((MassAnalyzer::ScanFunction) result.value(9).asInt());
			analyzer.setScanLaw((MassAnalyzer::ScanLaw) result.value(10).asInt());
			analyzer.setScanRate(result.value(11).toDouble());
			analyzer.setScanTime(result.value(12).toDouble());
			analyzer.setTandemScanMethod((MassAnalyzer::TandemScanningMethod) result.value(13).asInt());
			analyzer.setTOFTotalPathLength(result.value(14).toDouble());
			analyzer.setType((MassAnalyzer::AnalyzerType) result.value(15).asInt());
			loadMetaInfo_(result.value(16).asInt(), analyzer);
			
			analyzers.push_back(analyzer);
			result.next();
		}
		exp.getInstrument().setMassAnalyzers(analyzers);

		//id
		exp.setPersistenceId(id);

		// if we don't have to load the spectra, we're already done
		if (options.getMetadataOnly())
		{
			return;
		}

		//spectra
		query.str("");
		query << "SELECT id FROM DATA_Spectrum WHERE fid_MSExperiment=" << id;
		if (options.hasRTRange())
		{
			query << " AND RetentionTime > " << options.getRTRange().min() << " AND RetentionTime < " << options.getRTRange().max();
		}
		if (options.hasMSLevels())
		{
			const std::vector<int>& levels = options.getMSLevels();
			query << " AND (";
			for (std::vector<int>::const_iterator it = levels.begin(); it != levels.end(); it++)
			{
				query << "MSLevel=" << *it;
				if (it+1 != levels.end())
				{
					query << " OR ";
				}
			}
			query << ")";
		}
		query << " ORDER BY id ASC";
		
		db_con_.executeQuery(query.str(),result);
		exp.resize(result.size());
		UnsignedInt i = 0;
		result.first();
		while (result.isValid())
		{
			loadSpectrum(result.value(0).asInt(), exp[i], options);
			++i;
			result.next();
		}
	}

	template <class SpectrumType>
	void DBAdapter::loadSpectrum(UID id, SpectrumType& spec, PeakFileOptions options)
	{
		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION --------------------------------------- 
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;
		
		spec = SpectrumType();
		
		std::stringstream query;     // query to build
		QSqlQuery result;            // place to store the query results in
		InstrumentSettings settings; // stores settings that are read from DB
		UID parent_id;               // stores parent_id of Acquisition
		
		query << "SELECT Type-1,RetentionTime,MSLevel,Description,fid_MetaInfo,fid_File FROM DATA_Spectrum WHERE id='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		//Spectrum meta info
		spec.setType((SpectrumSettings::SpectrumType)(result.value(0).asInt()));
		spec.setRetentionTime(result.value(1).toDouble());
		spec.setMSLevel(result.value(2).asInt());
		spec.setComment(result.value(3).asString().ascii());
		loadMetaInfo_(result.value(4).asInt(),spec);
		loadFile_(result.value(5).asInt(),spec.getSourceFile());
		
		// Instrument settings
		query.str("");
		query << "SELECT MZRangeBegin, MZRangeEnd, Polarity-1, ScanMode-1, fid_MetaInfo FROM META_InstrumentSettings WHERE fid_Spectrum=" << id;
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		settings.setMzRangeStart(result.value(0).toDouble());
		settings.setMzRangeStop(result.value(1).toDouble());
		settings.setPolarity((IonSource::Polarity) (result.value(2).asInt()));
		settings.setScanMode((InstrumentSettings::ScanMode) (result.value(3).asInt()));
		spec.setInstrumentSettings(settings);
		loadMetaInfo_(result.value(4).asInt(),spec.getInstrumentSettings());

		// AcquisitionInfo
		query.str("");
		query << "SELECT id, MethodOfCombination FROM META_AcquisitionInfo WHERE fid_Spectrum=" << id;
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		spec.getAcquisitionInfo().setMethodOfCombination(result.value(1).asString().ascii());
		parent_id = result.value(0).asInt();
		
		// Acquisition
		query.str("");
		query << "SELECT Number,fid_MetaInfo FROM META_Acquisition WHERE fid_AcquisitionInfo='" << parent_id << "' ORDER BY id ASC";
		db_con_.executeQuery(query.str(),result);

		Acquisition acquisition;
		
		result.first();
		while(result.isValid())
		{
			acquisition.setNumber(result.value(0).asInt());
			loadMetaInfo_(result.value(1).toInt(), acquisition);
			spec.getAcquisitionInfo().push_back(acquisition);
			result.next();
		}
		
		// MetaInfoDescription
		query.str("");
		query << "SELECT Name, Description, fid_MetaInfo, fid_File FROM META_MetaInfoDescription WHERE fid_Spectrum=" << id;

		db_con_.executeQuery(query.str(), result);
		result.first();
		
		std::map<String,MetaInfoDescription> descriptions;
		MetaInfoDescription desc;
		
		while(result.isValid())
		{
			desc.setName(result.value(0).asString().ascii());
			desc.setComment(result.value(1).asString().ascii());
			loadMetaInfo_(result.value(2).asInt(), desc);
			loadFile_(result.value(3).asInt(),desc.getSourceFile());
			
			descriptions[result.value(0).asString().ascii()] = desc;
			result.next();
		}
		spec.setMetaInfoDescriptions (descriptions);

		//precursor
		if(spec.getMSLevel()>1)
		{
			query.str("");
			query << "SELECT mz,Intensity,Charge,ActivationMethod-1,ActivationEnergyUnit-1,ActivationEnergy,WindowSize,fid_MetaInfo FROM DATA_Precursor WHERE fid_Spectrum='" << id << "'";
			db_con_.executeQuery(query.str(),result);
			result.first();
			spec.getPrecursorPeak().getPosition()[0] = (result.value(0).toDouble());
			spec.getPrecursorPeak().setIntensity(result.value(1).toDouble());
			spec.getPrecursorPeak().setCharge(result.value(2).asInt());
			spec.getPrecursor().setActivationMethod((Precursor::ActivationMethod)(result.value(3).asInt()));
			spec.getPrecursor().setActivationEnergyUnit((Precursor::EnergyUnits)(result.value(4).asInt()));
			spec.getPrecursor().setActivationEnergy(result.value(5).toDouble());
			spec.getPrecursor().setWindowSize(result.value(6).toDouble());
			loadMetaInfo_(result.value(7).asInt(),spec.getPrecursor());
		}
		
		//Peaks
		query.str("");
		query << "SELECT mz,Intensity,fid_MetaInfo FROM DATA_Peak WHERE fid_Spectrum='" << id << "' ";
		if (options.hasMZRange())
		{
			query << " AND mz > " << options.getMZRange().min() << " AND mz < " << options.getMZRange().max();
		}
		if (options.hasIntensityRange())
		{
			query << " AND Intensity > " << options.getIntensityRange().min() << " AND Intensity < " << options.getIntensityRange().max();
		}
		query << " ORDER BY mz ASC";
		db_con_.executeQuery(query.str(),result);

		typename SpectrumType::PeakType p;
		result.first();
		while(result.isValid())
		{
			p.getPosition()[0] = result.value(0).toDouble();
			p.setIntensity(result.value(1).toDouble());
			loadMetaInfo_(result.value(2).toInt(), p);
			spec.push_back(p);
			result.next();
		}
		
		//id
		spec.setPersistenceId(id);
	}
}

#endif
