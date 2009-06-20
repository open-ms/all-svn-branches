// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
//#include <OpenMS/FORMAT/AnalysisXMLFile.h>

#include <ctime>
#include <vector>
#include <map>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_SpecLibSearcher SpecLibSearcher
 
 @brief Identifies peptides by spectral matching with a searchable spectral library.
 
 Hier noch Beschreibung

 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpecLibSearcher
	: public TOPPBase
	{
	public:
		TOPPSpecLibSearcher()
		: TOPPBase("SpecLibSearcher","Matches given spectra against a searchable spectral library.")
		{
		}
		
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input File ");
			setValidFormats_("in",StringList::create("mzData"));
			registerInputFile_("lib","<file>","","searchable spectral library(MSP formated)");
			registerOutputFile_("out","<file>","","Output File");
			setValidFormats_("out",StringList::create("IdXML"));
			registerDoubleOption_("precursor_mass_tolerance","<tolerance>",3,"Precursor mass tolerance, (Th)",false);
			registerIntOption_("precursor_mass_multiplier","<number>",100,"multipling this number with precursor_mass creates an integer",false,true);
		//	registerDoubleOption_("fragment_mass_tolerance","<tolerance>",0.3,"Fragment mass error",false);
			
      registerStringOption_("precursor_error_units", "<unit>", "Da", "parent monoisotopic mass error units", false);
     // registerStringOption_("fragment_error_units", "<unit>", "Da", "fragment monoisotopic mass error units", false);
      vector<String> valid_strings;
      valid_strings.push_back("Da");
      setValidStrings_("precursor_error_units", valid_strings);
     // setValidStrings_("fragment_error_units", valid_strings);
			registerIntOption_("min_precursor_charge", "<charge>", 1, "minimum precursor ion charge", false);
      registerIntOption_("max_precursor_charge", "<charge>", 3, "maximum precursor ion charge", false);
     	registerStringOption_("compare_function","<string>","ZhangSimilarityScore","function for similarity comparisson",false);
     PeakSpectrumCompareFunctor::registerChildren();
     setValidStrings_("compare_function",Factory<PeakSpectrumCompareFunctor>::registeredProducts());
     registerDoubleOption_("remove_peak_intensity_threshold","<threshold>",2.01,"All peaks of a query spectrum with intensities below <threshold> will be zeroed.",false);
     registerIntOption_("min_peak_number","<number>",5, "required mininum number of peaks for a query spectrum",false);
     registerDoubleOption_("filter_all_peaks_below_mz","<threshold>",520,"Remove spectra with almost no peaks above a certain m/z value. All query spectra with 95%+ of the total intensity below <m/z> will be removed.",false);	
     //registerStringOption_("fixed_modifications", "<mods>", "", "fixed modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048 currently no effect", false);
     // registerStringOption_("variable_modifications", "<mods>", "", "variable modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048", false);
		//	vielleicht des noch später??
			addEmptyLine_();
			addText_("");
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
		
			String in_spec = getStringOption_("in");
			String out = getStringOption_("out");
			String in_lib = getStringOption_("lib");
			String compare_function = getStringOption_("compare_function");
			Int precursor_mass_multiplier = getIntOption_("precursor_mass_multiplier");
			Real precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
			Int min_precursor_charge = getIntOption_("min_precursor_charge");
			Int max_precursor_charge = getIntOption_("max_precursor_charge");
			Real remove_peak_intensity_threshold = getDoubleOption_("remove_peak_intensity_threshold");
			UInt min_peak_number = getIntOption_("min_peak_number");
			Real filter_all_peaks_below_mz = getDoubleOption_("filter_all_peaks_below_mz");
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			//time counter
			time_t start_time = time(NULL);
			
			MzDataFile spectra;
			MSPFile spectral_library;
			
			RichPeakMap query, library;
			
			spectra.setLogType(log_type_);
			
			//spectrum which will be identified
			spectra.load(in_spec,query);
			
			//library containing already identified peptide spectra
			vector< PeptideIdentification > ids;			
			spectral_library.load(in_lib,ids,library);
			
			//*******
			//building map for precursor mass search
			//*******
			map<Size, vector<RichPeakSpectrum*> > MSLibrary;
			RichPeakMap::iterator s;
			vector<PeptideIdentification>::iterator i;
			for( s = library.begin(), i = ids.begin() ; s < library.end();++s,++i)
			{
				Size pm = (Size)(*i).getHits()[0].getSequence().getMonoWeight()*precursor_mass_multiplier;
				map<Size,vector<RichPeakSpectrum*> >::iterator found;
				found = MSLibrary.find(pm);				
				if(found != MSLibrary.end())
				{
					PeptideIdentification& translocate_id = (*i);
					(*s).getPeptideIdentifications().push_back(translocate_id);
					RichPeakSpectrum* translocate_s = &(*s);
					found->second.push_back(translocate_s);
				}
				else
				{
					vector<RichPeakSpectrum*> tmp;
					PeptideIdentification& translocate_id = (*i);
					RichPeakSpectrum* translocate_s = &(*s);
					(*s).getPeptideIdentifications().push_back(translocate_id);
					tmp.push_back(translocate_s);
					MSLibrary.insert(make_pair(pm,tmp));
				}
			}
			//compare function
			PeakSpectrumCompareFunctor* comparor = Factory<PeakSpectrumCompareFunctor>::create(compare_function);
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			DoubleReal score;
			//Will hold valuable hits
			vector<PeptideIdentification> peptide_ids;
			vector<ProteinIdentification> protein_ids;
			Real min_pm,max_pm;
			
			ProteinIdentification prot_id;
			//Parameters of identificaion
			prot_id.setIdentifier("test");
			prot_id.setSearchEngineVersion("SpecLibSearcher");
			prot_id.setDateTime(DateTime::now());
			prot_id.setScoreType(compare_function);
			ProteinIdentification::SearchParameters searchparam;
			searchparam.charges += (min_precursor_charge);
			searchparam.charges += " - ";
			searchparam.charges += (max_precursor_charge);
			searchparam.precursor_tolerance = precursor_mass_tolerance;
			prot_id.setSearchParameters(searchparam);
			for(UInt j = 0; j < query.size(); ++j)
			{
			//Set identifier for each identifications
				PeptideIdentification pep_id;
				pep_id.setIdentifier("test");
				pep_id.setScoreType(compare_function);
				ProteinHit pr_hit;
				pr_hit.setAccession(j);
				prot_id.insertHit(pr_hit);
				
				Real total_intensity = 0;
				Real intensity_below_mz = 0;
				
				//RichPeak1D to Peak1D transformation for the compare function query 
				PeakSpectrum quer;
				bool peak_ok = true;
				for(UInt k = 0; k < query[j].size(); ++k)
				{
					Peak1D peak;
					if(query[j][k].getIntensity() >  remove_peak_intensity_threshold)
					{
						total_intensity +=  query[j][k].getIntensity();
						peak.setIntensity(query[j][k].getIntensity());
						peak.setMZ(query[j][k].getMZ());
						if (peak.getMZ() < filter_all_peaks_below_mz) 
						{
        			intensity_below_mz += peak.getIntensity();
      			}
						peak.setPosition(query[j][k].getPosition());
						quer.push_back(peak);
					}
				}
				//if not enough peaks in the specturm pass that one out
				//Remove spectra with almost no peaks above a certain m/z value. 
				//All query spectra with 95%+ of the total intensity below <m/z> will be removed.
				if(quer.size() >= min_peak_number && intensity_below_mz < 0.95 * total_intensity)
				{
					peak_ok = true;
				}
				else
				{
					peak_ok = false;
				}
				if(peak_ok)
				{
					for(Int charge = min_precursor_charge; charge <= max_precursor_charge; ++charge)
					{
						min_pm = ((query[j].getPrecursors()[0].getMZ() - precursor_mass_tolerance)*charge-(charge*1.00728)) *precursor_mass_multiplier;
						max_pm = ((query[j].getPrecursors()[0].getMZ() + precursor_mass_tolerance)*charge-(charge*1.00728))*precursor_mass_multiplier;
				
					//SEARCH
						for(Size pm = (Size)min_pm; pm <= ((Size)max_pm)+1; ++pm)
						{
							map<Size, vector<RichPeakSpectrum*> >::iterator found;
							found = MSLibrary.find(pm);
							if(found != MSLibrary.end())
							{
								vector<RichPeakSpectrum*>& library = found->second;
								for(Size i = 0; i < library.size(); ++i)
								{
									Real this_pm  = library[i]->getPeptideIdentifications()[0].getHits()[0].getSequence().getMonoWeight()*precursor_mass_multiplier;
									if(this_pm >= min_pm && max_pm >= this_pm)
									{
										PeakSpectrum librar;
										//library entry transformation 
										for(UInt l = 0; l< library[i]->size(); ++l)
										{
											Peak1D peak;
											if((*library[i])[l].getIntensity() >  remove_peak_intensity_threshold)
											{
												peak.setIntensity((*library[i])[l].getIntensity());									
												peak.setMZ((*library[i])[l].getMZ());
												peak.setPosition((*library[i])[l].getPosition());
												librar.push_back(peak);
											}
										}
										score = (*comparor)(quer,librar);
											PeptideHit hit(library[i]->getPeptideIdentifications()[0].getHits()[0]);
											DataValue RT(library[i]->getRT());
											DataValue MZ(library[i]->getPrecursors()[0].getMZ());
											hit.setMetaValue("RT",RT);
											hit.setMetaValue("MZ",MZ);
											hit.setScore(score);
											hit.addProteinAccession(pr_hit.getAccession());
											pep_id.insertHit(hit);
									}
								}
							}
						}
					}
				}
				pep_id.setHigherScoreBetter(true);
				pep_id.sort();
				peptide_ids.push_back(pep_id);
			}
			protein_ids.push_back(prot_id);			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			IdXMLFile id_xml_file;
			id_xml_file.store(out,protein_ids,peptide_ids);
			time_t end_time = time(NULL);
			cout<<"Time needed: "<<difftime(end_time,start_time)<<" seconds"<<"\n";  
			return EXECUTION_OK;
		}
		
	};




int main( int argc, const char** argv )
{
    TOPPSpecLibSearcher tool;
    return tool.main(argc,argv);
}

/// @endcond