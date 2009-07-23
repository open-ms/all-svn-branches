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
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <ctime>
#include <vector>
#include <map>
#include <cmath>
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
			registerInputFile_("lib","<file>","","searchable spectral library(MSP format)");
			registerOutputFile_("out","<file>","","Output File");
			setValidFormats_("out",StringList::create("IdXML"));
			registerDoubleOption_("precursor_mass_tolerance","<tolerance>",3,"Precursor mass tolerance, (Th)",false);
		//	registerIntOption_("precursor_mass_multiplier","<number>",10,"multipling this number with precursor_mass creates an integer",false,true);
		//	registerDoubleOption_("fragment_mass_tolerance","<tolerance>",0.3,"Fragment mass error",false);
			
   //   registerStringOption_("precursor_error_units", "<unit>", "Da", "parent monoisotopic mass error units", false);
     // registerStringOption_("fragment_error_units", "<unit>", "Da", "fragment monoisotopic mass error units", false);
     // vector<String> valid_strings;
     // valid_strings.push_back("Da");
     // setValidStrings_("precursor_error_units", valid_strings);
     // setValidStrings_("fragment_error_units", valid_strings);
			//registerIntOption_("min_precursor_charge", "<charge>", 1, "minimum precursor ion charge", false);
    //  registerIntOption_("max_precursor_charge", "<charge>", 3, "maximum precursor ion charge", false);
     	registerStringOption_("compare_function","<string>","ZhangSimilarityScore","function for similarity comparisson",false);
     PeakSpectrumCompareFunctor::registerChildren();
     setValidStrings_("compare_function",Factory<PeakSpectrumCompareFunctor>::registeredProducts());

     addEmptyLine_();    
     addText_("Filtering options. Most are especially useful when the query spectra are raw.");
     registerIntOption_("min_peaks","<number>",5, "required mininum number of peaks for a query spectrum",false);     
     registerDoubleOption_("remove_peaks_below_threshold","<threshold>",2.01,"All peaks of a query spectrum with intensities below <threshold> will be zeroed.",false);    
     registerIntOption_("max_peaks","<number>",150,"Use only the top <number> of peaks.",false);
     registerIntOption_("cut_peaks_below","<number>",1000,"Remove all peaks which are lower than 1/<number> of the highest peaks. Default equals all peaks which are lower than 0.001 of the maximum intensity peak",false);
     
     //registerStringOption_("fixed_modifications", "<mods>", "", "fixed modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048 currently no effect", false);
     registerStringList_("variable_modifications", "<mods>", StringList::create(""), "variable modifications, specified using PSI-MOD terms, e.g. MOD:01214 MOD:00048", false);
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
			Int precursor_mass_multiplier = 10;//getIntOption_("precursor_mass_multiplier");
			Real precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
			//Int min_precursor_charge = getIntOption_("min_precursor_charge");
			//Int max_precursor_charge = getIntOption_("max_precursor_charge");
			Real remove_peaks_below_threshold = getDoubleOption_("remove_peaks_below_threshold");
			UInt min_peaks = getIntOption_("min_peaks");
			UInt max_peaks= getIntOption_("max_peaks");
			Int cut_peaks_below = getIntOption_("cut_peaks_below");
			StringList variable_modifications = getStringList_("variable_modifications");
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			time_t prog_time = time(NULL);
			MzDataFile spectra;
			MSPFile spectral_library;
			
			RichPeakMap query, library;
			
			spectra.setLogType(log_type_);
			
			//spectrum which will be identified
			spectra.load(in_spec,query);
			
			//library containing already identified peptide spectra
			vector< PeptideIdentification > ids;			
			spectral_library.load(in_lib,ids,library);
			
			//-------------------------------------------------------------
			// Write parameters to ProteinIdentifcation
			//-------------------------------------------------------------
			ProteinIdentification prot_id;
			//Parameters of identificaion
			prot_id.setIdentifier("test");
			prot_id.setSearchEngineVersion("SpecLibSearcher");
			prot_id.setDateTime(DateTime::now());
			prot_id.setScoreType(compare_function);
			ProteinIdentification::SearchParameters searchparam;
		//	searchparam.charges += (min_precursor_charge);
		//	searchparam.charges += " - ";
		//	searchparam.charges += (max_precursor_charge);
			searchparam.precursor_tolerance = precursor_mass_tolerance;
			prot_id.setSearchParameters(searchparam);
			time_t start_build_time = time(NULL);
			//-------------------------------------------------------------
			//building map for faster search
			//-------------------------------------------------------------
			map<Size, vector<PeakSpectrum> > MSLibrary;
			{
			RichPeakMap::iterator s;
			vector<PeptideIdentification>::iterator i;
			ModificationsDB* mdb = ModificationsDB::getInstance();
											
			for(s = library.begin(), i = ids.begin() ; s < library.end();++s,++i)
			{
				DoubleReal precursor_MZ = (*s).getPrecursors()[0].getMZ();
				Size MZ_multi = (Size)precursor_MZ*precursor_mass_multiplier;
				map<Size,vector<PeakSpectrum> >::iterator found;
				found = MSLibrary.find(MZ_multi);
				
				PeakSpectrum librar;
				bool variable_modifications_ok = true;
				const AASequence& aaseq= i->getHits()[0].getSequence();
				//variable Modification
				
				if(aaseq.isModified() && !variable_modifications.empty())
				{
					for(Size i = 0; i< aaseq.size(); ++i)
					{
						if(aaseq.isModified(i))
						{
							const	Residue& mod  = aaseq.getResidue(i);
							
							for(Size s = 0; s < variable_modifications.size();++s)
							{
								if(mod.getOneLetterCode() ==mdb->getModification(variable_modifications[s]).getOrigin() && variable_modifications[s] != mod.getModification())
								{
									variable_modifications_ok = false;
									break;
								}	
							}
						}
					}
				}
				if(variable_modifications_ok)
				{
					PeptideIdentification& translocate_pid = *i;
					librar.getPeptideIdentifications().push_back(translocate_pid);		
					librar.setPrecursors(s->getPrecursors());
					//library entry transformation 
					for(UInt l = 0; l< s->size(); ++l)
					{
						Peak1D peak;
						if((*s)[l].getIntensity() >  remove_peaks_below_threshold)
						{
						//	if((*s)[l].getMZ() > precursor_MZ || (*s)[l].getMZ() < precursor_MZ - 1 )
						//	{
								peak.setIntensity(sqrt((*s)[l].getIntensity()));
						//	}
						//	else
						//	{
						//		peak.setIntensity(sqrt(0.2*(*s)[l].getIntensity()));
						//	}
							peak.setMZ((*s)[l].getMZ());
							peak.setPosition((*s)[l].getPosition());
							librar.push_back(peak);
						}	
					}
					if(found != MSLibrary.end())
					{
						found->second.push_back(librar);
					}
					else
					{
						vector<PeakSpectrum> tmp;
						tmp.push_back(librar);
						MSLibrary.insert(make_pair(MZ_multi,tmp));
					}
				}
			}
			}
			time_t end_build_time = time(NULL);
			//compare function
			PeakSpectrumCompareFunctor* comparor = Factory<PeakSpectrumCompareFunctor>::create(compare_function);
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			DoubleReal score;
			//Will hold valuable hits
			vector<PeptideIdentification> peptide_ids;
			vector<ProteinIdentification> protein_ids;
			time_t start_time = time(NULL);
			for(UInt j = 0; j < query.size(); ++j)
			{
			//Set identifier for each identifications
				PeptideIdentification pid;
				pid.setIdentifier("test");
				pid.setScoreType(compare_function);
				ProteinHit pr_hit;
				pr_hit.setAccession(j);
				prot_id.insertHit(pr_hit);
				//RichPeak1D to Peak1D transformation for the compare function query 
				PeakSpectrum quer;
				bool peak_ok = true;
				query[j].sortByIntensity(true);
				DoubleReal min_high_intensity = 0;
				if(!query[j].empty())
				{
					min_high_intensity = (1/cut_peaks_below)*query[j][0].getIntensity();
				}
				query[j].sortByPosition();
				for(UInt k = 0; k < query[j].size() && k < max_peaks; ++k)
				{
					if(query[j][k].getIntensity() >  remove_peaks_below_threshold && query[j][k].getIntensity() >= min_high_intensity )
					{
						Peak1D peak;
			//			if(query[j][k].getMZ() > query[j].getPrecursors()[0].getMZ() ||query[j][k].getMZ() < query[j].getPrecursors()[0].getMZ()-14 )
				//		{
							peak.setIntensity(sqrt(query[j][k].getIntensity()));
			//			}
				//		else
				//		{
				//			peak.setIntensity(sqrt(0.2*query[j][k].getIntensity()));
					//	}
						peak.setMZ(query[j][k].getMZ());

						peak.setPosition(query[j][k].getPosition());
						quer.push_back(peak);
					}
				}
				quer.sortByPosition();	
				if(quer.size() >= min_peaks)
				{
					peak_ok = true;
				}
				else
				{
					peak_ok = false;
				}
				if(peak_ok)
				{
					Real min_MZ = (query[j].getPrecursors()[0].getMZ() - precursor_mass_tolerance) *precursor_mass_multiplier;
					Real max_MZ = (query[j].getPrecursors()[0].getMZ() + precursor_mass_tolerance) *precursor_mass_multiplier;
					//SEARCH
					for(Size mz = (Size)min_MZ; mz <= ((Size)max_MZ)+1; ++mz)
					{
						map<Size, vector<PeakSpectrum> >::iterator found;
							found = MSLibrary.find(mz);
							if(found != MSLibrary.end())
							{
								vector<PeakSpectrum>& library = found->second;
								for(Size i = 0; i < library.size(); ++i)
								{
									Real this_MZ  = library[i].getPrecursors()[0].getMZ()*precursor_mass_multiplier;
									if(this_MZ >= min_MZ && max_MZ >= this_MZ )
									{
										PeptideHit hit = library[i].getPeptideIdentifications()[0].getHits()[0];
										PeakSpectrum& librar = library[i];
										//Special treatment for SpectraST score as it computes a score based on the whole library
										if(compare_function == "SpectraSTSimilarityScore")
										{
											SpectraSTSimilarityScore* sp= static_cast<SpectraSTSimilarityScore*>(comparor);
											BinnedSpectrum quer_bin = sp->transform(quer);
											BinnedSpectrum librar_bin = sp->transform(librar);
											score = (*sp)(quer_bin,librar_bin);
											double dot_bias = sp->dot_bias(quer_bin,librar_bin,score);
											hit.setMetaValue("DOTBIAS",dot_bias);
										}
										else 
										{
											if(compare_function =="CompareFouriertransform")
											{
												CompareFouriertransform* ft = static_cast<CompareFouriertransform*>(comparor);
												ft->transform(quer);
												ft->transform(librar);
											}
											score = (*comparor)(quer,librar);
										}
										
										DataValue RT(library[i].getRT());
										DataValue MZ(library[i].getPrecursors()[0].getMZ());
										hit.setMetaValue("RT",RT);
										hit.setMetaValue("MZ",MZ);
										hit.setScore(score);
										hit.addProteinAccession(pr_hit.getAccession());
										pid.insertHit(hit);
									}
								}
							}
						}
					}
				pid.setHigherScoreBetter(true);
				pid.sort();
				if(compare_function =="SpectraSTSimilarityScore")
				{

					if(!pid.empty() && !pid.getHits().empty())
					{
						vector<PeptideHit> final_hits;					
						final_hits.resize(pid.getHits().size());
						SpectraSTSimilarityScore* sp= static_cast<SpectraSTSimilarityScore*>(comparor);
						double delta_D = sp->delta_D(pid.getHits()[0].getScore(), pid.getHits()[1].getScore());
						for(Size s = 0; s < pid.getHits().size();++s)
						{
							final_hits[s] = pid.getHits()[s];
							final_hits[s].setScore(sp->compute_F(pid.getHits()[s].getScore(),delta_D,pid.getHits()[s].getMetaValue("DOTBIAS")));
							final_hits[s].removeMetaValue("DOTBIAS");
						}
						pid.setHits(final_hits);
						pid.sort();
						pid.setMetaValue("MZ",query[j].getPrecursors()[0].getMZ());
						pid.setMetaValue("RT",query[j].getRT());
					}
				}
				peptide_ids.push_back(pid);
			}
			protein_ids.push_back(prot_id);			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			IdXMLFile id_xml_file;
			id_xml_file.store(out,protein_ids,peptide_ids);
			time_t end_time = time(NULL);
		  cout<<"Time needed for preprocessing data: "<<(end_build_time - start_build_time)<<"\n";
			cout<<"Search time: "<<difftime(end_time,start_time)<<" seconds"<<"\n";
			cout<<"Total time: "<<difftime(end_time,prog_time)<<" secconds"<<"\n";
			return EXECUTION_OK;
		}
		
	};




int main( int argc, const char** argv )
{
    TOPPSpecLibSearcher tool;
    return tool.main(argc,argv);
}

/// @endcond
