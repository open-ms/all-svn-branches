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
			//registerDoubleOption_("precursor_mass","<mass>",-1, "Precursor mass",false,true);
			registerDoubleOption_("precursor_mass_tolerance","<tolerance>",1.5,"Precursor mass tolerance",false);
			registerIntOption_("precursor_mass_multiplier","<number>",10,"multipling this number with precursor_mass creates an integer",false,true);
			registerDoubleOption_("fragment_mass_tolerance","<tolerance>",0.3,"Fragment mass error",false);
			
      registerStringOption_("precursor_error_units", "<unit>", "Da", "parent monoisotopic mass error units", false);
      registerStringOption_("fragment_error_units", "<unit>", "Da", "fragment monoisotopic mass error units", false);
      vector<String> valid_strings;
      valid_strings.push_back("ppm"); 
      valid_strings.push_back("Da");
      setValidStrings_("precursor_error_units", valid_strings);
      setValidStrings_("fragment_error_units", valid_strings);
			registerIntOption_("min_precursor_charge", "<charge>", 1, "minimum precursor ion charge", false);
      registerIntOption_("max_precursor_charge", "<charge>", 3, "maximum precursor ion charge", false);
     	registerStringOption_("compare_function","<string>","ZhangSimilarityScore","function for similarity comparisson",false);
     PeakSpectrumCompareFunctor::registerChildren();
     setValidStrings_("compare_function",Factory<PeakSpectrumCompareFunctor>::registeredProducts());
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
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			MzDataFile spectra;
			MSPFile spectral_library;
			vector< PeptideIdentification > ids;
			
			RichPeakMap query, library;
			
			spectra.setLogType(log_type_);
			
			//spectrum which will be identified
			spectra.load(in_spec,query);
			
			//library containing already identified peptide spectra
			spectral_library.load(in_lib,ids,library);
			
			//building hashmap for precursor mass search
			map<Size, RichPeakMap> MSLibrary;
			RichPeakMap::iterator s;
			vector<PeptideIdentification>::iterator i;
			for( s = library.begin(), i = ids.begin() ; s < library.end();++s,++i)
			{
				Size pm = precursor_mass_multiplier * (*s).getPrecursors()[0].getMZ();
				pm = floor(pm);
				//cout<<"pm: "<< pm<<endl;
				if(MSLibrary.find(pm)!= MSLibrary.end())
				{
					(*s).getPeptideIdentifications().push_back(*i);
					MSLibrary.find(pm)->second.push_back(*s);
				}
				else
				{
					RichPeakMap tmp;
					(*s).getPeptideIdentifications().push_back(*i);
					tmp.push_back(*s);
					MSLibrary.insert(make_pair(pm,tmp));
				}
			}
			//cout<<"Erstellung lief durch"<<endl;
			//cout<<"MSLibrary.size(): " <<(int) MSLibrary.size()<<endl;

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

			for(UInt j = 0; j < query.size(); ++j)
			{
			//cout<<"j: "<<j <<endl;
				PeptideIdentification pep_id;
				pep_id.setIdentifier(j);
							ProteinIdentification prot_id;
			prot_id.setIdentifier(j);
				
				min_pm = query[j].getPrecursors()[0].getMZ() - precursor_mass_tolerance;
				max_pm = query[j].getPrecursors()[0].getMZ() + precursor_mass_tolerance;
				min_pm *= precursor_mass_multiplier;
				max_pm *= precursor_mass_multiplier;
				//RichPeak1D to Peak1D transformation for the compare function query 
				PeakSpectrum quer;
				for(UInt k = 0; k < query[j].size(); ++k)
				{
				//	cout<<"k" <<k<<endl;
					Peak1D peak;
					peak.setIntensity(query[j][k].getIntensity());
					peak.setMZ(query[j][k].getMZ());
					peak.setPosition(query[j][k].getPosition());
					quer.push_back(peak);
				}
				//SEARCH
				for(Size pm = floor(min_pm); pm < ceil(max_pm)-1; ++pm)
				{
				//cout<<"pm: "<<pm<<endl;
					if(MSLibrary.find(pm)!= MSLibrary.end())
					{
						RichPeakMap library = MSLibrary.find(pm)->second;
						for(Size i = 0; i < library.size(); ++i)
						{
				//		cout<<"i: "<<i<<endl;
							PeakSpectrum librar;
							//library entry transformation 
				//			cout<<"library[i].size(): "<<library[i].size()<<endl;
							for(UInt l = 0; l< library[i].size(); ++l)
							{
			//				cout<<"l: "<<l<<endl;
								Peak1D peak;
								peak.setIntensity(library[i][l].getIntensity());
								peak.setMZ(library[i][l].getMZ());
								peak.setPosition(library[i][l].getPosition());
								librar.push_back(peak);
							}
			//				cout<<"1"<<endl;
							score = (*comparor)(quer,librar);
				//			cout<<"2"<<endl;
					//		cout<<"große des dings"<<library[i].getPeptideIdentifications().size()<<endl;
							PeptideHit hit(library[i].getPeptideIdentifications()[0].getHits()[0]);
				//			cout<<"3"<<endl;
							hit.setScore(score);
				//			cout<<"4"<<endl;
							pep_id.insertHit(hit);
				//			cout<<"5"<<endl;
						}
					}
				}
			//	cout<<"hinterm score"<<endl;
				pep_id.setHigherScoreBetter(true);
				pep_id.sort();
				peptide_ids.push_back(pep_id);
				protein_ids.push_back(prot_id);
			}

				//!!!old
				/*for(UInt i =0; i< library.size();++i)
				{
					PeakSpectrum quer;
					PeakSpectrum librar;
					//RichPeak1D to Peak1D transformation for the compare function
					//query 
					for(UInt k = 0; k < query[j].size(); ++k)
					{
						Peak1D peak;
						peak.setIntensity(query[j][k].getIntensity());
						peak.setMZ(query[j][k].getMZ());
						peak.setPosition(query[j][k].getPosition());
						quer.push_back(peak);
					}
					//library 
					for(UInt l = 0; l< library[i].size(); ++l)
					{
						Peak1D peak;
						peak.setIntensity(library[i][l].getIntensity());
						peak.setMZ(library[i][l].getMZ());
						peak.setPosition(library[i][l].getPosition());
						librar.push_back(peak);
					}
					//write hits down
					score = (*comparor)(quer,librar);
						PeptideHit hit(ids[i].getHits()[0]);
						hit.setScore(score);
						pep_id.insertHit(hit);
				}
				pep_id.setHigherScoreBetter(true);
				pep_id.sort();
				peptide_ids.push_back(pep_id);
				protein_ids.push_back(prot_id);
			}*/
			
			
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			IdXMLFile id_xml_file;
			id_xml_file.store(out,protein_ids,peptide_ids);
			
			return EXECUTION_OK;
		}
		
	};




int main( int argc, const char** argv )
{
    TOPPSpecLibSearcher tool;
    return tool.main(argc,argv);
}

/// @endcond