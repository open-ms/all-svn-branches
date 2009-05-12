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
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
//#include <OpenMS/FORMAT/AnalysisXMLFile.h>



#include <vector>
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
		: TOPPBase("SpectralFinder","Matches given spectra against a searchable spectral library.")
		{
		}
		
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("spec","<file>","","input raw data file ");
			setValidFormats_("spec",StringList::create("mzData"));
			registerInputFile_("lib","<file>","","searchable spectral library");
			//setValidFormats_("lib",StringList::create("msp"));
			registerOutputFile_("out","<file>","","output matches");
			setValidFormats_("out",StringList::create("IdXML"));

			addEmptyLine_();
			addText_("");
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
		
			String in_spec = getStringOption_("spec");
			String out = getStringOption_("out");
			String in_lib = getStringOption_("lib");
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
			
			//compare function
			PeakSpectrumCompareFunctor* comparor = ZhangSimilarityScore::create();
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			DoubleReal score;
			//Will hold valuable hits
			vector<PeptideIdentification> peptide_ids;
			vector<ProteinIdentification> protein_ids;
			for(UInt j = 0; j < query.size(); ++j)
			{
				PeptideIdentification pep_id;
				ProteinIdentification prot_id;
				pep_id.setIdentifier(j);
				prot_id.setIdentifier(j);
				for(UInt i =0; i< library.size();++i)
				{
					
					PeakSpectrum quer;
					PeakSpectrum librar;
					for(UInt k = 0; k < query[j].size(); ++k)
					{
						Peak1D peak;
						peak.setIntensity(query[j][k].getIntensity());
						peak.setMZ(query[j][k].getMZ());
						peak.setPosition(query[j][k].getPosition());
						quer.push_back(peak);
					}
					for(UInt l = 0; l< library[i].size(); ++l)
					{
						Peak1D peak;
						peak.setIntensity(library[i][l].getIntensity());
						peak.setMZ(library[i][l].getMZ());
						peak.setPosition(library[i][l].getPosition());
						librar.push_back(peak);
					}
					score = (*comparor)(quer,librar);
					if(score > 0.5)
					{
						PeptideHit hit(ids[i].getHits()[0]);
						hit.setScore(score);

						pep_id.insertHit(hit);
					}
				}
				peptide_ids.push_back(pep_id);
				protein_ids.push_back(prot_id);
			}
			
			
			
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