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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/ROCCurve.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <map>
#include <list>
#include <algorithm>
#include <cmath> // needed to check whether a score is nan  =not a number
using namespace OpenMS;
using namespace std;
using namespace OpenMS::Math;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_DataForROCCurve DataForROCCurve
 
 @brief generates an File prepared for ROC curves. At the moment only for spectral matching
 
 Hier noch Beschreibung

 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDataForROCCurve
	: public TOPPBase
	{
	public:
		TOPPDataForROCCurve()
		: TOPPBase("DataForROCCurve","cret")
		{
		}
		
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("fasta","<file>","","fasta file with already identified peptides");
		//	registerInputFile_("MSP","<file>","","spectral library which was used for matching");
			registerInputFile_("in","<file>","","file which stores scores of the last search,IdXML or pepXML");
			registerOutputFile_("out","<file>","","Output file");			
			registerIntOption_("curve_resolution","<number>",10,"points which represent the ROCCurve",false);
			registerIntOption_("highest","<number>",1,"1 if only the highest value should be taken into account",false);
			addEmptyLine_();
			addText_("");
		}
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			String fasta = getStringOption_("fasta");
		//	String MSP = getStringOption_("MSP");
			String in = getStringOption_("in");
			String out = getStringOption_("out");
			Int curve_resolution = getIntOption_("curve_resolution");
			Int highest = getIntOption_("highest");
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			FASTAFile fasta_file;
			vector<FASTAFile::FASTAEntry> fasta_entries;
			fasta_file.load(fasta,fasta_entries);
			
			vector<ProteinIdentification> idxml_proteinIdentification;
			vector<PeptideIdentification> idxml_peptideIdentification;
			FileHandler fh;
			FileTypes::Type in_type = fh.getType(in);
			if(in_type == FileTypes::UNKNOWN)
			{
				writeLog_("Warning: Could not determine input file type!");
			}
			else if(in_type == FileTypes::PEPXML)
			{
				cout<<"Spectrast ROC-Curve\n";
				PepXMLFile pepxml_file;
				String name;
				MSExperiment<> exp;
				ProteinIdentification pi;
				pepxml_file.load(in,pi,idxml_peptideIdentification,name,exp);
			
			}
			else if(in_type == FileTypes::IDXML)
			{
				cout<<"SpecLibSearcher ROC-Curve\n";
				IdXMLFile id_xml_file;
				id_xml_file.load(in, idxml_proteinIdentification, idxml_peptideIdentification);
			}

			
		//	MSPFile msp_file;
		//	RichPeakMap msp_mspexperiment;
		//	vector<PeptideIdentification> msp_peptideIdentification;
		//	msp_file.load(MSP, msp_peptideIdentification,msp_mspexperiment);
			
			
			//-------------------------------------------------------------
			// prepare for calculations
			//-------------------------------------------------------------
			//map containing peptides with highest score
			map<String, DoubleReal> IdXMLPeptides;
			//list containing peptides which are in MSPFile but not in IdXML
			//list<String> MSPPeptides;
			//first these which wher scored
			for(vector<PeptideIdentification>::iterator it = idxml_peptideIdentification.begin(); it< idxml_peptideIdentification.end();++it)
			{
				if(!highest)
				{
					if(!it->empty())
					{
						for(vector<PeptideHit>::const_iterator iter = it->getHits().begin(); iter < it->getHits().end();++iter)
						{
							map<String, DoubleReal>::iterator sequence;
							sequence = IdXMLPeptides.find(iter->getSequence().toUnmodifiedString());
							if(sequence != IdXMLPeptides.end() || isnan(iter->getScore()) || sequence->second > iter->getScore())
							{
								//DO NOTHING
							}
							else
							{
								IdXMLPeptides.insert(pair<String,DoubleReal>(iter->getSequence().toUnmodifiedString(),iter->getScore()));
							}
						}
					}
				}
				else
				{
					if(!it->empty())
					{
						
						
						it->setHigherScoreBetter(true);
						it->sort();
						map<String, DoubleReal>::iterator sequence;
						sequence = IdXMLPeptides.find(it->getHits()[0].getSequence().toUnmodifiedString());
						if(sequence != IdXMLPeptides.end()) 
						{
							if(sequence->second < it->getHits()[0].getScore())
							{
								IdXMLPeptides[it->getHits()[0].getSequence().toUnmodifiedString()] = it->getHits()[0].getScore();		
							}
						}
						else
						{
							IdXMLPeptides.insert(pair<String,DoubleReal>(it->getHits()[0].getSequence().toUnmodifiedString(),it->getHits()[0].getScore()));
						}
					}
				}
			}
			//now those which are in the library but nowhere else
	/*		for(vector<PeptideIdentification>::iterator it = msp_peptideIdentification.begin(); it < msp_peptideIdentification.end(); ++it)
			{
				if(!it->empty())
				{
					for(vector<PeptideHit>::iterator iter = it->getHits().begin(); iter < it->getHits().end();++it)
					{
						if(IdXMLPeptides.find(iter->getSequence().toString()) == IdXMLPeptides.end())
						{
							MSPPeptides.push_back(iter->getSequence().toString());
						}
					}
				}
			}*/
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			//generate ROC Curve
			
			ROCCurve ROC;
			Int true_count = 0;
			Int false_count = 0;
	//		cout<<"False:"<<endl;
			for(map<String,DoubleReal>::iterator idxmliter = IdXMLPeptides.begin(); idxmliter != IdXMLPeptides.end(); ++idxmliter)
			{
				
				bool exists  = false;
				for(vector<FASTAFile::FASTAEntry>::iterator fastaiter = fasta_entries.begin(); fastaiter < fasta_entries.end();++fastaiter)
				{
					if(fastaiter->sequence.hasSubstring(idxmliter->first))
					{
						exists  = true;
						continue;
					}
				}
				if(exists)
				{
					ROC.insertPair(idxmliter->second,true);
					++true_count;
				}
				else
				{
					ROC.insertPair(idxmliter->second,false);
					//cout<<idxmliter->first<<"\n";
					++false_count;
				}
			}

			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
					TextFile file;
					String temp;
					temp = "All counts: ";
					temp += String(true_count+false_count);
					temp += " , true: ";
					temp += String(true_count);
					temp += " , false: ";
					temp += String(false_count);
					String temp2;
					temp2 = "AUC: ";
					temp2 += ROC.AUC(); 
					file.push_back("Unique Peptides with highest score: ");
					file.push_back(temp);
					file.push_back(temp2);
					
					String temp21(ROC.cutoffPos(0.95));
					file.push_back(temp21);
					
					String temp3( ROC.cutoffNeg(0.95));
					file.push_back(temp3);
					file.push_back(" ");
			
			vector<pair<double,double> > curve = ROC.curve(curve_resolution);
			
			for(vector<pair<double,double> >::iterator it = curve.begin(); it < curve.end(); ++it)
			{
				String temp1(it->first);
				String temp2(it->second);
				temp1.append("\t");
				temp1.append(temp2);
				file.push_back(temp1);
				//file.push_back(temp2);
			}

			file.store(out);
			cout<<"Unique Peptides with highest score: "<<endl;
			cout<<"All counts: "<<true_count+false_count<<", true: "<<true_count<<", false: "<<false_count<<endl;
			cout<<"AUC: "<<ROC.AUC()<<endl;
			
			ROCCurve ROC1;
			true_count = 0;
			false_count = 0;
		//	cout<<"False:"<<endl;
			for(vector<PeptideIdentification>::iterator idxmliter = idxml_peptideIdentification.begin(); idxmliter != idxml_peptideIdentification.end(); ++idxmliter)
			{
				bool exists  = false;
				for(vector<FASTAFile::FASTAEntry>::iterator fastaiter = fasta_entries.begin(); fastaiter < fasta_entries.end();++fastaiter)
				{
					if(fastaiter->sequence.hasSubstring(idxmliter->getHits()[0].getSequence().toUnmodifiedString()))
					{
						exists  = true;
						continue;
					}
				}
				if(exists)
				{
					ROC1.insertPair(idxmliter->getHits()[0].getScore(),true);
					++true_count;
				}
				else
				{
					ROC1.insertPair(idxmliter->getHits()[0].getScore(),false);
			//		cout<<idxmliter->first<<"\n";
					++false_count;
				}
			}
					
					String temp5;
					temp5 = "All counts: ";
					temp5 += String(true_count+false_count);
					temp5 += " , true: ";
					temp5 += String(true_count);
					temp5 += " , false: ";
					temp5 += String(false_count);
					String temp6;
					temp6 = "AUC: ";
					temp6 += ROC1.AUC(); 
					file.push_back("\nUnique Peptides with highest score: ");
					file.push_back(temp5);
					file.push_back(temp6);
					
					String temp8(ROC1.cutoffPos(0.95));
					file.push_back(temp8);
					
					String temp9( ROC1.cutoffNeg(0.95));
					file.push_back(temp9);
					file.push_back(" ");
			
			vector<pair<double,double> > curve2 = ROC1.curve(curve_resolution);
			
			for(vector<pair<double,double> >::iterator it = curve2.begin(); it < curve2.end(); ++it)
			{
				String temp1(it->first);
				String temp2(it->second);
				temp1.append("\t");
				temp1.append(temp2);
				file.push_back(temp1);
				//file.push_back(temp2);
			}

			file.store(out);			
			cout<<"All top hits are count:"<<endl;
			cout<<"All counts: "<<true_count+false_count<<", true: "<<true_count<<", false: "<<false_count<<endl;
			cout<<"AUC: "<<ROC1.AUC()<<endl;
			
			
			return EXECUTION_OK;
		}
	};			
			
int main( int argc, const char** argv )
{
    TOPPDataForROCCurve tool;
    return tool.main(argc,argv);
}

/// @endcond
