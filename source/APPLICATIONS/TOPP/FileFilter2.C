// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FileFilter2 FileFilter2

	@brief Extracts portions of the data from featureXML files by their annotation.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileFilter \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinder </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool that profits @n on annotated features </td>
		</tr>

	</table>
</CENTER>

This tool filters featureXML or consensusXML files by their annotation. You can choose the following options:
	- filter unassigned peptide identifications
	- filter annotated features
	- filter non-annotated features
	- filter annotated features by sequence



<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileFilter2.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileFilter2
  : public TOPPBase
{

public:
	TOPPFileFilter2()
		: TOPPBase("FileFilter2","Filters featureXML files by ids.")
	{
	}

protected:
	static bool checkPeptideIdentification(const BaseFeature& feature, const bool annotated, const bool not_annotated, const StringList& sequence)
	{
		//flag: annotated and non-empty peptideIdentifications
		if (annotated && !feature.getPeptideIdentifications().empty())
		{
			if (sequence.size() > 0)
			{
				//loop over all peptideIdentifications
				for (vector<PeptideIdentification>::const_iterator pep_id_it = feature.getPeptideIdentifications().begin(); pep_id_it != feature.getPeptideIdentifications().end(); ++pep_id_it)
				{
					//loop over all peptideHits
					for (vector<PeptideHit>::const_iterator pep_hit_it = pep_id_it->getHits().begin(); pep_hit_it != pep_id_it->getHits().end(); ++pep_hit_it)
					{
						//loop over all sequence entries of the StringList
						for (StringList::ConstIterator seq_it = sequence.begin(); seq_it != sequence.end(); ++seq_it)
						{
							if (pep_hit_it->getSequence().toString().hasSubstring(*seq_it)
								|| pep_hit_it->getSequence().toUnmodifiedString().hasSubstring(*seq_it))
							{
								return true;
							}
						}
					}
				}
			}else
			{
				return true;
			}
		}
		//flag: not_annotated and no peptideIdentifications
		if (not_annotated && feature.getPeptideIdentifications().empty())
		{
			return true;
		}
		return false;
	}

	void registerOptionsAndFlags_()
	{
		String formats("featureXML,consensusXML");

		registerInputFile_("in", "<file>", "", "input file");
		setValidFormats_("in", StringList::create(formats));

		registerOutputFile_("out", "<file>", "", "output file");
		setValidFormats_("out", StringList::create(formats));

		addEmptyLine_();
		registerFlag_("unassigned","keep unassigned Peptide Identifications");
		registerFlag_("annotated","filter features with annotations");
		registerFlag_("not_annotated","filter features without annotations");
		registerStringList_("sequence","<sequence>",StringList(),"sequences to filter, i.e. Oxidation or LYSNLVER", false);
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		FileHandler fh;

		//input file name and type
		String in = getStringOption_("in");
		FileTypes::Type in_type = fh.getType(in);

		//output file name and type
		String out = getStringOption_("out");

		//id-filtering parameters
		bool unassigned = getFlag_("unassigned");
		bool annotated = getFlag_("annotated");
		bool not_annotated = getFlag_("not_annotated");
		StringList sequence = getStringList_("sequence");

		if (in_type == FileTypes::FEATUREXML)
		{
			FeatureXMLFile f;
			FeatureMap<> feature_map;
			f.load(in, feature_map);
	
			//copy all properties
			FeatureMap<> new_map = feature_map;
			//but delete feature information
			new_map.clear(false);
			//loop over all features
			for (FeatureMap<>::ConstIterator fm_it = feature_map.begin(); fm_it != feature_map.end(); ++fm_it)
			{
				if (checkPeptideIdentification(*fm_it, annotated, not_annotated, sequence))
				{
					new_map.push_back(*fm_it);
				}
			}
			//delete unassignedPeptideIdentifications
			if (!unassigned)
			{
				new_map.getUnassignedPeptideIdentifications().clear();
			}
			//update minimum and maximum position/intensity
			new_map.updateRanges();

			//annotate output with data processing info
			addDataProcessing_(new_map, getProcessingInfo_(DataProcessing::FILTERING));

			f.store(out, new_map);

		}else if (in_type == FileTypes::CONSENSUSXML)
		{
			ConsensusXMLFile f;
			ConsensusMap consensus_map;
			f.load(in, consensus_map);
			
			//copy all properties
			ConsensusMap new_map = consensus_map;
			//but delete feature information
			new_map.clear(false);
			//loop over all features
			for (ConsensusMap::ConstIterator cm_it = consensus_map.begin(); cm_it != consensus_map.end(); ++cm_it)
			{
				if (checkPeptideIdentification(*cm_it, annotated, not_annotated, sequence))
				{
					new_map.push_back(*cm_it);
				}
			}
			//delete unassignedPeptideIdentifications
			if (!unassigned)
			{
				new_map.getUnassignedPeptideIdentifications().clear();
			}
			//update minimum and maximum position/intensity
			new_map.updateRanges();

			//annotate output with data processing info
			addDataProcessing_(new_map, getProcessingInfo_(DataProcessing::FILTERING));

			f.store(out, new_map);
		}else
		{
			writeLog_("Unknown input file type given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPFileFilter2 tool;
  return tool.main(argc,argv);
}

/// @endcond
