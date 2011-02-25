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

This tool filters featureXML files by their annotation. You can choose the following options:
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
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "input file");
		setValidFormats_("in", StringList::create("featureXML"));
		registerOutputFile_("out", "<file>", "", "output file");
		setValidFormats_("out", StringList::create("featureXML"));
		addEmptyLine_();
		registerFlag_("unassigned","keep unassigned Peptide Identifications");
		registerFlag_("annotated","filter features with annotations");
		registerFlag_("not_annotated","filter features without annotations");
		registerStringList_("sequence","<sequence>",StringList(),"sequences to filter, i.e. Oxidation or LYSNLVER", false);
	}

	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in");

		FeatureXMLFile infile;
		FeatureMap<> map;
		infile.load(in, map);
		bool unassigned = getFlag_("unassigned");
		bool annotated = getFlag_("annotated");
		bool not_annotated = getFlag_("not_annotated");
		StringList sequence = getStringList_("sequence");

		//copy all properties
		FeatureMap<> new_map = map;
		//but delete feature information
		new_map.clear(false);
		//loop over all features
		for (FeatureMap<>::ConstIterator fm_it = map.begin(); fm_it != map.end(); ++fm_it)
		{
			//flag: annotated and non-empty peptideIdentifications
			if (annotated && !fm_it->getPeptideIdentifications().empty())
			{
				if (sequence.size() > 0)
				{
					//loop over all peptideIdentifications
					bool seq_found = false;
					for (vector<PeptideIdentification>::const_iterator pep_id_it = fm_it->getPeptideIdentifications().begin(); pep_id_it != fm_it->getPeptideIdentifications().end(); ++pep_id_it)
					{
						if (seq_found) break;
						//loop over all peptideHits
						for (vector<PeptideHit>::const_iterator pep_hit_it = pep_id_it->getHits().begin(); pep_hit_it != pep_id_it->getHits().end(); ++pep_hit_it)
						{
							if (seq_found) break;
							//loop over all sequence entries of the StringList
							for (StringList::ConstIterator seq_it = sequence.begin(); seq_it != sequence.end(); ++seq_it)
							{
								if (seq_found) break;
								if (pep_hit_it->getSequence().toString().hasSubstring(*seq_it)
									|| pep_hit_it->getSequence().toUnmodifiedString().hasSubstring(*seq_it))
								{
									new_map.push_back(*fm_it);
									seq_found = true;
								}
							}
						}
					}
				}else
				{
					new_map.push_back(*fm_it);
				}
			}
			//flag: not_annotated and no peptideIdentifications
			if (not_annotated && fm_it->getPeptideIdentifications().empty())
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

		String out = getStringOption_("out");
		infile.store(out, new_map);

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPFileFilter2 tool;
  return tool.main(argc,argv);
}

/// @endcond
