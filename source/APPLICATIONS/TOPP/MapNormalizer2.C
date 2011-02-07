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
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MapNormalizer2 MapNormalizer2

	@brief Normalizes maps of one consensusXML file.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapNormalizer2
  : public TOPPBase
{

public:
	TOPPMapNormalizer2()
		: TOPPBase("MapNormalizer2","Normalizes maps of one consensusXML file")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create("consensusXML"));
		registerOutputFile_("out","<file>","","Output file",true);
		setValidFormats_("out",StringList::create("consensusXML"));
	}

	double mean(const vector<double>& x)
	{
		double sum = 0.0;
		UInt N = x.size();
	
		for (UInt i = 0; i < N; ++i)
		{
			sum += x[i];
		}
		double mean = sum / (double)N;
	
		return mean;
	}

	vector<double> computeCorrelation(const ConsensusMap& map)
	{
		UInt number_of_features = map.size();
		UInt number_of_maps = map.getFileDescriptions().size();
		vector<vector<double> > feature_int(number_of_maps);
		UInt map_with_most_features = 0;
		for (UInt i = 0; i < number_of_maps; i++)
		{
			feature_int[i].resize(number_of_features);
			if (map.getFileDescriptions()[i].size > map.getFileDescriptions()[map_with_most_features].size)
			{
				map_with_most_features = i;
			}
		}
		
		ConsensusMap::ConstIterator cf_it;
		UInt idx = 0;
		for (cf_it = map.begin(); cf_it != map.end(); ++cf_it, ++idx)
		{
			ConsensusFeature::HandleSetType::const_iterator f_it;
			for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
			{
				feature_int[f_it->getMapIndex()][idx] = f_it->getIntensity();
			}
		}

		vector<double> ratio_vector(number_of_maps);
		for (UInt j = 0; j < number_of_maps; j++)
		{
			vector<double> ratios;
			for (UInt k = 0; k < number_of_features; ++k)
			{
				if (feature_int[map_with_most_features][k] != 0.0 && feature_int[j][k] != 0.0)
				{	
					double ratio = feature_int[map_with_most_features][k] / feature_int[j][k];
					//TODO ratio als log Parameter
					if (ratio > 0.67 && ratio < 1.5)
					{
						ratios.push_back(ratio);	
					}
				}
			}
			ratio_vector[j] = mean(ratios);
		}
		return ratio_vector;
	}

	void normalizeMaps(ConsensusMap& map, const vector<double>& ratios)
	{
		ConsensusMap::Iterator cf_it;
		for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
		{
			ConsensusFeature::HandleSetType::iterator f_it;
			for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
			{	
				f_it->asMutable().setIntensity(f_it->getIntensity() * ratios[f_it->getMapIndex()]);
			}
		}
	}

	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in");

		ConsensusXMLFile infile;
		ConsensusMap map;
		infile.load(in, map);

		map.sortBySize();

		//map normalization
		vector<double> results = computeCorrelation(map);
		normalizeMaps(map, results);

		String out = getStringOption_("out");
		infile.store(out,map);

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPMapNormalizer2 tool;
  return tool.main(argc,argv);
}

/// @endcond
