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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
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

Test

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

	/**
	 * computes simple stats (mean, stddev)
	 */
	pair<double, double> simpleStats(const vector<double>& x)
	{
		double sum_of_squares = 0.0;
		double sum = 0.0;
		UInt N = x.size();
	
		for (UInt i = 0; i < N; ++i)
		{
			sum += x[i];
			sum_of_squares += x[i] * x[i];
		}
	
		double mean = sum / (double)N;
		double stddev = 1./(double)N * sqrt(N * sum_of_squares - sum * sum);
	
		return make_pair(mean, stddev);
	}

	pair<Matrix<double>, Matrix<double> > computeMap2MapCorrelations(const ConsensusMap& map)
	{
		UInt number_of_features = map.size();
		UInt number_of_maps = map.getFileDescriptions().size();
		vector<vector<double> > feature_int(number_of_maps);
		for (UInt i = 0; i < number_of_maps; i++)
		{
			feature_int[i].resize(number_of_features);
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
		// info: mean, stddev for intensities of one map
		for (UInt i = 0; i < number_of_maps; ++i)
		{
			pair<double, double> stats = simpleStats(feature_int[i]);
			cout << "stats for map " << i << ": " << stats.first << " +/- " << stats.second << endl;
		}
	
		Matrix<double> ratio_matrix(number_of_maps, number_of_maps);
		Matrix<double> corr_matrix(number_of_maps, number_of_maps);
		//vector<double> zero_ratios(number_of_maps);
		for (UInt i = 0; i < number_of_maps; ++i)
		{
			for (UInt j = 0; j < number_of_maps; j++)
			{
				double corr = 0.0;
				double ss_x = 0.0;
				double ss_y = 0.0;
				vector<double> ratios;
				for (UInt k = 0; k < number_of_features; ++k)
				{
					corr += feature_int[i][k] * feature_int[j][k];
					ss_x += feature_int[i][k] * feature_int[i][k];
					ss_y += feature_int[j][k] * feature_int[j][k];
					if (feature_int[i][k] != 0.0 && feature_int[j][k] != 0.0)
					{	
						double ratio = feature_int[i][k] / feature_int[j][k];
						if (ratio > 0.67 && ratio < 1.5)
						{
							ratios.push_back(ratio);	
						}
					}
				}
				corr /= sqrt(ss_x * ss_y);
				pair<double, double> ratio = simpleStats(ratios);
				cout << "normalized correlation between " << i << " and " << j << " : " << corr << " ratio: " << ratio.first << " +/- " << ratio.second << " for " << ratios.size() << " shared pairs." << endl;
				ratio_matrix(i,j) = ratio.first;
				corr_matrix(i,j) = corr;
			}
		}
		return make_pair(ratio_matrix, corr_matrix);
	}

	void normalizeMaps(ConsensusMap& map, const vector<double>& ratios)
	{
		ConsensusMap::Iterator cf_it;
		UInt idx = 0;
		for (cf_it = map.begin(); cf_it != map.end(); ++cf_it, ++idx)
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
		cout << "Loading consensus map " << in << endl;
		infile.load(in, map);

		map.sortBySize();

		// info
		cout << "Number of consensus features: " << map.size() << endl;
		cout << "Number of maps: " << map.getFileDescriptions().size() << endl;

		//map normalization
		pair<Matrix<double>, Matrix<double> > results = computeMap2MapCorrelations(map);
		normalizeMaps(map, results.first.row(0));
		cout << "ratio matrix: " << endl << results.first << endl << endl;
		cout << "corr matrix: " << endl << results.second << endl;

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
