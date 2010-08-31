// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/MISC/CrossCorrelationCalculator1.h>

using namespace OpenMS;

class TOPPAutocorrelation3
	: public TOPPBase
{
 public:
	TOPPAutocorrelation3()
		: TOPPBase("Autocorrelation3","Calculates the autocorrelation within given spectra.")
	{

	}

 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file");
		setValidFormats_("in",StringList::create("mzML"));
		registerOutputFile_("out","<file>","","output file");
		setValidFormats_("out",StringList::create("mzML"));
		registerDoubleOption_("stepwidth","<stepwidth>",0.001,"Stepwidth");
		registerIntOption_("spectrum_selection","<spectrum_selection>",0,"Spectrum selection for exact positions");
		registerDoubleList_("estimated_positions","<estimated_positions>",DoubleList::create(0.5),"Estimated positions. Autocorrelation2 prints the exact positions of these values");
		registerDoubleOption_("tolerance","<tolerance>",0.0,"Maximal possible deviation from estimated positions");
		registerDoubleOption_("gauss_mean","<gauss_mean>",0.0,"Mean position of the gaussian curve");
		registerDoubleOption_("gauss_width","<gauss_width>",1.0,"Width of the gaussian curve");
	}

	ExitCodes main_(int , const char**)
		{
			//input file names and types
			String in = getStringOption_("in");
			String out = getStringOption_("out");

			String debug_trunk = in;
			if (in.has('.'))
			{
				debug_trunk = in.substr(0,in.find_first_of('.'));
			}

			MzMLFile file;
			MSExperiment<Peak1D> exp_in;
			file.setLogType(log_type_);
			file.load(in,exp_in);
			exp_in.updateRanges();

			MSExperiment<Peak1D> exp_out;

			DoubleReal tolerance=getParam_().getValue("tolerance");
			DoubleReal gauss_mean=getParam_().getValue("gauss_mean");
			DoubleReal gauss_sigma=getParam_().getValue("gauss_width");
			DoubleReal stepwidth=getParam_().getValue("stepwidth");
			Size spectrum_selection=(Size) getParam_().getValue("spectrum_selection");
			DoubleList positions=getParam_().getValue("estimated_positions");
			CrossCorrelationCalculator1 ac(stepwidth,gauss_mean,gauss_sigma);
			ac.setLogType(log_type_);
			std::vector<DoubleReal> data=ac.calculate(exp_in,exp_out,spectrum_selection);
			std::vector<DoubleReal> position_vector(positions.begin(),positions.end());
			std::vector<DoubleReal> exact_positions=ac.getExactPositions(data,position_vector,tolerance);
			addDataProcessing_(exp_out, getProcessingInfo_(DataProcessing::PEAK_PICKING));
			file.store(out,exp_out);
			std::cout << "\nExact positions for spectrum " << spectrum_selection << ":\n\nExpected\tExact\n";
			for (Size i=0;i<exact_positions.size();++i)
			{
				std::cout << position_vector[i] << "\t\t" << exact_positions[i] << "\n";
			}
			std::cout << std::endl;
			return EXECUTION_OK;
		}
};

	int main( int argc, const char** argv )
	{
		TOPPAutocorrelation3 tool;
		return tool.main(argc,argv);
	}
