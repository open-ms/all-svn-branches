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
#include <OpenMS/MATH/MISC/CrossCorrelationCalculator.h>

using namespace OpenMS;

class TOPPAutocorrelation2
	: public TOPPBase
{
 public:
	TOPPAutocorrelation2()
		: TOPPBase("Autocorrelation","Calculates the autocorrelation within given spectra.")
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
		registerIntOption_("spectrum_selection","<spectrum_selection>",0,"Spectrum to fit by a gaussian curve");
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

			DoubleReal stepwidth=getParam_().getValue("stepwidth");
			DoubleReal gauss_mean=getParam_().getValue("gauss_mean");
			DoubleReal gauss_sigma=getParam_().getValue("gauss_width");
			Size spectrum_selection=(Size) getParam_().getValue("spectrum_selection");
			CrossCorrelationCalculator ac(stepwidth,gauss_mean,gauss_sigma);
			ac.setLogType(log_type_);
			ac.calculate(exp_in,exp_out,spectrum_selection);

			addDataProcessing_(exp_out, getProcessingInfo_(DataProcessing::PEAK_PICKING));
			file.store(out,exp_out);

			return EXECUTION_OK;
		}
};

	int main( int argc, const char** argv )
	{
		TOPPAutocorrelation2 tool;
		return tool.main(argc,argv);
	}
