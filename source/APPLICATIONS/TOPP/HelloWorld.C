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
// $Maintainer: Holger Plattfaut $
// $Authors: Holger Plattfaut $
// --------------------------------------------------------------------------
// #include <OpenMS/FORMAT/ConsensusXMLFile.h>
// #include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
// #include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_HelloWorld HelloWorld

	@brief Test.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPHelloWorld
  : public TOPPBase
{

public:
	TOPPHelloWorld()
		: TOPPBase("HalloWelt","Test.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","Input file",true);
		setValidFormats_("in",StringList::create("mzML"));
		registerOutputFile_("out","<file>","","Output file",true);
		setValidFormats_("out",StringList::create("mzML"));
		registerStringOption_("string","<name>","","String to print",true);
		// setValidStrings_("type", getToolList()[toolName_()]);

		// registerSubsection_("algorithm","Algorithm parameters section");
	}

	// Param getSubsectionDefaults_(const String& /*section*/) const
	/*{
		String type = getStringOption_("type");
    FeatureGroupingAlgorithm* algo = Factory<FeatureGroupingAlgorithm>::create(type);
		Param p = algo->getParameters();
    delete algo;
    return p;
	}*/

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");

		String out = getStringOption_("out");

		String str = getStringOption_("string");

		MSExperiment<Peak1D> exp;		
		MzMLFile f;
		f.load(in,exp);
		f.store(out,exp);

		cout << str << endl;

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPHelloWorld tool;
  return tool.main(argc,argv);
}

/// @endcond
