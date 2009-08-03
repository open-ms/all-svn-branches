// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2009
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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/SettingsParser.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FileSettings FileSettings

	@brief Add meta infos to data.

	This tool can show basic information about the data in several peak, feature and consensus feature files. It can
	- show information about the data range of a file (m/z, RT, intensity)
	- show a statistical summary for intensities and qualities
	- show an overview of the metadata
	- validate several XML formats against their XML schema
	- check for corrupt data in a file (e.g. dupliacte spectra)

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileInfo.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileSettings
	: public TOPPBase
{
 public:
	TOPPFileSettings()
		: TOPPBase("FileSettings","Add meta infos to data.")
	{
	}

 protected:
	virtual void registerOptionsAndFlags_()
	{
  	registerInputFile_("in","<file>","","input raw data file ");
		setValidFormats_("in", StringList::create("mzML"));
		registerOutputFile_("out","<file>","","output raw data file ");
  	setValidFormats_("out", StringList::create("mzML"));

  	registerInputFile_("settings_file","<file>","","input settings text file ");
		setValidFormats_("settings_file", StringList::create("txt"));				
		registerStringOption_("settings_command","<command>","","input settings command line.",false);
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");		
		String out = getStringOption_("out");
		String settingsFile = getStringOption_("settings_file");		
		String settingsCommand = getStringOption_("settings_command");
		
		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------
		MzMLFile file;
		MSExperiment<Peak1D> exp;
		file.setLogType(log_type_);
		file.load(in, exp);
		
		// add meta infos from text file
		if(settingsFile != "")
      SettingsParser().settingsFile(settingsFile, exp);

    // add meta infos from command line
    if(settingsCommand != "")
      SettingsParser().settingsCommand(settingsCommand, exp);
		
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		file.store(out, exp);

		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPFileSettings tool;
	return tool.main(argc,argv);
}

/// @endcond
