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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandlerThread.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

namespace OpenMS
{
  FileHandlerThread::FileHandlerThread(
      String filename,
      bool show_options, 
      bool add_to_recent, 
      String caption, 
      UInt window_id, 
      Size spectrum_id)
  : QThread(),
    filename_(filename),
    abs_filename_(""),
    show_options_(show_options),
    add_to_recent_(add_to_recent),
    caption_(caption),    
    window_id_(window_id),
    spectrum_id_(spectrum_id),
    is_2D_(false),
		is_feature_(false),
		msg_("")
  {
  }
  
  void FileHandlerThread::run()
  {
    abs_filename_ = File::absolutePath(filename_);

  	//check if the file exists
    if (!File::exists(abs_filename_))
    {
    	msg_ = String("The file '") + abs_filename_ + String("' does not exist!");
      emit(fileLoaded());    	
      return;
    }

		//determine file type		    
  	FileHandler fh;
		FileTypes::Type file_type = fh.getType(abs_filename_);

		if(FileTypes::UNKNOWN == file_type)
		{
			msg_ = String("Could not determine file type of '") + abs_filename_ + String("'!");
      return;
		}
		
		//abort if file type unsupported
		if(FileTypes::INI==file_type || FileTypes::IDXML==file_type)
		{
			msg_ = String("The type '") + fh.typeToName(file_type) + String("' is not supported!");
      emit(fileLoaded());			
      return;
		}

		//try to load data and determine if it is 1D or 2D data
    try
    {
	    if(FileTypes::FEATUREXML == file_type)
	    {
        FeatureXMLFile().load(abs_filename_, feature_map_);
        is_2D_ = true;
        is_feature_ = true;
      }
      else if(FileTypes::CONSENSUSXML == file_type)
	    {
        ConsensusXMLFile().load(abs_filename_, consensus_map_);
        is_2D_ = true;
        is_feature_ = true;
      }
      else
      {     
      	fh.loadExperiment(abs_filename_, peak_map_, file_type);
      	UInt ms1_scans = 0;
      	for(Size iSpectrum=0; iSpectrum<peak_map_.size(); ++iSpectrum)
      	{
      		if(peak_map_[iSpectrum].getMSLevel() == 1) 
      		  ++ms1_scans;
      		if(ms1_scans > 1)
      		{
      			is_2D_ = true;
      			break;
      		}
      	}
      	if (peak_map_.getChromatograms().size() > 1) 
      	  is_2D_ = true;
      }
    }
    catch(Exception::BaseException& e)
    {
    	msg_ = e.what();
    	emit(fileLoaded());
      return;
    }
    
    //try to add the data
		if("" == caption_)
		{
			caption_ = File::basename(abs_filename_);
    }
    else
    {
    	abs_filename_ = "";
    }
    
    emit(fileLoaded(true));
  }
} // namespace OpenMS

