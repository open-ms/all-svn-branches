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

#ifndef OPENMS_FORMAT_FILEHANDLERTHREAD_H
#define OPENMS_FORMAT_FILEHANDLERTHREAD_H

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <QThread>

namespace OpenMS
{
	/**
		@brief Facilitates file handling by file type recognition.

		This class provides file type recognition from the file name and
		from the file content.

		It also offer a common interface to load MSExperiment data
		and allows querying for supported file types.

		@see FileTypes

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI FileHandlerThread
    : public QThread  
  {
    Q_OBJECT
	  
	  public:
    	///@name Type definitions
    	//@{
    	//Feature map type
    	typedef LayerData::FeatureMapType FeatureMapType;
    	//Consensus feature map type
    	typedef LayerData::ConsensusMapType ConsensusMapType;
    	//Peak map type
    	typedef LayerData::ExperimentType ExperimentType;
    	//@}
    		  
	    FileHandlerThread(
	      String filename,
	      bool show_options, 
	      bool add_to_recent, 
	      String caption, 
	      UInt window_id, 
	      Size spectrum_id);
	    void run();

    signals:
      void fileLoaded(bool success = false);	    
    
    public:
      String filename_;
      String abs_filename_;
      bool show_options_;
      bool add_to_recent_;
      String caption_; 
      UInt window_id_;
      Size spectrum_id_;
      bool is_2D_;
      bool is_feature_;      
		  FeatureMapType feature_map_;
		  ExperimentType peak_map_;
		  ConsensusMapType consensus_map_;
		  String msg_;
	};

} //namespace

#endif //OPENMS_FORMAT_FILEHANDLERTHREAD_H
