// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FEATUREXMLFILE_H
#define OPENMS_FORMAT_FEATUREXMLFILE_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS
{
	/**
  	@brief This class provides Input/Output functionality for feature maps

		A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/. 
		
  	@note This format will eventually be replaced by the HUPO-PSI AnalysisXML format!
  	
  	@note If a feature meta value named 'id' is present, it is stored in the feature 'id' attribute in the file.
  	      If meta value 'id' is not set, incrementing numbers are used.

  	@ingroup FileIO
  */
  class OPENMS_DLLAPI FeatureXMLFile
  	: public Internal::XMLFile
  {
	
		public:
		
			/** @name Constructors and Destructor */
			//@{
			
			///Default constructor
			FeatureXMLFile();
			///Destructor
			~FeatureXMLFile();
			//@}
			
			/**
				@brief loads the file with name @p filename into @p map.
			
				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(String filename, FeatureMap<>& feature_map);
					
			/**
				@brief stores the map @p feature_map in file with name @p filename.
				
				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(String filename, const FeatureMap<>& feature_map) const;
			
      /// Mutable access to the options for loading/storing 
      PeakFileOptions& getOptions();

      /// Non-mutable access to the options for loading/storing 
      const PeakFileOptions& getOptions() const;
		
		protected:
		
			/// options for reading / writing
			PeakFileOptions options_;

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_FeatureXMLFile_H
