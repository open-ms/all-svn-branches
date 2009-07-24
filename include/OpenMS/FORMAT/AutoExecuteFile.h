// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    AutoExecute file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_AUTOEXECUTEFILE_H
#define OPENMS_FORMAT_AUTOEXECUTEFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/FORMAT/XMassFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <fstream>
using namespace std;

namespace OpenMS
{
 	/**
 		@brief File adapter for autoExecute file or list file.
  	
    For example, to merge fid files "find source | grep fid > destination"
  
  	@ingroup FileIO
  */
  
  class OPENMS_DLLAPI AutoExecuteFile
  {
    public:
      /// Default constructor
      AutoExecuteFile();
      ~AutoExecuteFile();     
      
      template <typename MapType>
		  void load(const String& filename, MapType& map)
      {
        StringList file_list;
        try
				{
					file_list = getFileList(filename);
cout << "list size : " << file_list.size() << endl;					
				}
				catch(Exception::FileNotFound)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("'AutoExecute' file not writable."));
				}
				map.reset();
				map.resize(file_list.size());
				for(Size index = 0; index < file_list.size(); ++index)
				{
				  typename MapType::SpectrumType spectrum;
					XMassFile().load(file_list[index], spectrum);
          spectrum.setRT(index);
          spectrum.setMSLevel(1);					
				  map.push_back(spectrum);
				}
				XMassFile().importExperimentalSettings(file_list[0], map);
      }

      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("'AutoExecute' file not writable."));
      }      
      
 	      
	  private:
						
      /// Import file list
      StringList getFileList(
	      const String & filename,  
	      const bool isAutoExecute = false);
	    String autoExecuteToFilename(const String & line);
      void readAutoExecuteHeader(std::ifstream & is_);

	    unsigned int PosOnScout_; // position on chip, format [A-P]:[1-24]
	    unsigned int SpectrumDirectory_; // spectrum directory
	    unsigned int SpectrumFilename_; // spectrum filename
	    unsigned int ChipOnScout_; // normal or calibrant position, "0" or "1"
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_AUTOEXECUTEFILE_H

