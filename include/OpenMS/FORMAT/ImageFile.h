// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
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

#ifndef OPENMS_FORMAT_IMAGEFILE_H
#define OPENMS_FORMAT_IMAGEFILE_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <QImage>

namespace OpenMS
{
 	/**
 		@brief File adapter for images files.

 		Open image as map. X coords become MZ and Y coords become RT. Intensity was calculed as sum of RGB channels.<br />
 		Supported formats are : BMP (Windows Bitmap), GIF (Graphic Interchange Format), JPG and JPEG (Joint Photographic Experts Group), PNG (Portable Network Graphics), PBM (Portable Bitmap), PGM (Portable Graymap), PPM (Portable Pixmap), TIFF (Tagged Image File Format), XBM (X11 Bitmap), XPM (X11 Pixmap).
  	
  	@ingroup FileIO
  */
  
  class OPENMS_DLLAPI ImageFile
  {
    public:
      /// Default constructor
      ImageFile();

			/**
				@brief Loads a map from a XMass file.

				@p map has to be a MSSpectrum or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be opened
			*/      #include <QImage>
      template <class MapType>
      void load(const String& filename, MapType& map)
      {        
		    //try to open file
		    QImage img(filename.toQString());
		    if(img.isNull())
		    {
			    throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		    }

		    int w = img.width();
		    int h = img.height();
		    
		    map.reset();
		    map.clear();
				map.reserve(h);
        typename MapType::SpectrumType spectrum;
          
        for(int y=0; y<h; ++y)
        {
          spectrum.clear();
          spectrum.reserve(w);
          spectrum.setRT(1.0 + 1.0 * y);
          
          for(int x=0; x<w; ++x)
          {
            typename MapType::SpectrumType::PeakType p;
					  p.setIntensity(256.0 - 1.0 * qGray(img.pixel(x, h-y-1)));
					  p.setMZ(1.0 + 1.0 * x);
					  spectrum.push_back(p);            
          }   
          map.push_back(spectrum);
        }
      }

			/**
				@brief Stores a map in a image file (not avaible)

				@exception Exception::FileNotWritable is thrown
			*/
      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Images files not writable."));
      }
      
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_IMAGEFILE_H

