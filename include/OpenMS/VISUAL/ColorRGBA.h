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

#ifndef OPENMS_VISUAL_COLORRGBA_H
#define OPENMS_VISUAL_COLORRGBA_H

#include <QColor>
    
namespace OpenMS
{
	/**
		@brief OpenGL Canvas for 3D-visualization of map data
		
		@note Do not use this class directly. Use Spectrum3DCanvas instead!
		
		@ingroup SpectrumWidgets
	*/
  
  class ColorRGBA
  {
    public:
	    ColorRGBA();
	    ColorRGBA(double rr, double gg, double bb, double aa = 1.0);
	    ColorRGBA(QColor color);
	    ColorRGBA(const ColorRGBA& color);
	    
      QColor toQColor() const;
      void glColor() const;
      
    private:
      double r, g, b, a;
  };
  
}
#endif
