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

#ifndef OPENMS_VISUAL_ARROW3D_H
#define OPENMS_VISUAL_ARROW3D_H


#include <QtOpenGL/QGLWidget>
#include <OpenMS/MATH/MISC/Struct3d.h>
#include <OpenMS/VISUAL/ColorRGBA.h>

namespace OpenMS
{
	/**
		@brief OpenGL Canvas for 3D-visualization of map data
		
		@note Do not use this class directly. Use Spectrum3DCanvas instead!
		
		@ingroup SpectrumWidgets
	*/
	
  class Arrow3d
  {
    public:
	    Arrow3d();
	    ~Arrow3d();

      void draw(const Struct3d& pos, const Struct3d& dir, const double sizeFactor = 1.0) const;
      void setColor(ColorRGBA color);
      void setColor(double rr, double gg, double bb, double aa);
      void setColor(QColor color);

    private:
	    double rotate_(Struct3d& axis, const Struct3d& vector) const;
	    
	    GLUquadricObj *hat;
	    GLUquadricObj *disk;
	    GLUquadricObj *base;
	    GLUquadricObj *bottom;

      ColorRGBA color_;     
  };
}
#endif
