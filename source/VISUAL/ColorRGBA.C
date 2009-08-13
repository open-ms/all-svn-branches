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

#include <OpenMS/VISUAL/ColorRGBA.h>
#include <QtOpenGL/QGLWidget>

namespace OpenMS
{
	ColorRGBA::ColorRGBA()
		: r(0.0), g(0.0), b(0.0), a(1.0)
	{
	}
	
	ColorRGBA::ColorRGBA(double rr, double gg, double bb, double aa)
		: r(rr), g(gg), b(bb), a(aa)
	{
	}

  ColorRGBA::ColorRGBA(QColor color)
   : r(color.redF()), g(color.greenF()), b(color.blueF()), a(color.alphaF())
  {	
  }
  	
	ColorRGBA::ColorRGBA(const ColorRGBA& color)
	{
    if (&color == this)
       return;  
    r = color.r;
    g = color.g;
    b = color.b;	  
    a = color.a;
	}
	
	QColor ColorRGBA::toQColor() const
  {
    QColor color;
    color.setRedF(r);
    color.setGreenF(g);
    color.setBlueF(b);
    color.setAlphaF(a);
	  return color;	
  }
  
  void ColorRGBA::glColor() const
  {
    glColor4d(r, g, b, a);
  }
	
}//Namespace

