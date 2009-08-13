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

#include <OpenMS/MATH/MISC/Struct3d.h>
#include <QtOpenGL/QGLWidget>

namespace OpenMS
{

  Struct3d::Struct3d(double xx, double yy, double zz)
    : x(xx), y(yy), z(zz)
  {
  }

  Struct3d::Struct3d(const Struct3d& value)
  {
    if (&value == this)
       return;
    x = value.x;
    y = value.y;
    z = value.z;
  }
  
  Struct3d::Struct3d(const Struct3d& p1, const Struct3d& p2)
    : x(p2.x-p1.x), y(p2.y-p1.y), z(p2.z-p1.z)        
  {
  }
  
  const Struct3d& Struct3d::operator=(const Struct3d& value)
  {
    if (&value == this)
      return *this;
    x = value.x;
    y = value.y;
    z = value.z;
    return *this;
  }

  Struct3d& Struct3d::operator+=(Struct3d value)
  {
    x += value.x;
    y += value.y;
    z += value.z;
    return *this;
  }

  Struct3d& Struct3d::operator-=(Struct3d value)
  {
    x -= value.x;
    y -= value.y;
    z -= value.z;
    return *this;
  }
  
  Struct3d& Struct3d::operator*=(double value)
  {
    x *= value;
    y *= value;
    z *= value;
    return *this;
  }
  
  Struct3d& Struct3d::operator/=(double value)
  {
    if(!approximatelyNul(value))
    {
      x /= value;
      y /= value;
      z /= value;
    }
    else
      *this = Struct3d();
    return *this;
  }
  
  Struct3d& Struct3d::operator*=(Struct3d value)
  {
    x *= value.x;
    y *= value.y;
    z *= value.z;
    return *this;
  }

  bool Struct3d::operator==(Struct3d value) const
  {
    return approximatelyEqual(x, value.x) && approximatelyEqual(y, value.y) && approximatelyEqual(z, value.z);
  }
  
  bool Struct3d::operator!=(Struct3d value) const
  {
    return !approximatelyEqual(x, value.x) || !approximatelyEqual(y, value.y) || !approximatelyEqual(z, value.z);
  }

  double Struct3d::length() const
  {
    double l = x*x + y*y + z*z;
    return (approximatelyNul(l)) ? 0.0 : sqrt(l);
  }
  	
  void Struct3d::normalize()
  {
    double l = length();
    if(!approximatelyNul(l))
    {
      x /= l;
      y /= l;
      z /= l;
    }
    else
      *this = Struct3d();
  }

  Struct3d Struct3d::normalize() const
  {
    Struct3d value = *this;
    value.normalize();
    return value;
  }
  
  void Struct3d::glTranslate() const
  {
    glTranslatef(x, y, z);
  }
  
  void Struct3d::glRotate(const double phi) const
  {
    glRotatef(phi, x, y, z);
  }
    
} // namespace OpenMS

