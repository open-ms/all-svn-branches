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

  Struct3d::Struct3d(double mzValue, double rtValue, double intensityValue)
    : mz(mzValue), rt(rtValue), intensity(intensityValue)
  {
  }

  Struct3d::Struct3d(const Struct3d& value)
  {
    if (&value == this)
       return;
    mz = value.mz;
    rt = value.rt;
    intensity = value.intensity;
  }
  
  Struct3d::Struct3d(const Struct3d& p1, const Struct3d& p2)
    : mz(p2.mz-p1.mz), rt(p2.rt-p1.rt), intensity(p2.intensity-p1.intensity)        
  {
  }
  
  const Struct3d& Struct3d::operator=(const Struct3d& value)
  {
    if (&value == this)
      return *this;
    mz = value.mz;
    rt = value.rt;
    intensity = value.intensity;
    return *this;
  }

  Struct3d& Struct3d::operator+=(Struct3d value)
  {
    mz += value.mz;
    rt += value.rt;
    intensity += value.intensity;
    return *this;
  }

  Struct3d& Struct3d::operator-=(Struct3d value)
  {
    mz -= value.mz;
    rt -= value.rt;
    intensity -= value.intensity;
    return *this;
  }
  
  Struct3d& Struct3d::operator*=(double value)
  {
    mz *= value;
    rt *= value;
    intensity *= value;
    return *this;
  }
  
  Struct3d& Struct3d::operator/=(double value)
  {
    if(!approximatelyNul(value))
    {
      mz /= value;
      rt /= value;
      intensity /= value;
    }
    else
      *this = Struct3d();
    return *this;
  }
  
  Struct3d& Struct3d::operator*=(Struct3d value)
  {
    mz *= value.mz;
    rt *= value.rt;
    intensity *= value.intensity;
    return *this;
  }

  bool Struct3d::operator==(Struct3d value) const
  {
    return approximatelyEqual(mz, value.mz) && approximatelyEqual(rt, value.rt) && approximatelyEqual(intensity,value.intensity);
  }
  
  bool Struct3d::operator!=(Struct3d value) const
  {
    return !approximatelyEqual(mz, value.mz) || !approximatelyEqual(rt, value.rt) || !approximatelyEqual(intensity,value.intensity);
  }

  double Struct3d::length() const
  {
    double l = mz*mz + rt*rt + intensity*intensity;
    return (approximatelyNul(l)) ? 0.0 : sqrt(l);
  }
  	
  void Struct3d::normalize()
  {
    double l = length();
    if(!approximatelyNul(l))
    {
      mz /= l;
      rt /= l;
      intensity /= l;
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
    glTranslatef(mz, rt, intensity);
  }
  
  void Struct3d::glRotate(const double phi) const
  {
    glRotatef(phi, mz, rt, intensity);
  }
    
} // namespace OpenMS

