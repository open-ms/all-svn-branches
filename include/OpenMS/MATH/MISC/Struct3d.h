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

#ifndef OPENMS_MATH_MISC_STRUCT3D_H
#define OPENMS_MATH_MISC_STRUCT3D_H

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <vector>
    
namespace OpenMS
{
	/**
		@brief rounds @p x up to the next decimal power 10 ^ @p decPow

		@verbatim
		e.g.: (123.0 , 1)  => 130
		      (123.0 , 2)  => 200
				  (0.123 ,-2)  => 0.13 ( 10^-2 = 0.01 )
		@endverbatim

		@ingroup MathFunctionsMisc
	*/
	
	using namespace Math;
	using namespace std;
	
  class Struct3d
  {    
    public:    
      Struct3d(double xx = 0.0, double yy = 0.0, double zz = 0.0);
      Struct3d(const Struct3d& value);
      Struct3d(const Struct3d& p1, const Struct3d& p2);
      
      const Struct3d& operator=(const Struct3d& value);
      Struct3d& operator+=(Struct3d value);
      Struct3d& operator-=(Struct3d value);
      Struct3d& operator*=(double value);
      Struct3d& operator/=(double value);
      Struct3d& operator*=(Struct3d value);
      bool operator==(Struct3d value) const;
      bool operator!=(Struct3d value) const;
      
      double length() const;
      void normalize();
      Struct3d normalize() const;
      
      void glTranslate() const;
      void glRotate(const double phi) const;
        
	    double x,y,z;       
  };
  
  inline Struct3d normalCrossProduct(Struct3d const& u, Struct3d const& v)
  {
	  Struct3d n;
    n.x = u.y * v.z - u.z * v.y;
    n.y = u.z * v.x - u.x * v.z;
    n.z = u.x * v.y - u.y * v.x;	
    n.normalize();
	  return n;
  }

  inline double dotProduct(Struct3d const& u, Struct3d const& v)
  {
	  return u.x*v.x + u.y*v.y + u.z*v.z;
  }  
  
  typedef vector<Struct3d> Vector3d;
  typedef Vector3d::iterator Iterator3d;
      
} // namespace OpenMS

#endif // OPENMS_MATH_MISC_STRUCT3D_H
