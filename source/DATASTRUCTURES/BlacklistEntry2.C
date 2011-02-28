// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/BlacklistEntry2.h>

namespace OpenMS
{
  	
  BlacklistEntry2::BlacklistEntry2()
  {
    blackArea = DRange<2>(0,0,0,0);
 	  charge = 1;
	  mass_separations = std::vector<DoubleReal>(0);
    relative_peak_position = 0;
  }
  	
  BlacklistEntry2::BlacklistEntry2(const BlacklistEntry2 &copyin) : GridElement(copyin)
  {
    blackArea = copyin.blackArea;
	  charge = copyin.charge;
    mass_separations = copyin.mass_separations;
    relative_peak_position = copyin.relative_peak_position;
  }

  BlacklistEntry2& BlacklistEntry2::operator=(const BlacklistEntry2 &rhs)
  {
    this->blackArea = rhs.blackArea;
    this->charge = rhs.charge;
    this->mass_separations = rhs.mass_separations;
    this->relative_peak_position = rhs.relative_peak_position;
	  return *this;
  }

  bool BlacklistEntry2::operator==(const BlacklistEntry2 &rhs) const
  {
	  if( this->blackArea != rhs.blackArea) return false;
	  if( this->charge != rhs.charge) return false;
	  if( this->mass_separations != rhs.mass_separations) return false;
	  if( this->relative_peak_position != rhs.relative_peak_position) return false;
	  return true;
  }

  bool BlacklistEntry2::operator!=(const BlacklistEntry2 &rhs) const
  {
	  return !(*this==rhs);
  }
  
  bool BlacklistEntry2::operator<(const BlacklistEntry2 &rhs) const // required for built-in STL functions like sort
  {
	  if ( *this == rhs) return false;

	  return this->charge < rhs.charge;
  }

	int BlacklistEntry2::getID() const
  {
	  return id;
  }

}


