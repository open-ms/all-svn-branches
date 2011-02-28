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


#ifndef OPENMS_DATASTRUCTURES_BLACKLISTENTRY2_H
#define OPENMS_DATASTRUCTURES_BLACKLISTENTRY2_H

#include <OpenMS/DATASTRUCTURES/GridElement.h>

namespace OpenMS
{

/**
 * @brief blacklisted range in the rt-mz plane which can be stored in a HashGrid
 * @see HashGrid
 * @ingroup Datastructures
 */

class OPENMS_DLLAPI BlacklistEntry2 : public GridElement {

public:
	/**
	 * @brief blacklisted area in rt-mz plane
	 */
	DRange<2> blackArea;
  
	/**
	 * @brief charge of the filter causing the blacklisting
	 */
	Int charge;
  
	/**
	 * @brief mass separations [Da] of the filter causing the blacklisting
	 */	
	std::vector<DoubleReal> mass_separations;
  
	/**
	 * @brief position of the blacklisted area relative to the mono-isotopic peak of the light peptide
	 */	
	DoubleReal relative_peak_position;

	/**
	 * @brief unique ID of blacklisted area
	 */
	Int id;
  
  /**
	 * @brief default constructor
	 */
	BlacklistEntry2();
	
	/**
	 * @brief copy constructor
	 * @param this BlacklistEntry2 will be copied
	 */
	BlacklistEntry2(const BlacklistEntry2 &copyin);
	
	// destructor
	~BlacklistEntry2(){};
  
  BlacklistEntry2& operator = (const BlacklistEntry2 &rhs);
  bool operator == (const BlacklistEntry2 &rhs) const;
  bool operator != (const BlacklistEntry2 &rhs) const;
  bool operator < (const BlacklistEntry2 &rhs) const;

	/**
	 * @brief gets the unique ID of the blacklisted area
	 */
	Int getID() const;
  
};
}


#endif /* DATAPOINT_H_ */
