// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DOUBLELIST_H
#define OPENMS_DATASTRUCTURES_DOUBLELIST_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	/**
		@brief DoubleReal list
		
		This class is based on std::vector<DoubleReal> but adds some methods for convenience.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI DoubleList
		: public std::vector<DoubleReal>
	{
		public:

			///@name Constructors and assignment operators
			//@{
			/// Default constructor
			DoubleList();
			/// Copy constructor
			DoubleList(const DoubleList& rhs);
			/// Constructor from vector<Real>
			DoubleList(const std::vector<Real>& rhs);
			/// Constructor from vector<DoubleReal>
			DoubleList(const std::vector<DoubleReal>& rhs);
			///  Assignment operator
			DoubleList& operator=(const DoubleList& rhs);
			///  Assignment operator from vector<DoubleReal>
			DoubleList& operator=(const std::vector<DoubleReal>& rhs);
			///  Assignment operator from vector<Real>
			DoubleList& operator=(const std::vector<Real>& rhs);
			//@}
			
			///Operator for appending entries with less code
			template<typename DoubleType>
			DoubleList& operator<<(DoubleType value)
			{
				this->push_back(value);
				return *this;
			}

			/// Returns a list that is created by splitting the given comma-separated string (String are not trimmed!)
			static DoubleList create(const String& list);
			/// Returns if a string is contains in the list
			bool contains(DoubleReal s) const;
			
			
			/// output stream operator
			friend std::ostream& operator<<(std::ostream& os, const DoubleList& p);
			
	};


} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DOUBLELIST_H
