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

#ifndef OPENMS_DATASTRUCTURES_INTLIST_H
#define OPENMS_DATASTRUCTURES_INTLIST_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	/**
		@brief Int list
		
		This class is based on std::vector<Int> but adds some methods for convenience.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI IntList
		: public std::vector<Int>
	{
		public:

			///@name Constructors and assignment operators
			//@{
			/// Default constructor
			IntList();
			/// Copy constructor
			IntList(const IntList& rhs);
			/// Constructor from vector<UInt>
			IntList(const std::vector<UInt>& rhs);
			/// Constructor from vector<Int>
			IntList(const std::vector<Int>& rhs);
			///  Assignment operator
			IntList& operator=(const IntList& rhs);
			///  Assignment operator from vector<Int>
			IntList& operator=(const std::vector<Int>& rhs);
			///  Assignment operator from vector<UInt>
			IntList& operator=(const std::vector<UInt>& rhs);
			//@}
			
			///Operator for appending entries with less code
			template<typename IntType>
			IntList& operator<<(IntType value)
			{
				this->push_back(value);
				return *this;
			}

			/// Returns a list that is created by splitting the given comma-separated string (String are not trimmed!)
			static IntList create(const String& list);
			/// Returns if a string is contains in the list
			bool contains(Int s) const;
			
			
			/// output stream operator
			friend std::ostream& operator<<(std::ostream& os, const IntList& p);
			
	};


} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_INTLIST_H
