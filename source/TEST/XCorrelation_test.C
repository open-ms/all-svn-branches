X-Powered-By: PHP/5.1.2
Content-type: text/html

// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TODO/XCorrelation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(XCorrelation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

XCorrelation* ptr = 0;
START_SECTION(XCorrelation())
{
	ptr = new XCorrelation();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~XCorrelation())
{
	delete ptr;
}
END_SECTION

START_SECTION((XCorrelation(const XCorrelation &source)))
{
  // TODO
}
END_SECTION

START_SECTION((~XCorrelation()))
{
  // TODO
}
END_SECTION

START_SECTION((XCorrelation& operator=(const XCorrelation &source)))
{
  // TODO
}
END_SECTION

START_SECTION((void getXCorrelation(const SpectrumType &s1, const SpectrumType &s2, DoubleReal &best_score1, DoubleReal &best_score2, DoubleReal &best_shift, std::list< std::pair< Size, Size > > &best_matches) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void maxSparseMatches(const SpectrumType &s1, const SpectrumType &s2, std::vector< std::pair< Size, Size > > &all_matches, std::list< std::pair< Size, Size > > &best_matches) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool matchPredicate(std::pair< Size, Size > first, std::pair< Size, Size > second)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST