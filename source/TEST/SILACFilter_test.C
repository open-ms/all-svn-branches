// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SILACFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////






SILACFilter* ptr1 = 0;
SILACFilter* ptr2 = 0;

START_SECTION((SILACFilter(DoubleReal mass_separation_light_heavy, Int charge_)))
{
	ptr1 = new SILACFilter(1.0,1);
	TEST_NOT_EQUAL(ptr1, 0);
}
END_SECTION

START_SECTION((SILACFilter(DoubleReal mass_separation_light_medium, DoubleReal mass_separation_light_heavy, Int charge_)))
{
	ptr2 = new SILACFilter(1.0,2.0,1);
	TEST_NOT_EQUAL(ptr2, 0);
}
END_SECTION

START_SECTION((Int getSILACType()))
{
  TEST_EQUAL(ptr1->getSILACType(),2);
  TEST_EQUAL(ptr2->getSILACType(),3);
}
END_SECTION

START_SECTION(~SILACFilter())
{
	delete ptr;
}
END_SECTION

START_SECTION((void blockValue(DoubleReal value)))
{
  ptr1->blockValue(5.0);
  TEST_EQUAL(ptr1->blacklisted(5.0),true);
}
END_SECTION

START_SECTION((bool blacklisted(DoubleReal value)))
{
	TEST_EQUAL(ptr1->blacklisted(4.995),true);
	TEST_EQUAL(ptr1->blacklisted(5.0),true);
	TEST_EQUAL(ptr1->blacklisted(5.005),true);
}
END_SECTION

START_SECTION((void reset()))
{
  ptr1->reset();
  TEST_NOT_EQUAL(ptr1->blacklisted(4.995),true);
  	TEST_NOT_EQUAL(ptr1->blacklisted(5.0),true);
  	TEST_NOT_EQUAL(ptr1->blacklisted(5.005),true);
}
END_SECTION

START_SECTION((bool isFeature(DoubleReal act_rt, DoubleReal act_mz)))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<DoubleReal> getPeakValues()))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getEnvelopeDistanceLightMedium()))
{
	TEST_REAL_SIMILAR(ptr1->getEnvelopeDistanceLightMedium(),0.0);
	TEST_REAL_SIMILAR(ptr2->getEnvelopeDistanceLightMedium(),1.0);
}
END_SECTION

START_SECTION((DoubleReal getEnvelopeDistanceLightHeavy()))
{
	TEST_REAL_SIMILAR(ptr1->getEnvelopeDistanceLightMedium(),1.0);
	TEST_REAL_SIMILAR(ptr2->getEnvelopeDistanceLightMedium(),2.0);
}
END_SECTION

START_SECTION((DoubleReal getIsotopeDistance()))
{
	TEST_REAL_SIMILAR(ptr1->getIsotopeDistance(),1.0);
	TEST_REAL_SIMILAR(ptr2->getEnvelopeDistanceLightMedium(),2.0);
}
END_SECTION

START_SECTION((std::vector<DataPoint> getElements()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



