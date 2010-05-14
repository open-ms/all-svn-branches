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
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DataPoint, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DataPoint* ptr = 0;
START_SECTION(DataPoint())
{
	ptr = new DataPoint();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DataPoint())
{
	delete ptr;
}
END_SECTION

START_SECTION((DataPoint(const DataPoint &copyin)))
{
  DataPoint p(12.123,124.12,11.22,1);
  DataPoint copy(p);
  TEST_REAL_SIMILAR(copy.rt,12.123);
  TEST_REAL_SIMILAR(copy.mz,124.12);
  TEST_REAL_SIMILAR(copy.intensity,11.22);
  TEST_EQUAL(copy.getID(),1);
}
END_SECTION

START_SECTION((DataPoint(DoubleReal rt_, DoubleReal mz_, DoubleReal intensity_, Int feature_id_, DoubleReal real_rt_)))
{
	DataPoint p(12.123,124.12,11.22,1);
	  TEST_REAL_SIMILAR(p.rt,12.123);
	  TEST_REAL_SIMILAR(p.mz,124.12);
	  TEST_REAL_SIMILAR(p.intensity,11.22);
	  TEST_EQUAL(p.getID(),1);
}
END_SECTION


START_SECTION((DataPoint& operator=(const DataPoint &rhs)))
{
	DataPoint p(12.123,124.12,11.22,1);
	  DataPoint copy=p;
	  TEST_REAL_SIMILAR(copy.rt,12.123);
	  TEST_REAL_SIMILAR(copy.mz,124.12);
	  TEST_REAL_SIMILAR(copy.intensity,11.22);
	  TEST_EQUAL(copy.getID(),1);
}
END_SECTION

START_SECTION((bool operator==(const DataPoint &rhs) const ))
{
	DataPoint p1(12.123,124.12,11.22,1);
	DataPoint p2(12.123,124.12,11.22,1);
	TEST_EQUAL(p1==p2,true);
}
END_SECTION

START_SECTION((bool operator!=(const DataPoint &rhs) const ))
{
	DataPoint p1(12.123,124.12,11.22,1);
		DataPoint p2(112.123,124.12,11.22,1);
		TEST_EQUAL(p1!=p2,true);
		DataPoint p3(12.123,1214.12,11.22,1);
		TEST_EQUAL(p1!=p3,true);
		DataPoint p4(12.123,124.12,131.22,1);
		TEST_EQUAL(p1!=p4,true);
		DataPoint p5(12.123,124.12,11.22,2);
		TEST_EQUAL(p1!=p5,true);
}
END_SECTION

START_SECTION((bool operator<(const DataPoint &rhs) const ))
{
	DataPoint p1(12.123,124.12,131.22,1);
	DataPoint p2(12.123,124.12,11.22,2);
	TEST_EQUAL(p1<p2,true);
	DataPoint p3(12.123,124.12,11.22,1);
	TEST_NOT_EQUAL(p1<p3,true);
}
END_SECTION

START_SECTION((Int getID()))
{
	DataPoint p(12.123,124.12,11.22,2);
	TEST_EQUAL(p.getID(),2);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



