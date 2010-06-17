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
  DataPoint p1;
  p1.mz=124.12;
  p1.rt=12.123;
  p1.feature_id=1;
  DataPoint copy(p1);
  TEST_REAL_SIMILAR(copy.rt,12.123);
  TEST_REAL_SIMILAR(copy.mz,124.12);
  TEST_EQUAL(copy.getID(),1);
}
END_SECTION

START_SECTION((DataPoint()))
{
	DataPoint p1;
	  p1.mz=124.12;
	  p1.rt=12.123;
	  p1.feature_id=1;
	  TEST_REAL_SIMILAR(p1.rt,12.123);
	  TEST_REAL_SIMILAR(p1.mz,124.12);
	  TEST_EQUAL(p1.getID(),1);
}
END_SECTION


START_SECTION((DataPoint& operator=(const DataPoint &rhs)))
{
	DataPoint p1;
	  p1.mz=124.12;
	  p1.rt=12.123;
	  p1.feature_id=1;
	  DataPoint copy=p1;
	  TEST_REAL_SIMILAR(copy.rt,12.123);
	  TEST_REAL_SIMILAR(copy.mz,124.12);
	  TEST_EQUAL(copy.getID(),1);
}
END_SECTION

START_SECTION((bool operator==(const DataPoint &rhs) const ))
{
	DataPoint p1;
	  p1.mz=124.12;
	  p1.rt=12.123;
	  p1.feature_id=1;
	  DataPoint p2;
	    p2.mz=124.12;
	    p2.rt=12.123;
	    p2.feature_id=1;
	TEST_EQUAL(p1==p2,true);
}
END_SECTION

START_SECTION((bool operator!=(const DataPoint &rhs) const ))
{
	DataPoint p1;
	p1.mz=124.12;
	p1.rt=12.123;
	p1.feature_id=1;
	DataPoint p2;
	p2.mz=124.12;
	p2.rt=112.123;
	p2.feature_id=1;
	TEST_EQUAL(p1!=p2,true);
	DataPoint p3;
	p3.mz=1124.12;
	p3.rt=12.123;
	p3.feature_id=1;
	TEST_EQUAL(p1!=p3,true);
	DataPoint p5;
	p5.mz=124.12;
	p5.rt=12.123;
	p5.feature_id=2;
	TEST_EQUAL(p1!=p5,true);
}
END_SECTION

START_SECTION((bool operator<(const DataPoint &rhs) const ))
{
	DataPoint p1;
	p1.mz=124.12;
	p1.rt=12.123;
	p1.feature_id=1;
	DataPoint p2;
	p2.mz=124.12;
	p2.rt=112.123;
	p2.feature_id=2;
	TEST_EQUAL(p1<p2,true);
}
END_SECTION

START_SECTION((Int getID()))
{
	DataPoint p1;
		p1.mz=124.12;
		p1.rt=12.123;
		p1.feature_id=2;
	TEST_EQUAL(p1.getID(),2);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



