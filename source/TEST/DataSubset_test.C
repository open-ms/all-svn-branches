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
#include <OpenMS/DATASTRUCTURES/DataSubset.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DataSubset, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DataSubset* ptr = 0;
START_SECTION(DataSubset())
{
	DataPoint p;
		p.mz=444.2;
		p.rt=23.13;
		p.feature_id=1;
	ptr = new DataSubset(p);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DataSubset())
{
	delete ptr;
}
END_SECTION


START_SECTION((DataSubset(DataPoint &data_point)))
{
  DataPoint p;
  		p.mz=444.2;
  		p.rt=23.13;
  		p.feature_id=1;
  DataSubset subset(p);
  TEST_REAL_SIMILAR(subset.rt,23.13);
  TEST_REAL_SIMILAR(subset.mz,444.2);
  TEST_EQUAL(subset.getID(),1);
}
END_SECTION

START_SECTION((DataSubset(const DataSubset &copy)))
{
	DataPoint p;
			p.mz=444.2;
			p.rt=23.13;
			p.feature_id=1;
	DataSubset subset(p);
	DataSubset copy(subset);
	TEST_REAL_SIMILAR(copy.rt,23.13);
	TEST_REAL_SIMILAR(copy.mz,444.2);
	TEST_EQUAL(copy.getID(),1);
}
END_SECTION

START_SECTION((DataSubset(const DataSubset *copy_ptr)))
{
	DataPoint p;
			p.mz=444.2;
			p.rt=23.13;
			p.feature_id=1;
	DataSubset subset(p);
	DataSubset copy(&subset);
	TEST_REAL_SIMILAR(copy.rt,23.13);
	TEST_REAL_SIMILAR(copy.mz,444.2);
	TEST_EQUAL(copy.getID(),1);
}
END_SECTION

START_SECTION((Int getID()))
{
	DataPoint p1;
	p1.mz=444.2;
	p1.rt=23.13;
	p1.feature_id=4;
	DataPoint p2;
	p2.mz=144.2;
	p2.rt=223.13;
	p2.feature_id=5;
	DataPoint p3;
	p3.mz=144.2;
	p3.rt=223.13;
	p3.feature_id=6;
	DataSubset subset(p1);
	subset.data_points.push_back(&p2);
	subset.data_points.push_back(&p3);
	TEST_EQUAL(subset.getID(),4);
}
END_SECTION

START_SECTION((Int operator<(const DataSubset &el) const ))
{
		DataPoint p1;
		p1.mz=444.2;
		p1.rt=23.13;
		p1.feature_id=4;
		DataPoint p2;
		p2.mz=144.2;
		p2.rt=223.13;
		p2.feature_id=5;
		DataPoint p3;
		p3.mz=144.2;
		p3.rt=223.13;
		p3.feature_id=6;
	DataSubset subset(p1);
	subset.data_points.push_back(&p2);
	subset.data_points.push_back(&p3);
	DataSubset subset1(p1);
	TEST_EQUAL(subset1<subset,true);
}
END_SECTION

START_SECTION((Int size()))
{
		DataPoint p1;
		p1.mz=444.2;
		p1.rt=23.13;
		p1.feature_id=4;
		DataPoint p2;
		p2.mz=144.2;
		p2.rt=223.13;
		p2.feature_id=5;
		DataPoint p3;
		p3.mz=144.2;
		p3.rt=223.13;
		p3.feature_id=6;
	DataSubset subset(p1);
	subset.data_points.push_back(&p2);
	subset.data_points.push_back(&p3);
	TEST_EQUAL(subset.size(),3);
}
END_SECTION

START_SECTION((bool operator!=(const DataSubset &el) const ))
{
	DataPoint p1;
	p1.mz=444.2;
	p1.rt=23.13;
	p1.feature_id=4;

	DataSubset subset(p1);
	DataPoint p2;
	p2.mz=144.2;
	p2.rt=223.13;
	p2.feature_id=5;

	subset.data_points.push_back(&p2);
	DataPoint p3;
	p3.mz=144.2;
	p3.rt=223.13;
	p3.feature_id=6;
	DataSubset subset1(p3);
	subset1.data_points.push_back(&p2);
	TEST_EQUAL(subset!=subset1,true);
}
END_SECTION

START_SECTION((bool operator==(const DataSubset &el) const ))
{
	DataPoint p1;
	p1.mz=444.2;
	p1.rt=23.13;
	p1.feature_id=4;
	DataPoint p2;
	p2.mz=144.2;
	p2.rt=223.13;
	p2.feature_id=5;
	DataSubset subset(p1);
	subset.data_points.push_back(&p2);
	DataSubset subset1(p1);
	subset1.data_points.push_back(&p2);
	TEST_EQUAL(subset==subset1,true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



