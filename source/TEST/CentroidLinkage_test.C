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
#include <OpenMS/COMPARISON/CLUSTERING/CentroidLinkage.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CentroidLinkage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CentroidLinkage* ptr = 0;
START_SECTION(CentroidLinkage())
{
	ptr = new CentroidLinkage(0.1);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~CentroidLinkage())
{
	delete ptr;
}
END_SECTION



START_SECTION((CentroidLinkage(DoubleReal rt_scaling_)))
{
	CentroidLinkage cl_temp(0.05);
  TEST_EQUAL(cl_temp.rt_scaling,0.05)
}
END_SECTION

CentroidLinkage cl(0.05);
vector<DataPoint> points;
points.push_back(DataPoint(10,30,5,0));
points.push_back(DataPoint(40,20,5,1));

START_SECTION((DoubleReal getDistance(DataSubset &subset1, DataSubset &subset2)))
{
  DataSubset subset1(points[0]);
  DataSubset subset2(points[1]);
  TEST_REAL_SIMILAR(cl.getDistance(subset1,subset2),10.111874);
}
END_SECTION

START_SECTION((DoubleReal getDistance(DataPoint &point1, DataPoint &point2)))
{
	TEST_REAL_SIMILAR(cl.getDistance(points[0],points[1]),10.111874);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



