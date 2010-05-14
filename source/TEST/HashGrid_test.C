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
#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HashGrid, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HashGrid* ptr = 0;
START_SECTION(HashGrid())
{
	ptr = new HashGrid(1.,1.);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~HashGrid())
{
	delete ptr;
}
END_SECTION

HashGrid grid(80,50);

START_SECTION((HashGrid(DoubleReal rt_threshold_, DoubleReal mz_threshold_)))
{
	TEST_EQUAL(grid.getRTThreshold(),80)
	TEST_EQUAL(grid.getMZThreshold(),50)
}
END_SECTION

DataPoint* el1=new DataPoint(550.5,270.2,1.,1);
DataPoint* el2=new DataPoint(130.7,460.6,1.,2);
DataPoint* el3=new DataPoint(435.2,160.6,1.,3);
DataPoint* el4=new DataPoint(32.3,960.1,1.,4);
DataPoint* el5=new DataPoint(551.5,271.2,1.,5);


START_SECTION((void insert(GridElement *element_)))
{
	grid.insert(el1);
	grid.insert(el2);
	grid.insert(el3);
	grid.insert(el4);
	grid.insert(el5);
	GridCells::iterator pos=grid.find(make_pair(5,6));
	TEST_EQUAL(pos->second.size(),2)
	 pos=grid.find(make_pair(9,1));
	TEST_EQUAL(pos->second.size(),1)
	 pos=grid.find(make_pair(3,5));
	TEST_EQUAL(pos->second.size(),1)
	 pos=grid.find(make_pair(19,0));
	TEST_EQUAL(pos->second.size(),1)
}
END_SECTION


START_SECTION((void removeElement(GridElement *element_, Int x, Int y)))
{
	grid.removeElement(el1,5,6);
	//This should not work
	grid.removeElement(el2,0,0);
	GridCells::iterator pos=grid.find(make_pair(5,6));
	TEST_EQUAL(pos->second.size(),1)
	pos=grid.find(make_pair(9,1));
	TEST_EQUAL(pos->second.size(),1)
}
END_SECTION

START_SECTION((void removeElement(GridElement *element_)))
{
  grid.removeElement(el2);
  grid.removeElement(el3);
  GridCells::iterator pos=grid.find(make_pair(9,1));
  TEST_EQUAL(pos==grid.end(),true);
  pos=grid.find(make_pair(3,5));
  TEST_EQUAL(pos==grid.end(),true)
}
END_SECTION


START_SECTION((int size()))
{
	TEST_EQUAL(grid.size(),2)
}
END_SECTION

START_SECTION((DoubleReal getRTThreshold() const ))
{
	TEST_EQUAL(grid.getRTThreshold(),80)
}
END_SECTION

START_SECTION((DoubleReal getMZThreshold() const ))
{
	TEST_EQUAL(grid.getMZThreshold(),50)
}
END_SECTION

START_SECTION((Int getGridSizeX()))
{
	TEST_EQUAL(grid.getGridSizeX(),19)
}
END_SECTION

START_SECTION((Int getGridSizeY()))
{
	TEST_EQUAL(grid.getGridSizeY(),6)
}
END_SECTION

START_SECTION((Int getNumberOfElements()))
{
	TEST_EQUAL(grid.getNumberOfElements(),2)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



