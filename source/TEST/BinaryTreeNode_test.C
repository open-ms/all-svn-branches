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
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
///////////////////////////


using namespace OpenMS;
using namespace std;

START_TEST(BinaryTreeNode, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinaryTreeNode* ptr = 0;
START_SECTION(BinaryTreeNode())
{
	ptr = new BinaryTreeNode();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~BinaryTreeNode())
{
	delete ptr;
}
END_SECTION

DataPoint* p1 = new DataPoint;
p1->mz=20;
p1->rt=30;
p1->feature_id=1;
DataPoint* p2 = new DataPoint;
p2->mz=15.5;
p2->rt=2.3;
p2->feature_id=2;
DataPoint* p3 = new DataPoint;
p3->mz=6.6;
p3->rt=8.8;
p3->feature_id=3;
BinaryTreeNode node(p1,p2,123.45);

START_SECTION((BinaryTreeNode(DataPoint *data1_, DataPoint *data2_, DoubleReal distance_)))
{
	TEST_EQUAL(node.data1->getID(),1);
	  TEST_EQUAL(node.data2->getID(),2);
	  TEST_REAL_SIMILAR(node.distance,123.45);
}
END_SECTION

START_SECTION((BinaryTreeNode(const BinaryTreeNode &source)))
{
  BinaryTreeNode cp(node);
  TEST_EQUAL(cp.data1->getID(),1);
  TEST_EQUAL(cp.data2->getID(),2);
  TEST_REAL_SIMILAR(cp.distance,123.45);
}
END_SECTION

START_SECTION((BinaryTreeNode& operator=(const BinaryTreeNode &source)))
{
	BinaryTreeNode cp;
	  cp=node;
	  TEST_EQUAL(node.data1->getID(),1);
	  TEST_EQUAL(node.data2->getID(),2);
	  TEST_REAL_SIMILAR(node.distance,123.45);
}
END_SECTION


START_SECTION((bool operator==(const BinaryTreeNode &cp) const ))
{
	BinaryTreeNode cp(p1,p2,123.45);
	TEST_EQUAL(node==cp,true);
}
END_SECTION

START_SECTION((bool operator!=(const BinaryTreeNode &cp) const ))
{
	BinaryTreeNode cp1(p1,p3,123.45);
	TEST_EQUAL(node!=cp1,true)
	BinaryTreeNode cp2(p2,p2,123.45);
	TEST_EQUAL(node!=cp2,true);
	BinaryTreeNode cp3(p1,p2,111.22);
	TEST_EQUAL(node!=cp3,true);
}
END_SECTION

START_SECTION((bool operator<(const BinaryTreeNode &cp) const ))
{
	BinaryTreeNode cp1(p1,p3,153.45);
	TEST_EQUAL(node<cp1,true);
	BinaryTreeNode cp2(p2,p3,123.45);
	TEST_EQUAL(node<cp2,true);
	BinaryTreeNode cp3(p1,p2,113.45);
	TEST_NOT_EQUAL(node<cp3,true);
	BinaryTreeNode cp4(p1,p2,123.45);
	TEST_NOT_EQUAL(node<cp4,true);;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



