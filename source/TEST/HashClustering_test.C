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
#include <OpenMS/COMPARISON/CLUSTERING/HashClustering.h>
#include <OpenMS/COMPARISON/CLUSTERING/CentroidLinkage.h>
#include <OpenMS/DATASTRUCTURES/DataSubset.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HashClustering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CentroidLinkage cl(1);

HashClustering* ptr = 0;
START_SECTION((HashClustering(std::vector< DataPoint > &data, int rt_threshold, int mz_threshold, ClusteringMethod &method_)))
{
	vector<DataPoint> v;
	ptr = new HashClustering(v,1,1,cl);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~HashClustering())
{
	delete ptr;
}
END_SECTION
START_SECTION((std::vector< Real > averageSilhouetteWidth(std::vector< BinaryTreeNode > &tree)))
{
	DataPoint* el0=new DataPoint(550.5,270.2,1.,0,1.);
  DataPoint* el1=new DataPoint(550.5,270.2,1.,1,1.);
  DataPoint* el2=new DataPoint(130.7,460.6,1.,2,1.);
  DataPoint* el3=new DataPoint(435.2,160.6,1.,3,1.);
  DataPoint* el4=new DataPoint(32.3,960.1,1.,4,1.);
  DataPoint* el5=new DataPoint(551.5,271.2,1.,5,1.);
	vector<BinaryTreeNode> tree;
	tree.push_back(BinaryTreeNode(el0,el1,0.707107));
	tree.push_back(BinaryTreeNode(el0,el5,1.41421));
tree.push_back(BinaryTreeNode(el0,el3,159.551));
tree.push_back(BinaryTreeNode(el0,el2,447.644));
tree.push_back(BinaryTreeNode(el0,el4,789.234));
vector<Real> result;
result.push_back(0.333333);
result.push_back(0.49705);
result.push_back(0.546109);
result.push_back(0.552759);
result.push_back(0);
vector<Real> asw = ptr->averageSilhouetteWidth(tree);
TEST_EQUAL(result.size(), asw.size());
for (int i=0;i<asw.size();++i)
{
	TEST_REAL_SIMILAR(result[i], asw[i]);
}

}
END_SECTION

START_SECTION((void cut(int cluster_quantity, std::vector< std::vector< DataPoint * > > &clusters, std::vector< BinaryTreeNode > &tree)))
{
	DataPoint* el0=new DataPoint(550.5,270.2,1.,0,1.);
  DataPoint* el1=new DataPoint(550.5,270.2,1.,1,1.);
  DataPoint* el2=new DataPoint(130.7,460.6,1.,2,1.);
  DataPoint* el3=new DataPoint(435.2,160.6,1.,3,1.);
  DataPoint* el4=new DataPoint(32.3,960.1,1.,4,1.);
  DataPoint* el5=new DataPoint(551.5,271.2,1.,5,1.);

	vector< vector<DataPoint*> > clusters;
	vector< vector<DataPoint*> > result;
	vector<DataPoint*> v1;
	v1.push_back(el0);
	v1.push_back(el1);
	v1.push_back(el2);
	vector<DataPoint*> v2;
	v2.push_back(el3);
	v2.push_back(el4);
	vector<DataPoint*> v3;
	v3.push_back(el5);

	result.push_back(v2);
	result.push_back(v1);
	result.push_back(v3);

  vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(el1,el2,0.3f));
	tree.push_back(BinaryTreeNode(el3,el4,0.4f));
	tree.push_back(BinaryTreeNode(el0,el1,0.5f));
	tree.push_back(BinaryTreeNode(el0,el3,0.6f));
	tree.push_back(BinaryTreeNode(el0,el5,0.7f));

  ptr->cut(3, clusters, tree);
	TEST_EQUAL(clusters.size(), result.size());
	for (Size i = 0; i < clusters.size(); ++i)
	{
		TEST_EQUAL(clusters[i].size(), result[i].size());
		for (Size j = 0; j < clusters[i].size(); ++j)
		{
			TEST_EQUAL(result[i][j]->getID(),clusters[i][j]->getID());
		}
	}

}
END_SECTION


START_SECTION((std::vector<std::vector<BinaryTreeNode> > performClustering()))
{
  DataPoint* el0=new DataPoint(1.5,1.2,1.,0,1.);
  DataPoint* el1=new DataPoint(2.5,2.2,1.,1,1.);
  DataPoint* el2=new DataPoint(10.7,10.6,1.,2,1.);
  DataPoint* el3=new DataPoint(10.2,10.6,1.,3,1.);
  DataPoint* el4=new DataPoint(20.3,20.1,1.,4,1.);
  DataPoint* el5=new DataPoint(20.5,20.2,1.,5,1.);
  vector<DataPoint> v;
  	v.push_back(*el0);
  	v.push_back(*el1);
  	v.push_back(*el2);
  	v.push_back(*el3);
  	v.push_back(*el4);
  	v.push_back(*el5);
  	vector< BinaryTreeNode > result;
  		result.push_back(BinaryTreeNode(el4,el5,0.223607));
  		result.push_back(BinaryTreeNode(el2,el3,0.5));
  		result.push_back(BinaryTreeNode(el0,el1,1.41421));
  		result.push_back(BinaryTreeNode(el0,el2,12.2724));
  		result.push_back(BinaryTreeNode(el0,el4,19.9231));
  	HashClustering hc(v,2,2,cl);
  	vector<vector<BinaryTreeNode> > trees;
  	hc.performClustering(trees);
  	TEST_EQUAL(trees[0].size(), result.size());
		for (Size j = 0; j < trees[0].size(); ++j)
		{
			TEST_EQUAL(trees[0][j].data1->getID(), result[j].data1->getID());
			TEST_EQUAL(trees[0][j].data2->getID(), result[j].data2->getID());
			TEST_REAL_SIMILAR(trees[0][j].distance, result[j].distance);
		}
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



