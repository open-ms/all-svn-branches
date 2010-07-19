// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/QTCluster.h>

namespace OpenMS {

QTCluster::QTCluster()
{

}
QTCluster::QTCluster(DataPoint* center_point_) : center_point(center_point_)
{
	cluster_members.insert(center_point);
}

QTCluster::~QTCluster() {
	// TODO Auto-generated destructor stub
}

DoubleReal QTCluster::getCenterRT()
{
	return center_point->rt;
}

DoubleReal QTCluster::getCenterMZ()
{
	return center_point->mz;
}

Size QTCluster::size() const
{
	return cluster_members.size();
}

bool QTCluster::operator<(const QTCluster &cluster) const
{
	return (this->size() < cluster.size());
}

void QTCluster::add(DataPoint* element)
{
	cluster_members.insert(element);
	traces[round(element->mz)].push_back(element);
	std::pair<DoubleReal,DoubleReal> diameters=getDiameters(element);
}

bool QTCluster::contains(DataPoint* element)
{
	std::set<DataPoint*>::iterator pos = cluster_members.find(element);
	return pos != cluster_members.end();
}

std::set<DataPoint*> QTCluster::getClusterMembers()
{
	return cluster_members;
}

std::pair<DoubleReal,DoubleReal> QTCluster::getDiameters(DataPoint* point)
{
	DoubleReal rt_diameter=std::numeric_limits<Real>::max();
	DoubleReal mz_diameter=0.0;
	DoubleReal point_mz=point->mz;
	DoubleReal point_rt=point->rt;
	for (std::set<DataPoint*>::iterator it=cluster_members.begin();it!=cluster_members.end();++it)
	{
		DoubleReal rt_dist=std::abs((*it)->rt-point_rt);
		if (rt_dist < rt_diameter)
			rt_diameter=rt_dist;
		DoubleReal mz_dist=std::abs((*it)->mz-point_mz);
		if (mz_dist > mz_diameter)
					mz_diameter=mz_dist;
	}
	return std::make_pair(rt_diameter,mz_diameter);
}

DoubleReal QTCluster::getMZstandardDeviation(DataPoint* point)
{
	DoubleReal deviation=std::numeric_limits<Real>::max();
	for (std::set<DataPoint*>::iterator it=cluster_members.begin();it!=cluster_members.end();++it)
	{
		DoubleReal act_mz=(*it)->mz;
		DoubleReal diff = std::abs(point->mz - act_mz);
		DoubleReal num = floor(diff);
		diff -= num;
		DoubleReal act_deviation=std::min(std::abs(diff-0.5),std::abs(diff));
		if (act_deviation < deviation)
			deviation=act_deviation;
	}
	return deviation;
}

void QTCluster::checkClusterShape()
{
	for (std::map<Int,std::list<DataPoint*> >::iterator trace_it=traces.begin();trace_it!=traces.end();++trace_it)
	{
		std::cout << "[" << trace_it->first << "]: " << trace_it->second.size() << std::endl;
	}
	std::cout << std::endl;
}

//DoubleReal QTCluster::getMZstandardDeviation()
//{
//
////	DoubleReal sd = 0.0;
////	DoubleReal mv=center_point->mz;
////	for (std::set<DataPoint*>::iterator it=cluster_members.begin();it!=cluster_members.end();++it)
////	{
////		sd += std::abs((*it)->mz - mv);
////	}
////	sd /= size();
////	return (sd);
//
//	DoubleReal mv=center_point->mz;
//	DoubleReal diff = std::abs((*it)->mz - mv);
//	DoubleReal num = floor(diff);
//	diff-=num;
//	return std::min(std::abs(diff-0.5),std::abs(diff));
//}

}
