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


#include <OpenMS/COMPARISON/CLUSTERING/QTClustering.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS {
QTClustering::QTClustering(std::vector<DataPoint>& data,DoubleReal rt_diameter_, DoubleReal mz_diameter_,DoubleReal isotope_distance_) : grid(HashGrid(rt_diameter_,mz_diameter_))
{
//	if(data.size()<2)
//	{
//		throw InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "The data set contains not enough elements");
//	}
	rt_diameter=rt_diameter_;
	mz_diameter=mz_diameter_;
	isotope_distance=isotope_distance_;
	for (std::vector<DataPoint>::iterator it=data.begin();it!=data.end();++it)
	{
		grid.insert(&(*it));
	}
}

bool QTClustering::isBelowDiameter(DataPoint* point1, DataPoint* point2)
{
	if (std::abs(point1->mz - point2->mz) > mz_diameter)
		return false;
	if (std::abs(point1->rt - point2->rt) > rt_diameter)
		return false;
	return true;
}

std::vector<std::vector<DataPoint*> > QTClustering::performClustering()
{
	Int number_of_elements=grid.getNumberOfElements();
	startProgress(0,number_of_elements,"clustering data");
	std::vector<std::vector<DataPoint*> > cluster_vector;
	//HashGrid act_grid=grid;
	while(grid.getNumberOfElements() > 0)
	{
		setProgress(number_of_elements - grid.getNumberOfElements());
		QTCluster act_cluster=QTClust(grid);
		clusters.push_back(act_cluster);
		const std::set<DataPoint*>& cluster_members=act_cluster.getClusterMembers();
		for (std::set<DataPoint*>::const_iterator it=cluster_members.begin();it!=cluster_members.end();++it)
		{
			grid.removeElement(*it);
		}
	}
	Size cluster_id=0;
	for (std::list<QTCluster>::iterator list_it=clusters.begin();list_it!=clusters.end();++list_it)
	{
//		std::cout << cluster_id << ": "<< std::endl;
//		list_it->checkClusterShape();
		std::vector<DataPoint*> act_vector;
		const std::set<DataPoint*>& cluster_members=list_it->getClusterMembers();
		for (std::set<DataPoint*>::const_iterator set_it=cluster_members.begin();set_it!=cluster_members.end();++set_it)
		{
			(*set_it)->cluster_id=cluster_id;
			act_vector.push_back(*set_it);
		}
		++cluster_id;
		cluster_vector.push_back(act_vector);
	}
	endProgress();
	return cluster_vector;
}

QTCluster QTClustering::QTClust(HashGrid& act_grid)
{
	if (act_grid.getNumberOfElements() <= 1)
	{
			DataPoint* element_ptr = dynamic_cast<DataPoint*> (*(act_grid.begin()->second.begin()));
			return QTCluster(element_ptr);
	}

	std::set<QTCluster> cluster_set;
		for (GridCells::iterator it=act_grid.begin();it!=act_grid.end();++it)
		{
			std::pair<int,int> act_coords=it->first;
			std::list<GridElement*>& elements=it->second;
			int x=act_coords.first;
			int y=act_coords.second;
			for (std::list<GridElement*>::iterator current_element=elements.begin();current_element!=elements.end();++current_element)
			{
				DataPoint* element_ptr = dynamic_cast<DataPoint*> (*current_element);
				QTCluster act_cluster(element_ptr);
				for (int i=x-1;i<=x+1;++i)
				{
					if (i<0 || i>act_grid.getGridSizeX())
						continue;
					for (int j=y-1;j<=y+1;++j)
					{
						if (j<0 || j>act_grid.getGridSizeY())
							continue;
						GridCells::iterator act_pos=grid.find(std::make_pair(i,j));
						if (act_pos==act_grid.end())
							continue;
						std::list<GridElement*>& neighbor_elements=act_pos->second;

						for (std::list<GridElement*>::iterator neighbor_element=neighbor_elements.begin();neighbor_element!=neighbor_elements.end();++neighbor_element)
						{
							if (*current_element==*neighbor_element)
								continue;
							DataPoint* neighbor_ptr = dynamic_cast<DataPoint*> (*neighbor_element);
	//						if (isBelowDiameter(element_ptr,neighbor_ptr))
							std::pair<DoubleReal,DoubleReal> diameters=act_cluster.getDiameters(neighbor_ptr);
	//						std::cout << diameters.first << " " << diameters.second << "\t";
							if (diameters.first <= 20 && diameters.second <= mz_diameter && act_cluster.getMZstandardDeviation(neighbor_ptr,isotope_distance)<0.02)
								act_cluster.add(neighbor_ptr);
						}
					}
				}
	//			std::cout << std::endl;
				cluster_set.insert(act_cluster);
			}
		}
		return *cluster_set.rbegin();
}

QTClustering::~QTClustering()
{
	// TODO Auto-generated destructor stub
}
}
