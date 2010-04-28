/*
 * GeometricHashing.cpp
 *
 *  Created on: 24.03.2010
 *      Author: steffen
 */

#include <OpenMS/DATASTRUCTURES/HashGrid.h>

namespace OpenMS
{
HashGrid::HashGrid()
{
}

HashGrid::HashGrid(DoubleReal rt_threshold_,DoubleReal mz_threshold_)
{
	rt_threshold=rt_threshold_;
	mz_threshold=mz_threshold_;
	grid_size_x=-1;
	grid_size_y=-1;
}
HashGrid::~HashGrid() {
	for (std::map<std::pair<Int,Int>, std::list<GridElement*> >::iterator it=elements.begin();it!=elements.end();++it)
	{
		std::list<GridElement*>& elements=it->second;
		for (std::list<GridElement*>::iterator lit=elements.begin();lit!=elements.end();++lit)
		{
			delete(*lit);
		}
	}

}
void HashGrid::setRTThreshold(DoubleReal threshold_)
{
	rt_threshold=threshold_;
}

void HashGrid::setMZThreshold(DoubleReal threshold_)
{
	mz_threshold=threshold_;
}

DoubleReal HashGrid::getDistance(DataSubset& element1,DataSubset& element2)
{
	DistanceSet::iterator pos=distances.find(boost::make_tuple(&element1,&element2));
	if (pos!=distances.end())
	{
		DistanceEntry el=*pos;
		return el.distance;
	}
	pos=distances.find(boost::make_tuple(&element2,&element1));
	if (pos!=distances.end())
	{
		DistanceEntry el=*pos;
		return el.distance;
	}
	return -1;
}

void HashGrid::removeElement(GridElement* element_,Int x,Int y)
{
	std::list<GridElement*>& subsets = elements[std::make_pair(x,y)];
	subsets.remove(element_);
	if (subsets.empty())
		elements.erase(std::make_pair(x,y));
}

void HashGrid::removeElement(GridElement* element_)
{
	int x = element_->mz / mz_threshold;
	int y = element_->rt / rt_threshold;
	removeElement(element_,x,y);
}



void HashGrid::insert(GridElement* element_)
{
	int x = element_->mz / mz_threshold;
	if (x>grid_size_x)
		grid_size_x=x;
	int y = element_->rt / rt_threshold;
	if (y>grid_size_y)
		grid_size_y=y;

	elements[std::make_pair(x,y)].push_back(element_);
}

void HashGrid::consoleOut()
{
	for (std::map<std::pair<Int,Int>, std::list<GridElement*> >::iterator it=elements.begin();it!=elements.end();++it)
	{
		std::pair<Int,Int> coords=it->first;
		std::list<GridElement*> act_elements= it->second;
		if (it->second.size()>0)
			std::cout << coords.first << "/" << coords.second<< ": ";
		for (std::list<GridElement*>::iterator list_it=act_elements.begin();list_it!=act_elements.end();++list_it)
		{
			std::cout << (*list_it)->getID() << " | ";
		}
		std::cout << std::endl;

	}
	std::cout << std::endl;
}

int HashGrid::size()
{
	return elements.size();
}

int HashGrid::distanceSize()
{
	return distances.size();
}


DoubleReal HashGrid::getRT_threshold() const
{
	return rt_threshold;
}

DoubleReal HashGrid::getMZ_threshold() const
{
	return mz_threshold;
}
}
