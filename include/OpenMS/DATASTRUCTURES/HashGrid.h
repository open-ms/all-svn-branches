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


#include<map>
#include<list>
#include<set>
#include<vector>
#include<iostream>

#include <OpenMS/DATASTRUCTURES/DataSubset.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>



#ifndef HASHGRID_H_
#define HASHGRID_H_

namespace OpenMS
{

class HashGrid {

private:

	DoubleReal rt_threshold;
	DoubleReal mz_threshold;
	DoubleReal hyp_threshold;
	DoubleReal rt_scaling;
	DoubleReal min_distance;
	std::pair<DataSubset*,DataSubset*> min_distance_subsets;

public:
	DistanceSet distances;
	int grid_size_x;
	int grid_size_y;
	std::map<std::pair<Int,Int>, std::list<GridElement*> > elements;
	HashGrid();
	HashGrid(DoubleReal rt_threshold_,DoubleReal mz_threshold_);
	~HashGrid();
	void setRTThreshold(DoubleReal threshold_);
	void setMZThreshold(DoubleReal threshold_);
	DoubleReal getDistance(DataSubset& element1,DataSubset& element2);
	void removeElement(GridElement* element_,Int x,Int y);
	void removeElement(GridElement* element_);
	void insert(GridElement* element_);
	void consoleOut();
	int size();
	int distanceSize();
	DoubleReal getRT_threshold() const;
	DoubleReal getMZ_threshold() const;

};
}




#endif /* HASHGRID_H_ */
