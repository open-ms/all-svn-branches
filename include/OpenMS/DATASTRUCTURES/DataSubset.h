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


#ifndef DATASUBSET_H_
#define DATASUBSET_H_

#include <OpenMS/DATASTRUCTURES/DataPoint.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>

namespace OpenMS
{
class DataSubset;

struct Dist{};

struct DistanceEntry
{
	DataSubset* data_point;
	DataSubset* owner;
	DoubleReal distance;
	DistanceEntry(DataSubset* owner_,DataSubset* data_point_,DoubleReal distance_)
	{
		owner=owner_;
		data_point=data_point_;
		distance=distance_;
	}
	Int operator<(const DistanceEntry &i) const
	{
		return distance < i.distance;
	}

};


typedef boost::multi_index::multi_index_container<
DistanceEntry,
boost::multi_index::indexed_by<
boost::multi_index::hashed_unique<
boost::multi_index::composite_key<
DistanceEntry,
boost::multi_index::member<DistanceEntry,DataSubset*,&DistanceEntry::owner>,
boost::multi_index::member<DistanceEntry,DataSubset*,&DistanceEntry::data_point> > >,
boost::multi_index::ordered_non_unique< boost::multi_index::tag<Dist>,
boost::multi_index::member<DistanceEntry,DoubleReal,&DistanceEntry::distance> > >
> DistanceSet;


class DataSubset : public GridElement
{

public:
	std::map<GridElement*,DistanceSet::iterator> distance_iterators;
	std::list<DataPoint*> data_points;
	std::vector<BinaryTreeNode> tree;

	DataSubset(DataPoint& data_point);
	DataSubset(const DataSubset& copy);
	DataSubset(const DataSubset* copy_ptr);
	Int operator<(const DataSubset &el) const;
	Int size();
	Int getID();
	bool operator !=(const DataSubset &el) const;
	bool operator ==(const DataSubset & el) const;
};
}

#endif /* DATASUBSET_H_ */
