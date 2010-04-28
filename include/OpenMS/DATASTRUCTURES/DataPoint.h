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


#ifndef DATAPOINT_H_
#define DATAPOINT_H_

#include <OpenMS/DATASTRUCTURES/GridElement.h>

namespace OpenMS
{

class DataPoint : public GridElement {

public:
	DoubleReal intensity; // intensity at RT and m/z
	Int cluster_id; // ID number of the cluster the data point belongs to
	Int cluster_size; // number of points in cluster 'cluster_id'
	Int feature_id;
	DataPoint();
	DataPoint(const DataPoint &copyin);
	DataPoint(DoubleReal rt_, DoubleReal mz_, DoubleReal intensity_, Int feature_id_);
	~DataPoint(){};
	DataPoint& operator=(const DataPoint &rhs);
	bool operator==(const DataPoint &rhs) const;
	bool operator!=(const DataPoint &rhs) const;
	bool operator<(const DataPoint &rhs) const;
	Int getID();

};
}


#endif /* DATAPOINT_H_ */
