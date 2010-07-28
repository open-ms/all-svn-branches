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


#ifndef OPENMS_DATASTRUCTURES_QTSUBSET_H
#define OPENMS_DATASTRUCTURES_QTSUBSET_H

#include<OpenMS/DATASTRUCTURES/DataPoint.h>

namespace OpenMS {
class OPENMS_DLLAPI QTCluster {
private:
	DataPoint* center_point;
	std::set<DataPoint*> cluster_members;
	std::map<Int,std::list<DataPoint*> > traces;
public:
	QTCluster();
	QTCluster(DataPoint* center_point_);
	virtual ~QTCluster();
	DoubleReal getCenterRT();
	DoubleReal getCenterMZ();
	Size size() const;
	bool operator<(const QTCluster &cluster) const;
	void add(DataPoint* element);
	std::set<DataPoint*> getClusterMembers();
	std::pair<DoubleReal,DoubleReal> getDiameters(DataPoint* point);
	DoubleReal getMZstandardDeviation(DataPoint* point,DoubleReal isotope_distance);
	bool contains(DataPoint* element);
	bool checkClusterShape(DataPoint* point);
};
}
#endif /* QTSUBSET_H_ */
