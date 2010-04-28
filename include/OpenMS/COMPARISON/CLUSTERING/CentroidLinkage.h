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

#ifndef CENTROIDLINKAGE_H_
#define CENTROIDLINKAGE_H_

#include <OpenMS/COMPARISON/CLUSTERING/ClusteringMethod.h>

namespace OpenMS
{

class CentroidLinkage : public ClusteringMethod {
public:
	CentroidLinkage(DoubleReal rt_scaling_);
	DoubleReal getDistance(DataSubset& subset1,DataSubset& subset2);
	DoubleReal getDistance(DataPoint& point1,DataPoint& point2);
};
}

#endif /* CENTROIDLINKAGE_H_ */
