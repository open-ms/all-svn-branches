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


#include <OpenMS/DATASTRUCTURES/DataPoint.h>

namespace OpenMS
{
DataPoint::DataPoint()
{
	rt = 0;
	mz = 0;
	intensity = 0;
	feature_id= 0;
	cluster_id = 0;
	cluster_size = 0;
}
DataPoint::DataPoint(const DataPoint &copyin) : GridElement(copyin)
{
	rt = copyin.rt;
	mz = copyin.mz;
	feature_id=copyin.feature_id;
	intensity = copyin.intensity;
	cluster_id = copyin.cluster_id;
	cluster_size = copyin.cluster_size;
}
DataPoint::DataPoint(DoubleReal rt_, DoubleReal mz_, DoubleReal intensity_, Int feature_id_) :GridElement(rt_,mz_)
{
	feature_id=feature_id_;
	intensity = intensity_;
	cluster_id = 0;
	cluster_size = 0;
}

DataPoint& DataPoint::operator=(const DataPoint &rhs)
{
	this->rt = rhs.rt;
	this->mz = rhs.mz;
	this->intensity = rhs.intensity;
	this->feature_id=rhs.feature_id;
	this->cluster_id = rhs.cluster_id;
	this->cluster_size = rhs.cluster_size;
	return *this;
}

bool DataPoint::operator==(const DataPoint &rhs) const
{
	if( this->feature_id != rhs.feature_id) return false;
	if( this->rt != rhs.rt) return false;
	if( this->mz != rhs.mz) return false;
	if( this->intensity != rhs.intensity) return false;
	//if( this->cluster_id != rhs.cluster_id) return false;
	//if( this->cluster_size != rhs.cluster_size) return false;
	return true;
}


bool DataPoint::operator!=(const DataPoint &rhs) const
{
	return !(*this==rhs);
}
bool DataPoint::operator<(const DataPoint &rhs) const // required for built-in STL functions like sort
{
	if ( *this == rhs) return false;

	return this->feature_id < rhs.feature_id;
}

int DataPoint::getID()
{
	return feature_id;
}
}


