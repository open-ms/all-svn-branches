// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Janett Köppen $
// $Authors: Janett Köppen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureFileOptions.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
	FeatureFileOptions::FeatureFileOptions()
		: loadConvexhull_(true),
			loadSubordinates_(true)
	{
	}
	
	FeatureFileOptions::~FeatureFileOptions()
	{
	}
	
	void FeatureFileOptions::setLoadConvexHull(bool convex)
	{
		loadConvexhull_ = convex;
	}
	
	bool FeatureFileOptions::getLoadConvexHull() const
	{
		return loadConvexhull_;
	}

	void FeatureFileOptions::setLoadSubordinates(bool sub)
	{
		loadSubordinates_ = sub;
	}

	bool FeatureFileOptions::getLoadSubordinates() const
	{
		return loadSubordinates_;
	}
	
	
} // namespace OpenMS
