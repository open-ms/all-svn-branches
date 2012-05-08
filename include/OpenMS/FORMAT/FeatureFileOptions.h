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

#ifndef OPENMS_FORMAT_FEATUREFILEOPTIONS_H
#define OPENMS_FORMAT_FEATUREFILEOPTIONS_H

#include <OpenMS/DATASTRUCTURES/DRange.h> 

#include <vector>

namespace OpenMS
{
	/**
		@brief Options for loading feature files.
	*/
	class OPENMS_DLLAPI FeatureFileOptions
	{
	public:
		///Default constructor
		FeatureFileOptions();
		///Destructor
		~FeatureFileOptions();

		///@name convex hull option
		///sets whether or not to load convex hull
		void setLoadConvexHull(bool convex);
		///returns whether or not to load convex hull
		bool getLoadConvexHull() const;

		///@name subordinate option
		///sets whether or not load subordinates
		void setLoadSubordinates(bool sub);
		///returns whether or not to load subordinates
		bool getLoadSubordinates() const;
		

	private:
		bool loadConvexhull_;
		bool loadSubordinates_;

		

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_FEATUREFILEOPTIONS_H
