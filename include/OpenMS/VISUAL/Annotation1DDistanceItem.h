// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ANNOTATION1DDISTANCEITEM_H
#define OPENMS_VISUAL_ANNOTATION1DDISTANCEITEM_H

#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

namespace OpenMS
{
	/** blablabla Abstract class blabla
	*/
	class Annotation1DDistanceItem
		: public Annotation1DItem
	{
		
		public:
			
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Type of the bounding boxes
			typedef DBoundingBox<2> BoundingBoxType;
			/// Intensity type
			typedef Real IntensityType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;

			/// Constructor
			Annotation1DDistanceItem(Annotation1DDistanceItem::PointType start, Annotation1DDistanceItem::PointType end);
			
			/// Returns the bounding box of the item
			const BoundingBoxType boundingBox();
			/// Sets whether this item is currently selected on the canvas or not
			void setSelected(bool selected);
			/// Returns true if this item is currently selected on the canvas, else false
			const bool isSelected();
			/// Sets the start point of the measured distance line
			void setStartPoint(PointType& start);
			/// Sets the peak index of the end peak of the measurement
			void setEndPoint(PointType& end);
			/// Returns the start point as (MZ,intensity)
			const PointType& getStartPoint();
			/// Returns the end point as (MZ,intensity)
			const PointType& getEndPoint();
				
		protected:
		
			/// Determines whether this item is currently selected on the canvas
			bool is_selected_;
			/// The start point of the measured distance line
			PointType start_point_;
			/// The end point of the measured distance line
			PointType end_point_;
			
	};
} // namespace OpenMS

#endif
