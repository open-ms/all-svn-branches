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
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/KERNEL/PeakIndex.h>

#include <QtCore/QRectF>

namespace OpenMS
{
	/** @brief An annotation item which represents a measured distance between two peaks.
			@see Annotation1DItem
	*/
	class Annotation1DDistanceItem
		: public Annotation1DItem
	{
		
		public:
			
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Intensity type
			typedef Real IntensityType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;

			/// Constructor
			Annotation1DDistanceItem(const PeakIndex& start_peak, const PeakIndex& end_peak,
				const Annotation1DDistanceItem::PointType start_point, const Annotation1DDistanceItem::PointType end_point);
			
			/// Returns the current bounding box of this item on the canvas where it has last been drawn
			const QRectF& boundingBox() const;
			/// Sets the bounding_box_ for this item
			void setBoundingBox(const QRectF& bbox);
			/// Sets whether this item is currently selected on the canvas or not
			void setSelected(bool selected);
			/// Returns true if this item is currently selected on the canvas, else false
			const bool isSelected() const;
			/// Sets the start peak index of the measurement
			void setStartPeak(const PeakIndex& start_peak);
			/// Sets the end peak index of the measurement
			void setEndPeak(const PeakIndex& end_peak);
			/// Sets the start point of the measured distance line
			void setStartPoint(const PointType& start);
			/// Sets the peak index of the end peak of the measurement
			void setEndPoint(const PointType& end);
			/// Returns the start peak index of the measurement
			const PeakIndex& getStartPeak() const;
			/// Returns the end peak index of the measurement
			const PeakIndex& getEndPeak() const;
			/// Returns the start point as (MZ,intensity)
			const PointType& getStartPoint() const;
			/// Returns the end point as (MZ,intensity)
			const PointType& getEndPoint() const;
				
		protected:
		
			/// The current bounding box of this item on the canvas where it has last been drawn
			QRectF bounding_box_;
			/// Determines whether this item is currently selected on the canvas
			bool is_selected_;
			/// PeakIndex of the start peak of the measurement
			PeakIndex start_peak_;
			/// PeakIndex of the end peak of the measurement
			PeakIndex end_peak_;
			/// The start point of the measured distance line
			PointType start_point_;
			/// The end point of the measured distance line
			PointType end_point_;
			
	};
} // namespace OpenMS

#endif
