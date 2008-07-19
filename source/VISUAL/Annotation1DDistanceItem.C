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

#include <OpenMS/VISUAL/Annotation1DDistanceItem.h>

namespace OpenMS
{	

	Annotation1DDistanceItem::Annotation1DDistanceItem(const PeakIndex& start_peak, const PeakIndex& end_peak,
		const Annotation1DDistanceItem::PointType start_point, const Annotation1DDistanceItem::PointType end_point)
		: bounding_box_(),
			is_selected_(true),
			start_peak_(start_peak),
			end_peak_(end_peak),
			start_point_(start_point),
			end_point_(end_point)
	{
	}
	
	const QRectF& Annotation1DDistanceItem::boundingBox() const
	{
		return bounding_box_;
	}
	
	void Annotation1DDistanceItem::setBoundingBox(const QRectF& bbox)
	{
		bounding_box_ = bbox;
	}
	
	void Annotation1DDistanceItem::setSelected(bool selected)
	{
		is_selected_ = selected;
	}
	
	const bool Annotation1DDistanceItem::isSelected() const
	{
		return is_selected_;
	}
	
	void Annotation1DDistanceItem::setStartPeak(PeakIndex& start_peak)
	{
		start_peak_ = start_peak;
	}
	
	void Annotation1DDistanceItem::setEndPeak(PeakIndex& end_peak)
	{
		end_peak_ = end_peak;
	}
	
	const PeakIndex& Annotation1DDistanceItem::getStartPeak() const
	{
		return start_peak_;
	}
	
	const PeakIndex& Annotation1DDistanceItem::getEndPeak() const
	{
		return end_peak_;
	}
	
	void Annotation1DDistanceItem::setStartPoint(PointType& p)
	{
		start_point_ = p;
	}

	void Annotation1DDistanceItem::setEndPoint(PointType& p)
	{
		end_point_ = p;
	}
	
	const Annotation1DDistanceItem::PointType& Annotation1DDistanceItem::getStartPoint() const
	{
		return start_point_;
	}
	
	const Annotation1DDistanceItem::PointType& Annotation1DDistanceItem::getEndPoint() const
	{
		return end_point_;
	}
	
}//Namespace




