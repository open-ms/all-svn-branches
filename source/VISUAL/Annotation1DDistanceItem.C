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

	Annotation1DDistanceItem::Annotation1DDistanceItem(PeakIndex& start_peak, PeakIndex& end_peak,
		Annotation1DDistanceItem::PointType start_point, Annotation1DDistanceItem::PointType end_point)
		: is_selected_(true),
			start_peak_(start_peak),
			end_peak_(end_peak),
			start_point_(start_point),
			end_point_(end_point)
	{
	}
	
	const Annotation1DDistanceItem::BoundingBoxType Annotation1DDistanceItem::boundingBox()
	{
		const Annotation1DDistanceItem::PointType min_point(start_point_.getX(), start_point_.getY());
		const Annotation1DDistanceItem::PointType max_point(end_point_.getX(), end_point_.getY());
		const Annotation1DDistanceItem::BoundingBoxType box(min_point, max_point);
		
		return box;
	}
	
	void Annotation1DDistanceItem::setSelected(bool selected)
	{
		is_selected_ = selected;
	}
	
	const bool Annotation1DDistanceItem::isSelected()
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
	
	PeakIndex& Annotation1DDistanceItem::getStartPeak()
	{
		return start_peak_;
	}
	
	PeakIndex& Annotation1DDistanceItem::getEndPeak()
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
	
	const Annotation1DDistanceItem::PointType& Annotation1DDistanceItem::getStartPoint()
	{
		return start_point_;
	}
	
	const Annotation1DDistanceItem::PointType& Annotation1DDistanceItem::getEndPoint()
	{
		return end_point_;
	}
	
}//Namespace




