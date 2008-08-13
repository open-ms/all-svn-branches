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

#include <OpenMS/VISUAL/Annotation1DPeakItem.h>

namespace OpenMS
{	

	Annotation1DPeakItem::Annotation1DPeakItem(const PointType& position, const PeakIndex& peak, const String& text)
		: Annotation1DItem(),
			bounding_box_(),
			is_selected_(true),
			position_(position),
			peak_(peak),
			text_(text)
	{
	}
	
	Annotation1DPeakItem::Annotation1DPeakItem(const Annotation1DPeakItem& rhs)
		: Annotation1DItem()
	{
		bounding_box_ = rhs.boundingBox();
		is_selected_ = rhs.isSelected();
		position_ = rhs.getPosition();
		peak_ = rhs.getPeak();
		text_ = rhs.getText();
	}
	
	Annotation1DPeakItem::~Annotation1DPeakItem()
	{
	}
	
	const QRectF& Annotation1DPeakItem::boundingBox() const
	{
		return bounding_box_;
	}
	
	void Annotation1DPeakItem::setBoundingBox(const QRectF& bbox)
	{
		bounding_box_ = bbox;
	}
	
	void Annotation1DPeakItem::setSelected(bool selected)
	{
		is_selected_ = selected;
	}
	
	const bool Annotation1DPeakItem::isSelected() const
	{
		return is_selected_;
	}
	
	void Annotation1DPeakItem::setPosition(const Annotation1DPeakItem::PointType& position)
	{
		position_ = position;
	}
	
 	const Annotation1DPeakItem::PointType& Annotation1DPeakItem::getPosition() const
 	{
 		return position_;
 	}
	
	void Annotation1DPeakItem::setPeak(const PeakIndex& peak)
	{
		peak_ = peak;
	}
	
	const PeakIndex& Annotation1DPeakItem::getPeak() const
	{
		return peak_;
	}
	
	void Annotation1DPeakItem::setText(const String& text)
	{
		text_ = text;
	}
	
	const String& Annotation1DPeakItem::getText() const
	{
		return text_;
	}
	
}//Namespace




