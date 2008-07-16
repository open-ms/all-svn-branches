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

#include <OpenMS/VISUAL/Annotations1D.h>

#include <algorithm>

namespace OpenMS
{	
	
	Annotations1D::Annotations1D()
	{
		empty_item_ = new Annotation1DEmptyItem();
	}
	
	Annotation1DItem* Annotations1D::getAnnotationAt(PointType& p)
	{
		return getAnnotationAt(p.getX(), p.getY());
	}
	
	Annotation1DItem* Annotations1D::getAnnotationAt(CoordinateType x, CoordinateType y)
	{
		// try to find an item whose bounding box encloses (x,y)
		for(Ann1DIterator it = annotation_items_.begin(); it != annotation_items_.end(); it++)
		{
			if ((*it)->boundingBox().encloses(x,y))
			{
				return *it;
			}
		}
		// if no such item exists, return the empty item
		return empty_item_;
	}
	
	UInt Annotations1D::size()
	{
		return annotation_items_.size();
	}
	
	bool Annotations1D::isEmpty()
	{
		return annotation_items_.empty();
	}
	
	void Annotations1D::add(Annotation1DItem* item)
	{
		annotation_items_.push_back(item);
	}
	
	bool Annotations1D::remove(Annotation1DItem* item)
	{
		Ann1DIterator it = find(annotation_items_.begin(), annotation_items_.end(), item);
		if (it == annotation_items_.end())
		{
			return false;
		}
		annotation_items_.erase(it);
		return true;
	}
	
	void Annotations1D::clear()
	{
		annotation_items_.clear();
	}
	
}//Namespace




