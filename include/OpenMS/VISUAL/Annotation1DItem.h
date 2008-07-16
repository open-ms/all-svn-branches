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

#ifndef OPENMS_VISUAL_ANNOTATION1DITEM_H
#define OPENMS_VISUAL_ANNOTATION1DITEM_H

#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

namespace OpenMS
{
	/** blablabla Abstract class blabla
	*/
	class Annotation1DItem
	{
		
		public:
			
			/// Type of the bounding boxes
			typedef DBoundingBox<2> BoundingBoxType;
			
			/// Returns the bounding box of the item
			virtual const BoundingBoxType boundingBox() = 0;
			
			/// Returns true if this item is currently selected on the canvas, else false
			virtual const bool isSelected() = 0;
			
			/// Sets whether this item is currently selected on the canvas or not
			virtual void setSelected(bool selected) = 0;
			
	};
} // namespace OpenMS

#endif
