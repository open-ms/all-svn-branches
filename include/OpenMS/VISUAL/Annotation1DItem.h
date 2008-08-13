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

#include <QtCore/QRectF>

namespace OpenMS
{
	/** @brief An abstract class acting as an interface for the different 1d annotation items.
	
			This is an abstract class which acts as an interface between its polymorphic
			subclasses and all containers and methods that contain or handle Annotation1DItem
			objects.
			
			If you want to add new kinds of annotation items, inherit this class,
			implement the pure virtual methods, and add everything else the item should
			have or be capable of. Then tell the Annotations1DManager how items of
			this new class are drawn and how their bounding boxes are computed
			
	*/
	class Annotation1DItem
	{
		
		public:
			
			/// Destructor
			virtual ~Annotation1DItem();
			
			/// Returns the current bounding box of this item on the canvas where it has last been drawn
			virtual const QRectF& boundingBox() const = 0;
			
			/// Returns true if this item is currently selected on the canvas, else false
			virtual const bool isSelected() const = 0;
			
			/// Sets whether this item is currently selected on the canvas or not
			virtual void setSelected(bool selected) = 0;
			
	};
} // namespace OpenMS

#endif
