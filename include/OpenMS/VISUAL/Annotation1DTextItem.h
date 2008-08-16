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

#ifndef OPENMS_VISUAL_ANNOTATION1DTEXTITEM_H
#define OPENMS_VISUAL_ANNOTATION1DTEXTITEM_H

#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtCore/QRectF>

namespace OpenMS
{
	/** @brief An annotation item which represents an arbitrary text on the canvas.
			@see Annotation1DItem
	*/
	class Annotation1DTextItem
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
			Annotation1DTextItem(const PointType& position, const String& text);
			/// Copy constructor
			Annotation1DTextItem(const Annotation1DTextItem& rhs);
			/// Destructor
			virtual ~Annotation1DTextItem();
			/// Returns the current bounding box of this item on the canvas where it has last been drawn
			virtual const QRectF& boundingBox() const;
			/// Sets the bounding_box_ for this item
			virtual void setBoundingBox(const QRectF& bbox);
			/// Sets whether this item is currently selected on the canvas or not
			virtual void setSelected(bool selected);
			/// Returns true if this item is currently selected on the canvas, else false
			virtual const bool isSelected() const;
			/// Sets the text
 			virtual void setText(const String& text);
 			/// Returns the text
			virtual const String& getText() const;
			/// Sets the position of the item (in MZ / intensity coordinates)
			void setPosition(const PointType& position);
			/// Returns the position of the item (in MZ / intensity coordinates)
 			const PointType& getPosition() const;
			
		protected:
		
			/// The current bounding box of this item on the canvas where it has last been drawn
			QRectF bounding_box_;
			/// Determines whether this item is currently selected on the canvas
			bool is_selected_;
			/// The position of the item (in MZ / intensity coordinates)
			PointType position_;
			/// The text to be drawn on the canvas
			String text_;
			
	};
} // namespace OpenMS

#endif
