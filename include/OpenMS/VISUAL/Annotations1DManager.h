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

#ifndef OPENMS_VISUAL_ANNOTATIONS1DMANAGER_H
#define OPENMS_VISUAL_ANNOTATIONS1DMANAGER_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/KERNEL/PeakIndex.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/VISUAL/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Annotation1DPeakItem.h>

#include <QtCore/QObject>
#include <QtCore/QRectF>
#include <QtGui/QPainter>

#include <list>

namespace OpenMS
{
	class SpectrumCanvas;

	/** @brief A manager for 1D annotations on the Spectrum1DCanvas
			
			All methods dealing with Annotation1DItem objects are located in this class.
			Every Spectrum1DCanvas object has an instance of this class as a member which
			is responsible for drawing the items on the canvas, calculating their
			bounding boxes, adding and deleting them, as well as selecting them
			on the canvas.
			
			This class has a pointer to the SpectrumCanvas of the corresponding
			Spectrum1DCanvas which it needs in order to translate the data points
			(MZ / intensity) contained in the annotation items to canvas coordinates.
			This is done by SpectrumCanvas::dataToWidget_() which is protected in
			SpectrumCanvas. For this reason, this class is a friend of SpectrumCanvas
			and is thus able to call this method.
			
			The bounding boxes of the items are computed while they are drawn, because 
			the current QPainter device is needed in order to compute the boundings of
			items containing text.
			
			If you want to add a new subclass "Annotation1DMyCustomItem" of Annotation1DItem,
			tell drawAnnotations() how this item is drawn and how its bounding box is computed.
			
	*/
	class Annotations1DManager
	{
		
		public:
			
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;
			
			/// Constructor
			Annotations1DManager(SpectrumCanvas* canvas);
			
			/// Destructor
			virtual ~Annotations1DManager();
			
			/** @brief Returns a pointer to the item at @p pos on the specified @p layer, or 0, if not existent
					
					If more than one item's bounding box encloses @p pos , the one in the
					foreground is returned.
			*/
			Annotation1DItem* getItemAt(const LayerData& layer, const QPoint& pos) const;
			/// Selects the item at @p pos on the canvas, if it exists. Searches only in the specified @p layer
			void selectItemAt(const LayerData& layer, const QPoint& pos);
			/// Deselects the item at @p pos on the canvas, if it exists. Searches only in the specified @p layer
			void deselectItemAt(const LayerData& layer, const QPoint& pos);
			/// Draws all annotation items of the specified @p layer on @p painter
			void drawAnnotations(const LayerData& layer, QPainter& painter);
			/// Draws the @p bounding_box on @p painter
			void drawBoundingBox(const QRectF& bounding_box, QPainter& painter);
			/// Draws a distance annotation @p distance_item on @p painter			
			void drawDistanceItem(const LayerData& layer, Annotation1DDistanceItem* distance_item, QPainter& painter);		
			/// Draws an arbitrary text annotation @p text_item on @p painter
			void drawTextItem(Annotation1DTextItem* text_item, QPainter& painter);
			/// Draws a peak annotation @p peak_item on @p painter
			void drawPeakItem(Annotation1DPeakItem* peak_item, QPainter& painter);
			/// Selects all items in the specified @p layer
			void selectAll(const LayerData& layer);
			/// Deselects all items in the specified @p layer
			void deselectAll(const LayerData& layer);
			/// Removes the selected items from the @p layer
			void removeSelectedItems(const LayerData& layer);
			/// Adds a new distance item between @p peak_1 and @p peak_2 at the position determined by @p start_point and @p end_point to the @p layer
			void addDistanceItem(const LayerData& layer, const PeakIndex& peak_1, const PeakIndex& peak_2, const PointType& start_point, const PointType& end_point);
			/// Adds a new text item displaying @p text at @p position on @p layer
			void addTextItem(const LayerData& layer, const PointType& position, const String& text);
			/// Adds an annotation for @p peak displaying @p text at @p position on @p layer
			void addPeakItem(const LayerData& layer, const PointType& position, const PeakIndex& peak, const String& text);
			
		protected:

			/// The parent canvas. Needed for translation of data points to positions on the widget and vice versa
			SpectrumCanvas* canvas_;
			
	};
} // namespace OpenMS

#endif
