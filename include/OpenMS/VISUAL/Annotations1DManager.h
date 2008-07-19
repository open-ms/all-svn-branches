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

#include <QtCore/QObject>
#include <QtCore/QRectF>
#include <QtGui/QPainter>

#include <list>

namespace OpenMS
{
	class Spectrum1DCanvas;

	/** blablabla
	*/
	class Annotations1DManager
	{
		
		public:
			
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;
			
			/// Constructor
			Annotations1DManager(Spectrum1DCanvas* canvas);
			
			/// Selects the item at @p pos on the canvas, if it exists. Searches only in the specified @p layer
			void selectItemAt(const LayerData& layer, const QPoint& pos);
			/// Draws all annotation items of the specified @p layer on @p painter
			void drawAnnotations(const LayerData& layer, QPainter& painter);
			/// Draws a distance annotation @p distance_item on @p painter
			void drawDistanceItem(Annotation1DDistanceItem* distance_item, QPainter& painter);		
			/// Selects all items in the specified @p layer
			void selectAll(const LayerData& layer);
			/// Deselects all items in the specified @p layer
			void deselectAll(const LayerData& layer);
			/// Removes the selected items from the @p layer
			void removeSelectedItems(const LayerData& layer);
			/// Adds a new distance item between @p peak_1 and @p peak_2 at the position determined by @p start_point and @p end_point to the @p layer
			void addDistanceItem(const LayerData& layer, const PeakIndex& peak_1, const PeakIndex& peak_2, const PointType& start_point, const PointType& end_point);
			
		protected:
			
			/// The parent canvas. Needed for translation of data points to positions on the widget and vice versa.
			Spectrum1DCanvas* canvas_1d_;
			
	};
} // namespace OpenMS

#endif
