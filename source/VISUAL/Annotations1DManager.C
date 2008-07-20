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

#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/Annotations1DManager.h>

#include <iostream>

namespace OpenMS
{	
	
	Annotations1DManager::Annotations1DManager(SpectrumCanvas* canvas)
	{
		canvas_ = canvas;
	}
	
	Annotation1DItem* Annotations1DManager::getItemAt(const LayerData& layer, const QPoint& pos) const
	{
		for (LayerData::Ann1DConstIterator it = layer.annotations_1d_.begin(); it != layer.annotations_1d_.end(); ++it)
		{
			if ((*it)->boundingBox().contains(pos))
			{
				return *it;
			}
		}
		return 0;
	}
	
	void Annotations1DManager::selectItemAt(const LayerData& layer, const QPoint& pos)
	{
		Annotation1DItem* item = getItemAt(layer, pos);
		if (item != 0)
		{
			item->setSelected(true);
		}
	}
	
	void Annotations1DManager::deselectItemAt(const LayerData& layer, const QPoint& pos)
	{
		Annotation1DItem* item = getItemAt(layer, pos);
		if (item != 0)
		{
			item->setSelected(false);
		}
	}
	
	
	void Annotations1DManager::drawAnnotations(const LayerData& layer, QPainter& painter)
	{
		for (LayerData::Ann1DConstIterator it = layer.annotations_1d_.begin(); it != layer.annotations_1d_.end(); ++it)
		{
			Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(*it);
			if (distance_item)
			{
				drawDistanceItem(distance_item, painter);
			}
		}
	}
	
	void Annotations1DManager::drawDistanceItem(Annotation1DDistanceItem* item, QPainter& painter)
	{
		//translate mz/intensity to pixel coordinates
		QPoint start_point, end_point;
		canvas_->dataToWidget_(item->getStartPoint().getX(), item->getStartPoint().getY(), start_point);
		canvas_->dataToWidget_(item->getEndPoint().getX(), item->getEndPoint().getY(), end_point);
		
		// compute bounding box of item on the specified painter
		QRectF bbox(QPointF(start_point.x(), start_point.y()), QPointF(end_point.x(), end_point.y()+4)); // +4 for lower half of arrow heads
		// bbox must enclose distance text:
		const SpectrumCanvas::ExperimentType::PeakType& peak_1 = item->getStartPeak().getPeak(canvas_->getCurrentLayer().peaks);
		const SpectrumCanvas::ExperimentType::PeakType& peak_2 = item->getEndPeak().getPeak(canvas_->getCurrentLayer().peaks);
		QString distance_string = QString("%1").arg(peak_2.getMZ()-peak_1.getMZ());
		// find out how much additional space is needed for the text:
		QRectF text_boundings = painter.boundingRect(QRectF(), Qt::AlignCenter, distance_string);
		bbox.setTop(bbox.top() - text_boundings.height());
		
		item->setBoundingBox(bbox);
		
		if (item->isSelected())
		{
			painter.setPen(Qt::green);
			
			// draw additional filled rectangles to highlight bounding box of selected item
			painter.fillRect(bbox.topLeft().x()-3, bbox.topLeft().y()-3, 3, 3, Qt::green);
			painter.fillRect(bbox.topRight().x(), bbox.topRight().y()-3, 3, 3, Qt::green);
			painter.fillRect(bbox.bottomRight().x(), bbox.bottomRight().y(), 3, 3, Qt::green);
			painter.fillRect(bbox.bottomLeft().x()-3, bbox.bottomLeft().y(), 3, 3, Qt::green);

		}
		else
		{
			painter.setPen(Qt::darkGreen);
		}
		
		// draw line
		painter.drawLine(start_point, end_point);
		// draw arrow heads and the ends
		painter.drawLine(start_point, QPoint(start_point.x()+5, start_point.y()-4));
		painter.drawLine(start_point, QPoint(start_point.x()+5, start_point.y()+4));
		painter.drawLine(end_point, QPoint(end_point.x()-5, end_point.y()-4));
		painter.drawLine(end_point, QPoint(end_point.x()-5, end_point.y()+4));
		// draw distance text
		painter.drawText(bbox, Qt::AlignHCenter, distance_string);

	}
	
	void Annotations1DManager::selectAll(const LayerData& layer)
	{
		for (LayerData::Ann1DIterator it = layer.annotations_1d_.begin(); it != layer.annotations_1d_.end(); ++it)
		{
			(*it)->setSelected(true);
		}
	}
	
	void Annotations1DManager::deselectAll(const LayerData& layer)
	{
		for (LayerData::Ann1DIterator it = layer.annotations_1d_.begin(); it != layer.annotations_1d_.end(); ++it)
		{
			(*it)->setSelected(false);
		}
	}
	
	// Predicate needed by removeSelectedItems()
	bool annotationItemSelected(const Annotation1DItem* item) { return item->isSelected(); }
	
	void Annotations1DManager::removeSelectedItems(const LayerData& layer)
	{
		layer.annotations_1d_.remove_if(annotationItemSelected);
	}
	
	void Annotations1DManager::addDistanceItem(const LayerData& layer, const PeakIndex& peak_1, const PeakIndex& peak_2, const PointType& start_point, const PointType& end_point)
	{
		Annotation1DItem* new_item = new Annotation1DDistanceItem(peak_1, peak_2, start_point, end_point);
		layer.annotations_1d_.push_front(new_item);
	}
	
}//Namespace




