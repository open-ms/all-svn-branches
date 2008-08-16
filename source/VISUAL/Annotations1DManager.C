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
#include <OpenMS/VISUAL/Annotations1DContainer.h>

#include <OpenMS/VISUAL/Annotations1DManager.h>

#include <iostream>

namespace OpenMS
{	
	
	Annotations1DManager::Annotations1DManager(SpectrumCanvas* canvas)
		: canvas_(canvas)
	{
	}
	
	Annotations1DManager::~Annotations1DManager()
	{
	}
	
	void Annotations1DManager::setPen(const QPen& pen)
	{
		pen_ = pen;
		selected_pen_ = pen;
		
		//make selected items a little brighter
		int sel_red = pen.color().red() + 50;
		int sel_green = pen.color().green() + 50;
		int sel_blue = pen.color().blue() + 50;
		//check if rgb out of bounds
		sel_red = sel_red > 255 ? 255 : sel_red;
		sel_green = sel_green > 255 ? 255 : sel_green;
		sel_blue = sel_blue > 255 ? 255 : sel_blue;
		
		selected_pen_.setColor(QColor(sel_red, sel_green, sel_blue));
	}
	
	Annotation1DItem* Annotations1DManager::getItemAt(const LayerData& layer, const QPoint& pos) const
	{
		for (Annotations1DContainer::ConstIterator it = layer.annotations_1d.begin(); it != layer.annotations_1d.end(); ++it)
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
		for (Annotations1DContainer::ConstIterator it = layer.annotations_1d.begin(); it != layer.annotations_1d.end(); ++it)
		{
			Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(*it);
			if (distance_item)
			{
				drawDistanceItem(distance_item, painter);
				continue;
			}
			
			Annotation1DTextItem* text_item = dynamic_cast<Annotation1DTextItem*>(*it);
			if (text_item)
			{
				drawTextItem(text_item, painter);
				continue;
			}
			
			Annotation1DPeakItem* peak_item = dynamic_cast<Annotation1DPeakItem*>(*it);
			if (peak_item)
			{
				drawPeakItem(peak_item, painter);
				continue;
			}
		}
	}
	
	void Annotations1DManager::drawBoundingBox(const QRectF& bounding_box, QPainter& painter)
	{
		// draw additional filled rectangles to highlight bounding box of selected distance_item
		painter.fillRect(bounding_box.topLeft().x()-3, bounding_box.topLeft().y()-3, 3, 3, selected_pen_.color());
		painter.fillRect(bounding_box.topRight().x(), bounding_box.topRight().y()-3, 3, 3, selected_pen_.color());
		painter.fillRect(bounding_box.bottomRight().x(), bounding_box.bottomRight().y(), 3, 3, selected_pen_.color());
		painter.fillRect(bounding_box.bottomLeft().x()-3, bounding_box.bottomLeft().y(), 3, 3, selected_pen_.color());
	}
	
	void Annotations1DManager::drawDistanceItem(Annotation1DDistanceItem* distance_item, QPainter& painter)
	{
		//translate mz/intensity to pixel coordinates
		QPoint start_point, end_point;
		canvas_->dataToWidget_(distance_item->getStartPoint().getX(), distance_item->getStartPoint().getY(), start_point);
		canvas_->dataToWidget_(distance_item->getEndPoint().getX(), distance_item->getEndPoint().getY(), end_point);
		
		// compute bounding box of distance_item on the specified painter
		QRectF bbox(QPointF(start_point.x(), start_point.y()), QPointF(end_point.x(), end_point.y()+4)); // +4 for lower half of arrow heads
		
		// find out how much additional space is needed for the text:
		QRectF text_boundings = painter.boundingRect(QRectF(), Qt::AlignCenter, distance_item->getText().toQString());
		bbox.setTop(bbox.top() - text_boundings.height());
		// if text doesn't fit between peaks, enlarge bounding box:
		if (text_boundings.width() > bbox.width())
		{
			float additional_space = (text_boundings.width() - bbox.width()) / 2;
			bbox.setLeft(bbox.left() - additional_space);
			bbox.setRight(bbox.right() + additional_space);
		}
		
		distance_item->setBoundingBox(bbox);
		
		if (distance_item->isSelected())
		{
			painter.setPen(selected_pen_);
			drawBoundingBox(bbox, painter);
		}
		else
		{
			painter.setPen(pen_);
		}
		
		// draw line
		painter.drawLine(start_point, end_point);
		// draw arrow heads and the ends
		painter.drawLine(start_point, QPoint(start_point.x()+5, start_point.y()-4));
		painter.drawLine(start_point, QPoint(start_point.x()+5, start_point.y()+4));
		painter.drawLine(end_point, QPoint(end_point.x()-5, end_point.y()-4));
		painter.drawLine(end_point, QPoint(end_point.x()-5, end_point.y()+4));
		// draw distance text
		painter.drawText(bbox, Qt::AlignHCenter, distance_item->getText().toQString());
	}
	
	void Annotations1DManager::drawTextItem(Annotation1DTextItem* text_item, QPainter& painter)
	{
		//translate mz/intensity to pixel coordinates
		QPoint position;
		canvas_->dataToWidget_(text_item->getPosition().getX(), text_item->getPosition().getY(), position);
		
		// compute bounding box of text_item on the specified painter
		QRectF bbox = painter.boundingRect(QRectF(position, position), Qt::AlignCenter, text_item->getText().toQString());
		text_item->setBoundingBox(bbox);
		
		if (text_item->isSelected())
		{
			painter.setPen(selected_pen_);
			drawBoundingBox(bbox, painter);
		}
		else
		{
			painter.setPen(pen_);
		}
		
		painter.drawText(bbox, Qt::AlignCenter, text_item->getText().toQString());

	}
	
	void Annotations1DManager::drawPeakItem(Annotation1DPeakItem* peak_item, QPainter& painter)
	{
		//translate mz/intensity to pixel coordinates
		QPoint position;
		canvas_->dataToWidget_(peak_item->getPosition().getX(), peak_item->getPosition().getY(), position);
		
		// compute bounding box of text_item on the specified painter
		QRectF bbox = painter.boundingRect(QRectF(position, position), Qt::AlignCenter, peak_item->getText().toQString());
		// shift position, annotation should be next to the peak and not overlap it
		bbox.translate(bbox.width()/2.0 + 10.0, -15.0); 
		peak_item->setBoundingBox(bbox);
		
		if (peak_item->isSelected())
		{
			painter.setPen(selected_pen_);
			drawBoundingBox(bbox, painter);
		}
		else
		{
			painter.setPen(pen_);
		}
		
		painter.drawText(bbox, Qt::AlignCenter, peak_item->getText().toQString());
		painter.drawLine(bbox.bottomLeft(), position);
	}
	
	void Annotations1DManager::selectAll(const LayerData& layer)
	{
		for (Annotations1DContainer::Iterator it = layer.annotations_1d.begin(); it != layer.annotations_1d.end(); ++it)
		{
			(*it)->setSelected(true);
		}
	}
	
	void Annotations1DManager::deselectAll(const LayerData& layer)
	{
		for (Annotations1DContainer::Iterator it = layer.annotations_1d.begin(); it != layer.annotations_1d.end(); ++it)
		{
			(*it)->setSelected(false);
		}
	}
		
	void Annotations1DManager::removeSelectedItems(const LayerData& layer)
	{
		for (Annotations1DContainer::Iterator it = layer.annotations_1d.begin(); it != layer.annotations_1d.end();)
		{
			if ((*it)->isSelected())
			{
				delete *it;
				it = layer.annotations_1d.erase(it);
			}
			else
			{
				++it;
			}
		}
	}
	
	void Annotations1DManager::addDistanceItem(const LayerData& layer, const String& text, const PointType& start_point, const PointType& end_point)
	{
		Annotation1DItem* new_item = new Annotation1DDistanceItem(text, start_point, end_point);
		layer.annotations_1d.push_front(new_item);
	}
	
	void Annotations1DManager::addTextItem(const LayerData& layer, const PointType& position, const String& text)
	{
		Annotation1DItem* new_item = new Annotation1DTextItem(position, text);
		layer.annotations_1d.push_front(new_item);
	}
	
	void Annotations1DManager::addPeakItem(const LayerData& layer, const PointType& position, const PeakIndex& peak, const String& text)
	{
		CoordinateType peak_height = peak.getPeak(layer.peaks).getIntensity();
		PointType actual_pos = position;

		if (position.getY() > peak_height)
		{
			actual_pos.setY(peak_height);
		}
		Annotation1DItem* new_item = new Annotation1DPeakItem(actual_pos, peak, text);
		layer.annotations_1d.push_front(new_item);
	}


}//Namespace
