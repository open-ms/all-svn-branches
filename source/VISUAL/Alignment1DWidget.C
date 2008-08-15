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
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

#include <QtCore/QPoint>
#include <QtGui/QPainter>
#include <QtGui/QWidget>
#include <QtGui/QPaintEvent>

#include <iostream>

#include <OpenMS/VISUAL/Alignment1DWidget.h>

using namespace std;

namespace OpenMS
{
	
	Alignment1DWidget::Alignment1DWidget(Spectrum1DWidget* spectrum_widget, QWidget* parent)
		: QWidget(parent),
			spectrum_widget_(spectrum_widget),
			alignment_is_set_(false)
	{
		setMinimumHeight(1);
		setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
	}
	
	Alignment1DWidget::~Alignment1DWidget()
	{
	}

	SpectrumCanvas* Alignment1DWidget::upperCanvas_()
	{
		return qobject_cast<SpectrumCanvas*>(spectrum_widget_->canvas());
	}
	
	SpectrumCanvas* Alignment1DWidget::lowerCanvas_()
	{
		return qobject_cast<SpectrumCanvas*>(spectrum_widget_->flippedCanvas());
	}
	
	void Alignment1DWidget::setAlignmentLines(const std::vector<std::pair<DoubleReal, DoubleReal > >& alignment_lines)
	{
		alignment_lines_ = alignment_lines;
		setMinimumHeight(10);
		update();
		alignment_is_set_ = true;
	}
	
	void Alignment1DWidget::clearAlignmentLines()
	{
		alignment_lines_.clear();
		setMinimumHeight(1);
		update();
		alignment_is_set_ = false;
	}
	
	bool Alignment1DWidget::alignmentIsSet()
	{
		return alignment_is_set_;
	}
	
	void Alignment1DWidget::paintEvent(QPaintEvent* e)
	{
		QPainter painter;
		painter.begin(this);
		
		if (!alignment_is_set_)
		{
			painter.setPen(Qt::black);
			painter.drawLine(0,0,width(),0);
		}
		else
		{
			painter.setPen(Qt::black);
			int h = height();
			int w = width();
			//draw x-axes
			painter.drawLine(0,0,w,0);
			painter.drawLine(0,h-1,w,h-1);
			//draw peak-connecting lines between the two spectra
			painter.setPen(Qt::red);
			QPoint begin_p, end_p;
			double dummy = 0.0;
			
			for (UInt i = 0; i < alignment_lines_.size(); ++i)
			{
				upperCanvas_()->dataToWidget_(alignment_lines_[i].first, dummy, begin_p);
				lowerCanvas_()->dataToWidget_(alignment_lines_[i].second, dummy, end_p);
				painter.drawLine(begin_p.x(), 0, end_p.x(), h-1);
			}
		}
		
		painter.end();
		e->accept();
	}
}//Namespace




