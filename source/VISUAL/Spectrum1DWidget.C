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

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>

#include <QtGui/QScrollBar>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum1DWidget::Spectrum1DWidget(const Param& preferences, QWidget* parent)
		: SpectrumWidget(preferences, parent),
			flipped_canvas_(0),
			has_second_canvas_(0)
	{
		//set the label mode for the axes  - side effect
		setCanvas_(new Spectrum1DCanvas(preferences, this));
		
		x_axis_->setLegend("m/z");
		x_axis_->setAllowShortNumbers(false);
		y_axis_->setLegend("Intensity");
		y_axis_->setAllowShortNumbers(true);
		y_axis_->setMinimumWidth(50);
	}
	
	void Spectrum1DWidget::setFlippedCanvas(Spectrum1DCanvas* flipped_canvas)
	{
		flipped_canvas_ = flipped_canvas;
		flipped_canvas_->setFlippedVertically(true);
		flipped_canvas_->setSpectrumWidget(this);
		// make sure canvasses don't overlap:
		canvas_->setMinimumHeight(0);
		flipped_canvas_->setMinimumHeight(0);
		
		grid_->removeWidget(x_axis_);
		grid_->removeWidget(x_scrollbar_);
		grid_->addWidget(x_axis_, 2, 2);
		grid_->addWidget(x_scrollbar_, 3, 2);
		grid_->addWidget(flipped_canvas_, 1, 2);
		flipped_y_axis_ = new AxisWidget(AxisWidget::LEFT,"Intensity",this);
		flipped_y_axis_->setInverseOrientation(true);
		flipped_y_axis_->setAllowShortNumbers(true);
		flipped_y_axis_->setMinimumWidth(50);
		grid_->addWidget(flipped_y_axis_, 1, 1);
		
		connect(flippedCanvas(), SIGNAL(visibleAreaChanged(DRange<2>)), canvas(), SLOT(setVisibleArea(DRange<2>)));
		connect(canvas(), SIGNAL(visibleAreaChanged(DRange<2>)), flippedCanvas(), SLOT(setVisibleArea(DRange<2>)));
		connect(flippedCanvas(), SIGNAL(zoomLevelAdded()), canvas(), SLOT(zoomAdd()));
		connect(canvas(), SIGNAL(zoomLevelAdded()), flippedCanvas(), SLOT(zoomAdd()));
		connect(flippedCanvas(), SIGNAL(zoomedForward()), canvas(), SLOT(zoomForward()));
		connect(canvas(), SIGNAL(zoomedForward()), flippedCanvas(), SLOT(zoomForward()));
		connect(flippedCanvas(), SIGNAL(zoomedBack()), canvas(), SLOT(zoomBack()));
		connect(canvas(), SIGNAL(zoomedBack()), flippedCanvas(), SLOT(zoomBack()));
		
		connect(flippedCanvas(), SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(updateAxes()));
		connect(flippedCanvas(), SIGNAL(recalculateAxes()), this, SLOT(updateAxes()));
		connect(flippedCanvas(), SIGNAL(changeLegendVisibility()), this, SLOT(changeLegendVisibility()));
		connect(flippedCanvas(), SIGNAL(updateHScrollbar(float,float,float,float)), this, SLOT(updateHScrollbar(float,float,float,float)));
		connect(flippedCanvas(), SIGNAL(updateVScrollbar(float,float,float,float)), this, SLOT(updateVScrollbar(float,float,float,float)));
		connect(x_scrollbar_, SIGNAL(valueChanged(int)), flippedCanvas(), SLOT(horizontalScrollBarChange(int)));
		connect(y_scrollbar_, SIGNAL(valueChanged(int)), flippedCanvas(), SLOT(verticalScrollBarChange(int)));
		connect(flippedCanvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)),this, SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)));
		connect(flippedCanvas(), SIGNAL(sendCursorStatus(double,double,double)), this, SIGNAL(sendCursorStatus(double,double,double)));
		
		has_second_canvas_ = true;
	}

	void Spectrum1DWidget::removeFlippedCanvas()
	{
		if (has_second_canvas_)
		{
			disconnect(flippedCanvas(), SIGNAL(visibleAreaChanged(DRange<2>)), canvas(), SLOT(setVisibleArea(DRange<2>)));
			disconnect(canvas(), SIGNAL(visibleAreaChanged(DRange<2>)), flippedCanvas(), SLOT(setVisibleArea(DRange<2>)));
			disconnect(flippedCanvas(), SIGNAL(zoomLevelAdded()), canvas(), SLOT(zoomAdd()));
			disconnect(canvas(), SIGNAL(zoomLevelAdded()), flippedCanvas(), SLOT(zoomAdd()));
			disconnect(flippedCanvas(), SIGNAL(zoomedForward()), canvas(), SLOT(zoomForward()));
			disconnect(canvas(), SIGNAL(zoomedForward()), flippedCanvas(), SLOT(zoomForward()));
			disconnect(flippedCanvas(), SIGNAL(zoomedBack()), canvas(), SLOT(zoomBack()));
			disconnect(canvas(), SIGNAL(zoomedBack()), flippedCanvas(), SLOT(zoomBack()));
			
			disconnect(flippedCanvas(), SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(updateAxes()));
			disconnect(flippedCanvas(), SIGNAL(recalculateAxes()), this, SLOT(updateAxes()));
			disconnect(flippedCanvas(), SIGNAL(changeLegendVisibility()), this, SLOT(changeLegendVisibility()));
			disconnect(flippedCanvas(), SIGNAL(updateHScrollbar(float,float,float,float)), this, SLOT(updateHScrollbar(float,float,float,float)));
			disconnect(flippedCanvas(), SIGNAL(updateVScrollbar(float,float,float,float)), this, SLOT(updateVScrollbar(float,float,float,float)));
			disconnect(x_scrollbar_, SIGNAL(valueChanged(int)), flippedCanvas(), SLOT(horizontalScrollBarChange(int)));
			disconnect(y_scrollbar_, SIGNAL(valueChanged(int)), flippedCanvas(), SLOT(verticalScrollBarChange(int)));
			disconnect(flippedCanvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)),this, SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)));
			disconnect(flippedCanvas(), SIGNAL(sendCursorStatus(double,double,double)), this, SIGNAL(sendCursorStatus(double,double,double)));
		
			grid_->removeWidget(flipped_canvas_);
			grid_->removeWidget(flipped_y_axis_);
			flipped_canvas_->close();
			flipped_y_axis_->close();
			delete flipped_canvas_;
			delete flipped_y_axis_;
			grid_->removeWidget(x_axis_);
			grid_->removeWidget(x_scrollbar_);
			grid_->addWidget(x_axis_, 1, 2);
			grid_->addWidget(x_scrollbar_, 2, 2);
			
			canvas()->setMinimumHeight(200);
			canvas()->setMirrorMode(false);
			has_second_canvas_ = false;
		}
	}
	
	void Spectrum1DWidget::calculateUnitedRanges(bool reset_zoom)
	{
		if (canvas() != 0 && flippedCanvas() != 0)
		{
			DRange<3> new_overall_range = canvas()->getDataRange().united(flippedCanvas()->getDataRange());
						
			DRange<3> canvas_range = new_overall_range;
			//don't change intensity range:
			canvas_range.setMinY(canvas()->getDataRange().minY());
			canvas_range.setMaxY(canvas()->getDataRange().maxY());
			
			DRange<3> fl_canvas_range = new_overall_range;
			//don't change intensity range:
			fl_canvas_range.setMinY(flippedCanvas()->getDataRange().minY());
			fl_canvas_range.setMaxY(flippedCanvas()->getDataRange().maxY());
			
			canvas()->setOverallDataRange(canvas_range);
			canvas()->setMirrorMode(true);
			flippedCanvas()->setOverallDataRange(fl_canvas_range);
			flippedCanvas()->setMirrorMode(true);
			
			if (reset_zoom)
			{
				canvas()->resetZoom();
				flippedCanvas()->resetZoom();
			}
		}
	}
	
	bool Spectrum1DWidget::hasSecondCanvas()
	{
		return has_second_canvas_;
	}
	
	void Spectrum1DWidget::setHasSecondCanvas(bool has_second_canvas)
	{
		has_second_canvas_ = has_second_canvas;
	}
	
	void Spectrum1DWidget::recalculateAxes_()
	{
		AxisWidget* mz_axis;
		AxisWidget* it_axis;
		
		//determine axes
		if (canvas()->isMzToXAxis())
		{
			mz_axis = x_axis_;
			it_axis = y_axis_;
		}
		else
		{
			mz_axis = y_axis_;
			it_axis = x_axis_;
		}
		
		recalculateAxes_(mz_axis, it_axis, canvas());
		
		if (has_second_canvas_)
		{
			recalculateAxes_(mz_axis, flipped_y_axis_, flippedCanvas());
		}
	}
	
	void Spectrum1DWidget::recalculateAxes_(AxisWidget* mz_axis, AxisWidget* it_axis, Spectrum1DCanvas* canvas)
	{
		// recalculate gridlines
		mz_axis->setAxisBounds(canvas->getVisibleArea().minX(), canvas->getVisibleArea().maxX());
		switch(canvas->getIntensityMode())
		{
			case SpectrumCanvas::IM_NONE:
				it_axis->setAxisBounds(canvas->getVisibleArea().minY(), canvas->getVisibleArea().maxY());
				break;
			case SpectrumCanvas::IM_PERCENTAGE:
				it_axis->setAxisBounds(canvas->getVisibleArea().minY() / canvas->getDataRange().maxY() * 100.0, canvas->getVisibleArea().maxY() / canvas->getDataRange().maxY() * 100.0);
				break;
			case SpectrumCanvas::IM_SNAP:
				it_axis->setAxisBounds(canvas->getVisibleArea().minY()/canvas->getSnapFactor(), canvas->getVisibleArea().maxY()/canvas->getSnapFactor());
				break;
			default:
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
	}
	
	Histogram<UInt,Real> Spectrum1DWidget::createIntensityDistribution_() const
	{
		Histogram<UInt,Real> tmp(canvas_->getCurrentMinIntensity(),canvas_->getCurrentMaxIntensity(),(canvas_->getCurrentMaxIntensity() - canvas_->getCurrentMinIntensity())/500.0);
	
		for (ExperimentType::SpectrumType::ConstIterator it = canvas_->getCurrentLayer().peaks[0].begin(); it != canvas_->getCurrentLayer().peaks[0].end(); ++it)
		{
			tmp.inc(it->getIntensity());
		}
		return tmp;
	}


	Histogram<UInt, Real> Spectrum1DWidget::createMetaDistribution_(const String& name) const
	{	
		Histogram<UInt,Real> tmp;
		const ExperimentType::SpectrumType::MetaDataArrays& meta_arrays = canvas_->getCurrentLayer().peaks[0].getMetaDataArrays();
		for(ExperimentType::SpectrumType::MetaDataArrays::const_iterator it = meta_arrays.begin(); it != meta_arrays.end(); it++)
		{
			if (it->getName()==name)
			{
				//determine min and max of the data
				Real min = numeric_limits<Real>::max(), max = -numeric_limits<Real>::max();
				for (UInt i=0; i<it->size(); ++i)
				{
					if ((*it)[i]<min) min = (*it)[i];
					if ((*it)[i]>max) max = (*it)[i];
				}
				if (min>=max) return tmp;
		
				//create histogram
				tmp.reset(min,max,(max-min)/500.0);
				for (UInt i=0; i<it->size(); ++i)
				{
					tmp.inc((*it)[i]);
				}
			}
		}
		//fallback if no array with that name exists
		return tmp;
	}
	
	Spectrum1DWidget::~Spectrum1DWidget()
	{
		
	}

	void Spectrum1DWidget::showGoToDialog()
	{
	  Spectrum1DGoToDialog goto_dialog(this);
	  goto_dialog.setRange(canvas()->getDataRange().minX(),canvas()->getDataRange().maxX());
	  if (goto_dialog.exec())
	  {
	  	canvas()->setVisibleArea(SpectrumCanvas::AreaType(goto_dialog.getMin(),0,goto_dialog.getMax(),0));
		}
	}

} //namespace


