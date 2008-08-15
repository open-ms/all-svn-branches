// -*- mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_VISUAL_SPECTRUM1DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM1DWIDGET_H

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

class QAction;

namespace OpenMS
{
	class Spectrum1DCanvas;
	class Alignment1DWidget;

	/**
		@brief Widget for visualization of several spectra
		
		The widget composes of a scoll bar, an AxisWidget and a Spectrum1DCanvas as central widget.
		
		@image html Spectrum1DWidget.png
		
		The example image shows %Spectrum1DWidget displaying a raw data layer and a peak data layer. 
		
		@ingroup SpectrumWidgets
	*/
	class Spectrum1DWidget 
		: public SpectrumWidget
	{
		Q_OBJECT
		
	public:
		/// Default constructor
		Spectrum1DWidget(const Param& preferences, QWidget* parent = 0);
		///Destructor
		virtual ~Spectrum1DWidget();
		
		/// This method is overwritten to make the class specific members accessable
		inline Spectrum1DCanvas* canvas()
		{
			return static_cast<Spectrum1DCanvas*>(canvas_);
		}
		
		/// Returns a pointer to the second (vertically flipped) canvas
		inline Spectrum1DCanvas* flippedCanvas()
		{
			return flipped_canvas_;
		}
		
		/// Returns a pointer to the alignment widget
		inline Alignment1DWidget* alignmentWidget()
		{
			return alignment_widget_;
		}
		
		/// Overwrites SpectrumWidget::setIntensityMode(); sets value also for second canvas, if existent
		virtual void setIntensityMode(SpectrumCanvas::IntensityModes mode);
		
		/// Sets the second (vertically flipped) canvas
		void setFlippedCanvas(Spectrum1DCanvas* flipped_canvas);
		
		/// Returns whether this widget has a second (vertically flipped) canvas
		bool hasSecondCanvas();
		
		/// Sets whether this widget has a second (vertically flipped) canvas
		void setHasSecondCanvas(bool has_second_canvas);
		
		/// Removes the second canvas
		void removeFlippedCanvas();
		
		/// Calculates and sets the united ranges for the canvas and the second (vertically flipped) canvas
		void calculateUnitedRanges(bool reset_zoom = false);
		
		/// Shows the alignment widget and sets its alignment_lines_ to @p alignment_lines
		void setAlignmentLines(const std::vector<std::pair<DoubleReal,DoubleReal > >& alignment_lines);
		
		/// Clears and hides the alignment widget
		void clearAlignmentLines();			
		
	signals:
		/// Is emitted whenever the visible area changes.		
		void visibleAreaChanged(double, double); 

	public slots:
		// Docu in base class
    virtual void showGoToDialog();

	protected:
		
		// Docu in base class
		virtual Math::Histogram<UInt, Real> createIntensityDistribution_() const;
		// Docu in base class
		virtual Math::Histogram<UInt, Real> createMetaDistribution_(const String& name) const;
		// Docu in base class
		virtual void recalculateAxes_();
		/// Recalculates the axis ticks of @p canvas from @p mz_axis and @p it_axis
		void recalculateAxes_(AxisWidget* mz_axis, AxisWidget* it_axis, Spectrum1DCanvas* canvas);
		
		/// The second (vertically flipped) canvas
		Spectrum1DCanvas* flipped_canvas_;
		/// The axis widget of the second canvas
		AxisWidget* flipped_y_axis_;
		/// The alignment widget between both canvasses
		Alignment1DWidget* alignment_widget_;
		/// Indicates whether this widget currently shows an additional (vertically flipped) canvas
		bool has_second_canvas_;
	
	};
} // namespace OpenMS

#endif
