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

#ifndef OPENMS_VISUAL_SPECTRUM1DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM1DCANVAS_H

// STL
#include <vector>
#include <utility>

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/VISUAL/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/Annotations1DManager.h>

//QT
class QAction;

namespace OpenMS
{
	/**
		@brief Canvas for visualization of one or several spectra.
		
		@image html Spectrum1DCanvas.png
		
		The example image shows %Spectrum1DCanvas displaying a raw data layer and a peak data layer. 
		
		@ref Spectrum1DCanvas_Parameters are explained on a separate page.
		
		@ingroup SpectrumWidgets
	*/
	
	class Spectrum1DCanvas 
		: public SpectrumCanvas
	{
		Q_OBJECT
		
		public:
			/// Label modes (percentage or absolut) of x axis and y axis
			enum LabelMode
			{
				LM_XABSOLUTE_YABSOLUTE,
				LM_XPERCENT_YABSOLUTE,
				LM_XABSOLUTE_YPERCENT,
				LM_XPERCENT_YPERCENT
			};
			
			/// Default constructor
			Spectrum1DCanvas(const Param& preferences, QWidget* parent = 0);
			/// Destructor
			virtual ~Spectrum1DCanvas();
		
			///Enumerate all avaiable paint styles
			enum DrawModes 
			{
				DM_PEAKS,						//< draw data as peak
				DM_CONNECTEDLINES		//< draw as connected lines
			};

			/// Returns the draw mode of the current layer
			DrawModes getDrawMode() const;
	
			/// Sets draw mode of the current layer
			void setDrawMode(DrawModes mode);
			
			/**
				@brief Sets the overall data range
				
				This method sets the overall data range. Normally, this is computed by the protected
				function recalculateRanges_(), but when there are two canvasses which should behave
				synchronously, the united range of both canvasses has to be computed from outside
				this class (where both canvasses and their ranges are known). Currently, this
				is done in Spectrum1DWidget.
			*/
			void setOverallDataRange(const DRange<3>& overall_range);
			
			/// Calls recalculateRanges_()
			void recalculateRanges();
			
			/// Sets whether we are currently in mirror mode or not (two 1d canvasses on one widget)
			void setMirrorMode(bool mode);
			
			/// Returns whether we are currently in mirror mode or not (two 1d canvasses on one widget)
			bool inMirrorMode();
			
			// Docu in base class
			virtual void showCurrentLayerPreferences();

			// Docu in base class
			virtual void saveCurrentLayer(bool visible);
	
		public slots:
			// Docu in base class
			void activateLayer(int layer_index);
			// Docu in base class
			void removeLayer(int layer_index);
			
			/**
				@brief Sets the visible area.
				
				Sets the visible area to a new value. Note that it does not emit visibleAreaChanged()
				@param range the new visible area
			*/
			void setVisibleArea(DRange<2> range); //Do not change this to AreaType the signal needs QT needs the exact type...
			// Docu in base class
			virtual void horizontalScrollBarChange(int value);
			
			/// Adds the current visible area to the zoom stack
			void zoomAdd();
		
			/// Calls zoomForward_()
			void zoomForward();
		
			/// Calls zoomBack_()
			void zoomBack();
		
		protected:
			// Docu in base class
			bool finishAdding_();
			
			/**
				@brief Changes visible area interval
				
				This method is for convenience only. It calls changeVisibleArea_(const AreaType&, bool, bool) .
			*/
			void changeVisibleArea_(double lo, double hi, bool repaint = true, bool add_to_stack = false);  
			
			/// Calls dataToWidget_(const PointType&, QPoint& point) but takes snap_factor_ and percentage_factor_ into account.
			void dataToWidget_(const PeakType& peak, QPoint& point);
			
			/// draws a highlighted peak; if draw_elongation is true, the elongation line is drawn (for measuring)
			void drawHighlightedPeak_(const PeakIndex& peak, QPainter& painter, bool draw_elongation = false);
								
			/**
				@brief Sets the visible area
				
				Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
				If parts of the area are outside of the data area, the new area will be adjusted.
				
				@param new_area The new visible area.
				@param repaint if repainting of the widget should ne triggered
				@param add_to_stack If the new area is to add to the zoom_stack_
			*/
			virtual void changeVisibleArea_(const AreaType& new_area, bool repaint = true, bool add_to_stack = false);
			// Docu in base class
			virtual void currentLayerParamtersChanged_();
			// Docu in base class
			virtual void recalculateSnapFactor_();
			// Docu in base class
			virtual void updateScrollbars_();
			
			/// Draw modes (for each spectrum)
			std::vector<DrawModes> draw_modes_; 
			/// Iterator on peak next to mouse position
			PeakIndex selected_peak_;
			/// start peak of measuring mode
      PeakIndex measurement_start_;
      /// start point of "ruler" for measure mode
      QPoint measurement_start_point_;
      /// The annotation manager
      Annotations1DManager annotation_manager_;
      /// Indicates whether this canvas is currently in mirror mode with another 1d canvas
      bool in_mirror_mode_;
      
			/// Find peak next to the given position
			PeakIndex findPeakAtPosition_(QPoint);
	
	    /** @name Reimplemented QT events */
	    //@{
			void paintEvent(QPaintEvent* e);
			void mousePressEvent(QMouseEvent* e);
			void mouseReleaseEvent(QMouseEvent* e);
			void mouseMoveEvent(QMouseEvent* e);
			void keyPressEvent(QKeyEvent* e);

			void contextMenuEvent(QContextMenuEvent* e);
	    //@}
			
			//docu in base class
			virtual void updateLayer_(UInt i);
			//docu in base class
			virtual void translateLeft_();
			//docu in base class
			virtual void translateRight_();

		protected slots:
			
			/**
				@brief Overwrites SpectrumCanvas::recalculateRanges_()
				
				This method overwrites SpectrumCanvas::recalculateRanges_(), in order
				to be able to prevent the canvas from recalculating its ranges. This
				is needed when two 1d canvasses are shown on a single 1d widget
				(mirror view) and the overall data range has to contain the ranges
				of both canvasses. The united overall range for both canvasses
				is set in Spectrum1DWidget.
			*/
			void recalculateRanges_(UInt mz_dim, UInt rt_dim, UInt it_dim);
	};
} // namespace OpenMS

#endif
