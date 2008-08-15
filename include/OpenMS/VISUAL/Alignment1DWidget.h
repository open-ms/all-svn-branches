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

#ifndef OPENMS_VISUAL_ALIGNMENT1DWIDGET_H
#define OPENMS_VISUAL_ALIGNMENT1DWIDGET_H

#include <QtGui/QWidget>

#include <vector>
#include <utility>

#include <OpenMS/CONCEPT/Types.h>

class QPaintEvent;

namespace OpenMS
{
	class Spectrum1DWidget;
	class SpectrumCanvas;
	
	/**
		@brief Widget visualizing an alignment of two spectra in mirror view
		
		This widget visualizes an alignment of two spectra in mirror view. It is
		located between the two canvasses. If no alignment has been performed yet, it
		has a fixed height of 1 pixel and draws only an x-axis in order to separate
		both spectra. If an alignment is performed, it changes its height to 10 pixels,
		and draws lines connecting the aligned peaks of the two spectra.
	*/
	class Alignment1DWidget
		: public QWidget
	{
		Q_OBJECT
		
		public:
			
			/// Default constructor
			Alignment1DWidget(Spectrum1DWidget* spectrum_widget, QWidget* parent = 0);
			/// Destructor
			virtual ~Alignment1DWidget();
		
			/// Sets alignment_lines_ to @p alignment_lines
			void setAlignmentLines(const std::vector<std::pair<DoubleReal,DoubleReal > >& alignment_lines);
			
			/// Clears alignment_lines_
			void clearAlignmentLines();
				
			/// Returns whether an alignment has been set
			bool alignmentIsSet();
	
		public slots:
		
		protected:
			
      /// Vector containing all the peak-connecting lines of the spectrum alignment (stored as m/z values)
      std::vector<std::pair<DoubleReal, DoubleReal > > alignment_lines_;
 			/// The enclosing spectrum widget
 			Spectrum1DWidget* spectrum_widget_;
			/// Stores whether an alignment has been set
			bool alignment_is_set_;
	    /** @name Reimplemented QT events */
	    //@{
			void paintEvent(QPaintEvent* e);
	    //@}
	    
	    /// Returns the SpectrumCanvas object of the upper canvas
	    SpectrumCanvas* upperCanvas_();
	    
	    /// Returns the SpectrumCanvas object of the lower canvas
			SpectrumCanvas* lowerCanvas_();
			
		protected slots:
			
	};
} // namespace OpenMS

#endif
