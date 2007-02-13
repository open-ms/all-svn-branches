// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_VISUAL_COLORSELECTOR_H
#define OPENMS_VISUAL_COLORSELECTOR_H

//QT
#include <qwidget.h>
//Added by qt3to4:
#include <QPaintEvent>
#include <QMouseEvent>

namespace OpenMS 
{
	
	/**
		@brief A widget for selecting a color.
		
		It represents a color (displayed as background color) and allows changing the color.
		
		\image html ColorSelector.png
		
		The above example image shows four ColorSelector instances on the right side.
		
		@ingroup Visual
	*/
	class ColorSelector: public QWidget
	{
		Q_OBJECT
		
		public:
			/// Constructor
			ColorSelector( QWidget * parent = 0, const char * name = 0 );
			/// Destructor
			~ColorSelector();
			
			/// Returns the selected color
			const QColor& getColor();
			/// Sets the selected color
			void setColor(const QColor&);
		
		protected:
			void paintEvent(QPaintEvent * e);
			void mousePressEvent(QMouseEvent* e);
			QColor color_;
	};

}
#endif // OPENMS_VISUAL_COLORSELECTOR_H

