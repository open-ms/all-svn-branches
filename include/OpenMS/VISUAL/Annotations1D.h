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

#ifndef OPENMS_VISUAL_ANNOTATIONS1D_H
#define OPENMS_VISUAL_ANNOTATIONS1D_H

#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/VISUAL/Annotation1DEmptyItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

namespace OpenMS
{
	/** blablabla
	*/
	class Annotations1D
	{
		
		public:
			
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;
			/// Container type
			typedef std::vector<Annotation1DItem*> ContainerType;
			/// Iterator
			typedef ContainerType::iterator Ann1DIterator;
			/// const iterator
			typedef ContainerType::const_iterator Ann1DConstIterator;
			/// Constructor
			Annotations1D();
			
			/// Provided for convenience; calls getAnnotationAt(p.getX(), p.getY())
			Annotation1DItem* getAnnotationAt(PointType& p);
			/// Returns the annotation at (@p x, @p y) on the canvas, or an Annotation1DEmptyItem if not existent
			Annotation1DItem* getAnnotationAt(CoordinateType x, CoordinateType y);
			/// Returns the number of annotation items
			UInt size();
			/// Returns if the annotation_items_ vector is empty
			bool isEmpty();
			/// Adds a new annotation @p item
			void add(Annotation1DItem* item);
			/// Removes the annotation @p item . Returns true, if item was contained and removed, else false.
			bool remove(Annotation1DItem* item);
			/// Clears the annotation_items_ vector
			void clear();
			
					
		protected:
			
			/// Vector containing all the Annotation1DItems
			ContainerType annotation_items_;
			
			/// Empty annotation
			Annotation1DEmptyItem* empty_item_;

	};
} // namespace OpenMS

#endif
