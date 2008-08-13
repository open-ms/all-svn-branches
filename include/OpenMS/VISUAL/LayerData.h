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

#ifndef OPENMS_VISUAL_LAYERDATA_H
#define OPENMS_VISUAL_LAYERDATA_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/VISUAL/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Annotation1DPeakItem.h>


#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>

#include <list>

namespace OpenMS 
{
	/**
		@brief Struct that stores the data for one layer
		
		@ingroup SpectrumWidgets
	*/
  struct LayerData
	{
		/**	@name Type definitions */
		//@{
		///Dataset types
		enum DataType
		{
			DT_PEAK,		      ///< Peak/Raw data
			DT_FEATURE,	      ///< Feature data
			DT_UNKNOWN			  ///< Undefined data type indicating an error
		};

		///Flags that determine which information is shown.
		enum Flags
		{
			F_HULL,       ///< Features: Overall convex hull
			F_HULLS,      ///< Features: Convex hulls of single mass traces 
			F_NUMBERS,    ///< Features: Number
			P_PRECURSORS, ///< Peaks: Mark precursor peaks of MS/MS scans
			P_PROJECTIONS ///< Peaks: Show projections
		};
				
		/// Main data type (experiment)
		typedef MSExperiment<Peak1D> ExperimentType;
		/// Main data type (features)
		typedef FeatureMap<> FeatureMapType;
		/// Container type for the 1D annotations
		typedef std::list<Annotation1DItem*> Annotations1DContainerType;
		/// Iterator for the 1D annotations
		typedef Annotations1DContainerType::iterator Ann1DIterator;
		/// Const iterator for the 1D annotations
		typedef Annotations1DContainerType::const_iterator Ann1DConstIterator;
		//@}
		
		/// Default constructor
		LayerData()
			: visible(true),
				type(DT_UNKNOWN),
				name(),
				filename(),
				peaks(),
				features(),
				f1(false),
				f2(false),
				f3(false),
				param(),
				gradient(),
				filters(),
				annotations_1d_()
		{
		}
		
		/// Copy constructor
		LayerData(const LayerData& rhs)
		{
			visible = rhs.visible;
			type = rhs.type;
			name = rhs.name;
			filename = rhs.filename;
			peaks = rhs.peaks;
			features = rhs.features;
			f1 = rhs.f1;
			f2 = rhs.f2;
			f3 = rhs.f3;
			param = rhs.param;
			gradient = rhs.gradient;
			filters = rhs.filters;
			//copy annotations
			Annotation1DItem* new_item = 0;
			for (Ann1DConstIterator it = rhs.annotations_1d_.begin(); it != rhs.annotations_1d_.end(); ++it)
			{
				Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(*it);
				if (distance_item)
				{
					new_item = new Annotation1DDistanceItem(*distance_item);
					annotations_1d_.push_back(new_item);
					continue;
				}
				Annotation1DTextItem* text_item = dynamic_cast<Annotation1DTextItem*>(*it);
				if (text_item)
				{
					new_item = new Annotation1DTextItem(*text_item);
					annotations_1d_.push_back(new_item);
					continue;
				}
				Annotation1DPeakItem* peak_item = dynamic_cast<Annotation1DPeakItem*>(*it);
				if (peak_item)
				{
					new_item = new Annotation1DPeakItem(*peak_item);
					annotations_1d_.push_back(new_item);
					continue;
				}
			}
		}
		
		/// Destructor
		virtual ~LayerData()
		{
			for (Ann1DIterator it = annotations_1d_.begin(); it != annotations_1d_.end(); ++it)
			{
				delete *it;
			}
		}
		
		/// if this layer is visible
		bool visible;
		/// data type (peak of feature data)
		DataType type;
		/// layer name
		String name;
		/// file name of the file the data comes from (if available)
		String filename;
		/// peak data
		ExperimentType peaks;
		/// feature data
		FeatureMapType features;
		
		/// Flag one (Features: convex hulls, Peak: precursors)
		bool f1;
		/// Flag two (Features: numbers, Peak: projections)
		bool f2;
		/// Flag tree (Features: convex hull, Peak: -)
		bool f3;
		
		///Layer parameters
		Param param;
		
		///Gradient for 2D and 3D views
		MultiGradient gradient;
		
		///Filters to apply before painting
		DataFilters filters;
				
		///Annotations for the 1D view
		mutable Annotations1DContainerType annotations_1d_;
		
	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
