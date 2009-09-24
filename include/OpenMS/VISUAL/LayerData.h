// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_LAYERDATA_H
#define OPENMS_VISUAL_LAYERDATA_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/Annotations1DContainer.h>
#include <OpenMS/VISUAL/MappingThread.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>

#include <vector>

namespace OpenMS 
{
	/**
		@brief Struct that stores the data for one layer
		
		@ingroup SpectrumWidgets
	*/

  class LayerData
	{
	  public:
		  /**	@name Type definitions */
		  //@{
		  ///Dataset types
		  enum DataType
		  {
			  DT_PEAK,		      ///< Spectrum profile or centroided data
			  DT_FEATURE,	      ///< Feature data
			  DT_CONSENSUS,     ///< Consensus feature data
			  DT_CHROMATOGRAM,  ///< Chromatogram data
			  DT_UNKNOWN			  ///< Undefined data type indicating an error
		  };

		  ///Flags that determine which information is shown.
		  enum Flags
		  {
			  F_HULL,       ///< Features: Overall convex hull
			  F_HULLS,      ///< Features: Convex hulls of single mass traces 
			  F_UNASSIGNED,///< Features: Unassigned peptide hits
			  P_PRECURSORS, ///< Peaks: Mark precursor peaks of MS/MS scans
			  P_PROJECTIONS,///< Peaks: Show projections
			  C_ELEMENTS    ///< Consensus features: Show elements
		  };

		  ///Enumerate all avaiable paint styles
		  enum PrimitiveModes 
		  {
			  PM_POINTS,
			  PM_LINES,
			  PM_LINESTRIP,
			  PM_TRIANGLES,
			  PM_TRIANGLESTRIP,
			  PM_QUADS,
			  PM_QUADSTRIP,
			  PM_POLYGON
		  };
		
		  ///Label used in visualization
		  enum LabelType
		  {
			  L_NONE,							///< No label is displayed
			  L_INDEX,						///< The element number is used
			  L_META_LABEL,				///< The 'label' meta information is used
			  L_ID,								///< The best peptide hit of the first identification run is used
			  L_ID_ALL,						///< All peptide hits of the first identification run are used
			  SIZE_OF_LABEL_TYPE
		  };
		  
		  //Label names
		  static const std::string NamesOfLabelType[SIZE_OF_LABEL_TYPE];
		
		  /// Main data type (experiment)
		  typedef MSExperiment<> ExperimentType;
		  /// Main data type (features)
		  typedef FeatureMap<> FeatureMapType;
		  /// Main data type (consensus features)
		  typedef ConsensusMap ConsensusMapType;
		  //@}
		
		  /// Default constructor.
		  LayerData();

		  /// Copy constructor.
		  LayerData(const LayerData& layer);
		  
		  /// Destructor.
		  ~LayerData();

			/// Assignment operator. 
			LayerData& operator= (const LayerData& layer);

		public:			
		  /// Returns a const reference to the current spectrum (1d view)
		  const ExperimentType::SpectrumType& getCurrentSpectrum() const;

		  /// Returns a const reference to the annotations of the current spectrum (1d view)
		  const Annotations1DContainer& getCurrentAnnotations() const;

		  /// Returns a mutable reference to the current spectrum (1d view)
		  ExperimentType::SpectrumType& getCurrentSpectrum();

		  /// Returns a mutable reference to the annotations of the current spectrum (1d view)
		  Annotations1DContainer& getCurrentAnnotations();
		  
		  /// Returns mutable reference to calcul thread
		  MappingThread* getMappingThread();
		  
		  /// Set mapping mode
		  void setMappingMode(const MappingThread::MappingModes mode);
		  
		  /// Set primitive mode
		  void setPrimitiveMode(const PrimitiveModes mode);
		  
		  /// Get mapping mode
		  MappingThread::MappingModes getMappingMode() const;
		  
		  /// Set primitive mode
		  PrimitiveModes getPrimitiveMode() const;
		  
		  /// Reset mapping data
		  void resetMapping();
		  
		  /// if this layer is visible
		  bool visible;
		  /// if this layer is flipped (1d mirror view)
		  bool flipped;
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
		  /// consensus feature data
		  ConsensusMapType consensus;
		  /// Index of the current spectrum (1d view)
		  Size current_spectrum;
		
		  /// Flag one (Features: convex hulls, Peak: precursors, Consensus: elements)
		  bool f1;
		  /// Flag two (Features: numbers, Peak: projections, Consensus: -)
		  bool f2;
		  /// Flag tree (Features: convex hull, Peak: -, Consensus: -)
		  bool f3;
		
		  ///Layer parameters
		  Param param;

		  ///Gradient for 2D and 3D views
		  MultiGradient gradient;
		
		  ///Filters to apply before painting
		  DataFilters filters;
				
		  ///Annotations of all spectra of the experiment (1D view)
		  std::vector<Annotations1DContainer> annotations_1d;
		
		  ///Flag that indicates if the layer data can be modified (so far used for features only)
		  bool modifiable;
		  ///Flag that indicates that the layer data was modified since loading it
		  bool modified;
		
		  ///Label type
		  LabelType label;
		
		private:
		  // mapping thread
		  MappingThread* mapping_thread_;
		  
			// mapping mode
			MappingThread::MappingModes mapping_mode_;

		  // primitive mode
		  PrimitiveModes primitive_mode_;
	};

	///Print the contents to a stream.
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
