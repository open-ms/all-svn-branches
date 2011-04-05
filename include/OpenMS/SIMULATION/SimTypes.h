// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_SIMTYPES_H
#define OPENMS_SIMULATION_SIMTYPES_H

#include <vector>
#include <utility>
#include <map>
#include <utility>

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace OpenMS 
{
  /// Coordinate type in mz and rt dimension
  typedef Peak2D::CoordinateType SimCoordinateType;
  
	/// Abundance of proteins/peptides
	typedef Peak2D::IntensityType SimIntensityType;
	
  /// Charge of a peptide
  typedef Feature::ChargeType SimChargeType;
  
  /// Raw data point
	typedef Peak1D SimPointType;
	
  /// Container for FASTAEntry & abundance information
  typedef std::vector< std::pair<FASTAFile::FASTAEntry, MetaInfoInterface> > SampleProteins;

  /// Container for multiple channels of SampleProteins
  typedef std::vector< SampleProteins > SampleChannels;

	/// Sim FeatureMap
	typedef FeatureMap<> FeatureMapSim;

  /// Sim FeatureMap Vector
  typedef std::vector<FeatureMapSim> FeatureMapSimVector;

  /// Sim MSExperiment type
  typedef MSExperiment< SimPointType > MSSimExperiment;

  /**
    @brief Wrapper class for random number generators used by the simulation classes

    The random numbers are separated two sources of randomness:

    <ul>
      <li><em>technical random numbers</em> which should represent technical
          sources of variability like instrument noise and </li>
      <li><em>biological random numbers</em> which should represent biological
          sources of variability (e.g. between two samples of the same composition)</li>
    </ul>

    @ingroup Simulation
  */
  struct SimRandomNumberGenerator
  {
    /// GSL random number generator for biological variability
    gsl_rng* biological_rng;
    /// GSL random number generator for technical variability
    gsl_rng* technical_rng;

    /// Default constructor
    SimRandomNumberGenerator()
      : biological_rng(NULL),
      technical_rng(NULL)
    {
    }

    /** @name Constructors and Destructors
      */
    //@{
    /// Copy constructor
    SimRandomNumberGenerator(const SimRandomNumberGenerator& other)
      : biological_rng(other.biological_rng),
      technical_rng(other.technical_rng)
    {
    }

    /// Destructor
    ~SimRandomNumberGenerator()
    {
      if(biological_rng != 0)
      {
        gsl_rng_free( biological_rng );
      }

      if(technical_rng != 0)
      {
        gsl_rng_free( technical_rng );
      }
    }
    //@}

    /// Assignment operator
    SimRandomNumberGenerator& operator = (const SimRandomNumberGenerator& source)
    {
      this->biological_rng = source.biological_rng;
      this->technical_rng = source.technical_rng;

      return *this;
    }
  };

}

#endif
