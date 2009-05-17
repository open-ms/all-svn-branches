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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_IONIZATIONSIMULATION_H
#define OPENMS_SIMULATION_IONIZATIONSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// STL includes
#include <set>

namespace OpenMS {

  /**
   @brief Simulates Protein ionization
   
   Supports ESI and MALDI. The abundance values are distributed among
   the charge states based on a binomial distribution for the ESI and 
   based on discrete distribution for MALDI.
	 In ESI mode, this class also supports different adduct types in addition to H+ 
	 (e.g. NH4+, K+) which can be specified by the user and influence
	 the mass and induce more charge variation.
   
   @htmlinclude OpenMS_IonizationSimulation.parameters
   
   @ingroup Simulation
  */
  class OPENMS_DLLAPI IonizationSimulation
    : public DefaultParamHandler
  {

  public: 
    /// possible ionization methods
    typedef enum 
		{
      MALDI,
      ESI
    } IonizationType;
    
    /** @name Constructors and Destructors
      */
    //@{
    /// 
    IonizationSimulation(const gsl_rng* );
    
    /// Copy constructor
    IonizationSimulation(const IonizationSimulation& source);

    /// Destructor
    virtual ~IonizationSimulation();
    //@}

    /// Assignment operator
    IonizationSimulation& operator = (const IonizationSimulation& source);
    
    /**
     @brief Ionize all peptide features inside the Feature-Map
     
     Depending on the parameters the passed peptide features are ionized by MALDI
     or by ESI.

     @param features FeatureMap which will be ionized
		 @param charge_consensus ConsensusMap which groups childs(=charge variants) of input-features
     */
    void ionize(FeatureMapSim & features, ConsensusMap & charge_consensus);

  private:  
    /// Default constructor
    IonizationSimulation();
    
		/// isonize using ESI
    void ionizeEsi_(FeatureMapSim &, ConsensusMap & charge_consensus);
    
		/// isonize using MALDI
    void ionizeMaldi_(FeatureMapSim &, ConsensusMap & charge_consensus);

    /// set defaults
    void setDefaultParams_();
    
    /// Synchronize members with param class
		void updateMembers_();        
    
    /**
     @brief counts all basic residues inside the amino acid sequence to give an upper bound on the maximal charge during ESI ionization
    */
    UInt countIonizedResidues_(const AASequence & ) const;
    
		
		// Members //

		/// ESI or MALDI ionization
    IonizationType ionization_type_;

		/*
     @brief List of residues that are counted as basic during execution of countBasicResidues_
    */
    std::set<String> basic_residues_;
    
    /**
     @brief Probability for the binomial distribution of ESI charge states
     */
    DoubleReal esi_probability_;

		/**
		 @brief Discrete distribution of impure charge adducts like Na+, K+, Ca++ etc besides the usual H+
		*/
		// important: leave that as vector<double> because gsl expects 'double' and not 'DoubleReal' (which might be something different)
		std::vector<double> esi_impurity_probabilities_;

   
		/**
		 @brief Corresponding table to @p esi_impurity_probabilities_ holding the actual element and its charge
		*/
		Adduct::AdductsType esi_adducts_;

		/**
		 @brief Maximal charge that any impure adduct from parameter list has
		*/
		Size max_adduct_charge_;

		/**
		 @brief Preprocessed table of discrete distribution (MALDI charges)
		*/
		// important: leave that as vector<double> because gsl expects 'double' and not 'DoubleReal' (which might be something different)
		std::vector<double> maldi_probabilities_;

  protected:
		/// Random number generator
		const gsl_rng* rnd_gen_;
  };

}

#endif
