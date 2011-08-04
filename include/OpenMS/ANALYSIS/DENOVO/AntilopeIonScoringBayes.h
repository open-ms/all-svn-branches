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
// $Maintainer: Sandro Andreotti $
// --------------------------------------------------------------------------

#ifndef BAYESSCORING_H
#define BAYESSCORING_H

#include <OpenMS/ANALYSIS/DENOVO/DeNovoIonScoring.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeBayesNetwork.h>

namespace OpenMS
{

  class BayesScoring
  {    
  private:

    typedef std::vector<DoubleReal>DRealVector;
    typedef std::vector<DRealVector> DRealMatrix;
    typedef DRealMatrix CondProbTable;
    typedef std::vector<UInt> CountTable;
    typedef std::vector<UInt>UIntVec;

//--------------------------------------------------------------------------------------------------------------
//----------------------------------------ScoringFunctionMB-----------------------------------------------------
//--------------------------------------------------------------------------------------------------------------

  protected:

    typedef IdSetup::ion_type ion_type;

    /// precision of algorithms  default: 0.1 DA
    DoubleReal precision_;

    ///allowed m/z error for each peak  default 0.5 DA
    DoubleReal delta_;

    /// the number of peak intensity levels
    UInt number_intensity_levels_;

    /// the intensity limits (for spectrum normalization)
    DRealVector intensity_limits_;

    /// the number of sectors each spectrum is split into
   Size number_of_sectors_;

   /// the set of selected ion_types
   std::vector<std::set<IdSetup::ion_type> >selected_types_;

   ///Bayesian Network models - One for each sector
   std::vector<BayesianNetwork>networks_;

   std::vector<int>node_map_;

   std::vector<int>node_map_rev_;   

  public:

    ///default constructor
    BayesScoring();

    ///constructor with single argument - filename of the learned parameters
    BayesScoring(String parameter_file);

    ///copy constructor
    BayesScoring(const BayesScoring &);

    ///destructor
    virtual ~BayesScoring()
    {
    };    

    ///assignment operator
    BayesScoring& operator=(const BayesScoring& in);

    ///load the statistic parameter file
    UInt loadFile(String parameter_file);

    ///return the log score of a given prefix mass given the MS/MS Spectrum
    DoubleReal getLogScore(const PeakSpectrum& spectrum, DoubleReal prefix_mass, std::vector<IdSetup::interpretation> &interprets);

    ///return the log score under Null-hypothesis
    DoubleReal getPrandScore(std::vector<UInt> &intensities, DoubleReal m, const PeakSpectrum &S) const;

    DoubleReal getSupportIonCount(const PeakSpectrum& spectrum, const DoubleReal prefix_mass) const;

    //normalize peak intensities as in PepNovo paper
    UInt normalizeIntensity(PeakSpectrum &S)const;

    //normalize peak intensities as in PepNovo paper
    static UInt normalizeIntensity(PeakSpectrum &S, Size number_intens_levels);

    //set delta
    void setDelta(DoubleReal d_in){delta_=d_in;}
  };



}//namespace


#endif // BAYESSCORING_H
