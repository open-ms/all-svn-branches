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
// '$'Maintainer: Sandro Andreotti '$'
// --------------------------------------------------------------------------


#ifndef RANKSCORING_H_
#define RANKSCORING_H_

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>


namespace OpenMS
{
  class RankScoringFunction : public DefaultParamHandler
  {
  private:

    typedef std::vector<DoubleReal>DVector;
    typedef std::vector<DVector> DMatrix;
    typedef DMatrix CondProbTable;
    typedef std::vector<UInt> CountTable;
    typedef std::vector<UInt>UIntVec;

  public:
    /// Default constructor
    RankScoringFunction();
    /// Custom constructor
    RankScoringFunction(String parameter_file);
    ///Copy contructor
    RankScoringFunction(RankScoringFunction & in);
    //assigment operator
    RankScoringFunction& operator=(const RankScoringFunction& in);
    ///Learn the rank scores (the input Map has to be annotated with the true sequences)
    void generateRankProbs(const PeakMap &specs);
    ///load a trained RankScoringFunction
    void loadModel();
    ///write the model to file file_name
    void writeModel();
    ///return the rank scores for being b or y ion for each peak of the spectrum
    DoubleReal getRankScores(DVector &b_scores, DVector &y_scores, const PeakSpectrum &spec);


  private:
    void updateMembers_();

    Size number_of_sectors_;
    DoubleReal delta_;
    Size max_rank_;
    //DRealMatrix rank_prob_;
    DMatrix log_odds_b_;
    DMatrix log_odds_y_;


  };
}
#endif /* RANKSCORING_H_ */
