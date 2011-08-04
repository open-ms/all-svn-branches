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


#ifndef _IDEVAL_H
#define	_IDEVAL_H

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringBayes.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>
//#include <Common/CommonTypes.h>
#include <iostream>
#include <gsl/gsl_fit.h>

/**This class shall contribute several functions to:
 *    -Perform a second stage scoring on the candidate Sequences obtained in the de novo sequencing stage
 *    -evaluate the quality of the derived sequences
 */

namespace OpenMS
{

class IdEval
{

public:

    //forward declaration
    struct IdAnnot;

protected:

    typedef std::vector<UInt>UIntVec;
    typedef std::vector<Int>IntVec;
    typedef std::vector<String>StringVec;
    typedef std::pair<DoubleReal, IdAnnot> ScoreAnnotPair;
    typedef std::pair<UInt, UInt> IndexPair;


    //compute all possible combinations of edge annotations for a single ID result (e set of edge lengths)
    static UInt getAllCombinations_(const std::vector<StringVec> &input, std::vector<IdAnnot> &output);

public:

    ///returns for an ordered list of prefix masses one possible annotation
    static UInt generateSampleAnnotation(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UInt> &prefix_masses, IdAnnot & annotation);

    ///returns a vector of all possible annotations (all combinations and permutations) for a given single ID result (vector of edge masses)
    static UInt generateAllAnnotations(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision,  const std::vector<UInt> &prefix_masses,  std::vector<IdAnnot> &all_annotations);

    ///returns the absolute number and percentage of correctly identified amino acids for some list of prefix
    ///masses and a given annotation
    static UInt predictionQuality(const IdAnnot & cand_annot, const IdAnnot & true_annot, UInt &correct_AA, DoubleReal & correct_percentage);

    ///returns the absolute number and percentage of correctly identified amino acids for some list of prefix
    ///masses and a given annotation
    static UInt predictionQuality(const std::vector<UInt> &prefix_masses, const DoubleReal &precision, const AASequence & true_seq, UInt &correct_AA, DoubleReal & correct_percentage);

    ///calculate the best possible quality given a set of edge weights and the true sequence
    static UInt possiblePredictionQuality(const IdSetup::AASeqVecMap &annot_map, const std::vector<UInt> &edge_masses, const DoubleReal &precision, const AASequence & true_seq);

    ///rescore gap edge using the given annotation, beginning at start_mz, ending at end_mz of the input spectrum spec
    ///returns the score of this annotation for the edge
    static DoubleReal rescoreGap(const AASequence & annot, const DoubleReal start_mz, const DoubleReal end_mz, const PeakSpectrum &spec, BayesScoring &scrfnc, DoubleReal delta);

    ///perform a brute force rescoring by generating the theoretical spectrum for every possible combination and permutation and perform spectrum similarity scoring
    static void brute_rescoring(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UIntVec> &result_masses, std::multimap<DoubleReal, IdAnnot> &ranked_candidates, const PeakSpectrum & orig_spec);

    ///perform a brute force rescoring by generating the theoretical spectrum for every possible combination and permutation and perform spectrum similarity scoring
    static void smart_rescoring(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UIntVec> &result_masses, std::multimap<DoubleReal, IdAnnot> &ranked_candidates, const PeakSpectrum & orig_spec);

    ///perform a brute force rescoring by generating the theoretical spectrum for every possible combination and permutation and perform spectrum similarity scoring
    static void smart_rescoring2(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UIntVec> &result_masses, std::set<ScoreAnnotPair> &ranked_candidates, const PeakSpectrum & orig_spec, const BayesScoring &scr_func);

    ///calibrate spectrum given annotation
    template <typename T>
    static void calibrateSpectrum(MSSpectrum<T> &spec, const AASequence & annot, DoubleReal delta);

    static DoubleReal scorePSM(const PeakSpectrum & spec, const AASequence & sequence, const std::set<IdSetup::ion_type>scored_types, /*const PriorTable & priors,*/ DoubleReal delta);

    static void getCIDSpectrum(PeakSpectrum& spec, const AASequence& sequence, Size charge, DoubleReal prefix, DoubleReal suffix);




}; //end of IdEval



///provisorial structure for an annotation as returned by an ID run
///where the gaps indicate which AA's in the sequence are related to a multi edge
struct IdEval::IdAnnot
{
    friend std::ostream& operator<< (std::ostream &out, IdEval::IdAnnot &annot);

    AASequence sequence;
    std::vector<Int> chunks;
    DoubleReal c_gap;
    DoubleReal n_gap;

    //default constructor (does nothing)
    IdAnnot():sequence(AASequence()),
              chunks(0),
              c_gap(0.0),
              n_gap(0.0)
              {}

    //constructor for call with both arguments
    IdAnnot(AASequence aa_seq, DoubleReal c_gap, DoubleReal n_gap, IntVec chunks_in = IntVec()):
            sequence(aa_seq),
            c_gap(c_gap),
            n_gap(n_gap),
            chunks(chunks_in)
            {}

    //constructor for call with sequence only
    IdAnnot(AASequence aa_seq):
            sequence(aa_seq),
            c_gap(0.0),
            n_gap(0.0)
    {
        //assume no multiple edges if no chunks are input
        chunks.assign(sequence.size(),1);
    }

    //copy contructor
    IdAnnot(const IdAnnot &rhs):
      sequence(rhs.sequence),
      chunks(rhs.chunks),
      c_gap(rhs.c_gap),
      n_gap(rhs.n_gap)
    {
    }

    //assignment operator
    IdAnnot& operator=(const IdAnnot &rhs)
    {
      if(this!=&rhs)
      {
        sequence = rhs.sequence;
        chunks = rhs.chunks;
        c_gap = rhs.c_gap;
        n_gap = rhs.n_gap;
      }
      return *this;
    }

    //comparison operator
    bool operator<(const IdAnnot &rhs) const
    {
      if(n_gap!=rhs.n_gap)
        return(n_gap<rhs.n_gap);

      else if(sequence!=rhs.sequence)
        return (sequence<rhs.sequence);

      else if(c_gap!=rhs.c_gap)
        return(c_gap<rhs.c_gap);
      else
        return false;
    }

    //destructor
    ~IdAnnot()
    {
    }

    //print the annotation in a human readable format with braces e.g AY[TY]LG
    void print()
    {
        std::vector<Int>::iterator chunk_it=  chunks.begin();

        AASequence::Iterator aa_seq_it = sequence.begin();
        AASequence::Iterator right_border = sequence.begin();

        if(n_gap!=0.0)
          std::cout<<"["<<n_gap<<"]";

        for(;chunk_it!=chunks.end(); ++chunk_it)
        {
            if(*chunk_it==-1)
            {
              //std::cout<<aa_seq_it->getOneLetterCode();
              std::cout<<"["<<aa_seq_it->getName()<<"]";
              ++aa_seq_it;
              ++right_border;
            }
            else if(*chunk_it==1)
            {
              ++right_border;
              std::cout<<aa_seq_it->getOneLetterCode();
              ++aa_seq_it;
            }

            else
            {
              right_border= right_border+ (*chunk_it);
              std::cout<<"[";
              while(aa_seq_it!=right_border && aa_seq_it!=sequence.end())
              {
                  std::cout<<aa_seq_it->getOneLetterCode();
                  ++aa_seq_it;
              }
              std::cout<<"]";
            }
        }
        if(c_gap!=0.0)
          std::cout<<"["<<c_gap<<"]";
        std::cout<<std::endl;
    }
}; //end of IdAnnot

    //print the annotation in a human readable format with braces e.g AY[TY]LG
    inline std::ostream& operator<<(std::ostream& os,  IdEval::IdAnnot& annot)
    {
        std::vector<Int>::const_iterator chunk_it=  annot.chunks.begin();

        AASequence::ConstIterator aa_seq_it = annot.sequence.begin();
        AASequence::ConstIterator right_border = annot.sequence.begin();


        if(annot.n_gap!=0.0)
          os<<"["<<annot.n_gap<<"]";

        for(;chunk_it!=annot.chunks.end(); ++chunk_it)
        {
            if(*chunk_it==-1)
            {
              //std::cout<<aa_seq_it->getOneLetterCode();
              os<<"["<<aa_seq_it->getName()<<"]";
              ++aa_seq_it;
              ++right_border;
            }
            else if(*chunk_it==1)
            {
              ++right_border;
              os<<aa_seq_it->getOneLetterCode();
              ++aa_seq_it;
            }

            else
            {
              right_border= right_border+ (*chunk_it);
              os<<"[";
              while(aa_seq_it!=right_border && aa_seq_it!=annot.sequence.end())
              {
                  os<<aa_seq_it->getOneLetterCode();
                  ++aa_seq_it;
              }
              os<<"]";
            }
        }
        if(annot.c_gap!=0.0)
          os<<"["<<annot.c_gap<<"]";
        os<<std::endl;

        return os;
    }

    template<typename T>
    void IdEval::calibrateSpectrum(MSSpectrum<T> &spec, const AASequence & annot, DoubleReal delta)
    {
      //generate the mz values of b- and -yions
      //these are used as reference points for calibration
      std::vector<DoubleReal>theoretical_mz,  observed_mz;
      for(Size i=1; i<annot.size(); ++i)
      {
        DoubleReal t_mass=annot.getPrefix(i).getMonoWeight(Residue::BIon,1);
        DoubleReal o_mass= spec[spec.findNearest(t_mass)].getMZ();
        if(fabs(t_mass-o_mass) < delta)
        {
          theoretical_mz.push_back(t_mass);
          observed_mz.push_back(o_mass);
        }

        t_mass=annot.getSuffix(i).getMonoWeight(Residue::YIon,1);
        o_mass= spec[spec.findNearest(t_mass)].getMZ();
        if(fabs(t_mass-o_mass) < delta)
        {
          theoretical_mz.push_back(t_mass);
          observed_mz.push_back(o_mass);
        }
      }

      DoubleReal slope, intercept, cov00, cov01, cov11, sumsq;
      gsl_fit_linear (&observed_mz[0], 1, &theoretical_mz[0],1, theoretical_mz.size(), &intercept, &slope, &cov00,  &cov01, &cov11, &sumsq);

      std::cout<<"slope: "<<slope<<"  intercept  "<<intercept<<std::endl;

      //now use the linear function to calibrate the spectrum
      for(Size peak=0; peak<spec.size(); ++peak)
      {
        DoubleReal new_mz= spec[peak].getMZ()*slope + intercept;
        spec[peak].setMZ(new_mz);
      }

      //DEBUG
      DoubleReal err_non_cal;
      DoubleReal err_cal;
      for(Size i=0; i<theoretical_mz.size(); ++i)
      {
        err_non_cal+=fabs(theoretical_mz[i]-observed_mz[i]);
        err_cal+=fabs(theoretical_mz[i]- (observed_mz[i]*slope+intercept));
      }
    }



}//namespace



#endif	/* _IDEVAL_H */

