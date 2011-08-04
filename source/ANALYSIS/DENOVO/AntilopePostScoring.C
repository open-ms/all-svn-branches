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

//#include <iostream>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>
#include <OpenMS/ANALYSIS/DENOVO/LagrangeProblemBase.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePostScoring.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringBayes.h>



#include <algorithm>
#include <map>
#include <sstream>
#include <numeric>

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>

namespace OpenMS
{

  //generates one possible (the shortest) annotation for the returned edge weights
  UInt IdEval::generateSampleAnnotation(const IdSetup::AASeqVecMap &mass_seq_map, const DoubleReal &precision, const std::vector<UInt> &prefix_masses, IdAnnot & annotation)
  {

    std::vector<UInt>::const_iterator p_mass_it;

    IdSetup::AASeqVecMap::const_iterator st_set_map_it;

    for (p_mass_it = prefix_masses.begin(); p_mass_it != prefix_masses.end(); ++p_mass_it)
    {
      //calculate the mass index
      UInt index_value = *p_mass_it;
      st_set_map_it = mass_seq_map.find(index_value);

      if (st_set_map_it != mass_seq_map.end())
      {
        annotation.sequence += (st_set_map_it->second)[0];
        annotation.chunks.push_back((st_set_map_it->second)[0].size());
      }
      else
      {
        std::stringstream tmp;
        tmp << "[" << *p_mass_it << "]";
        //annotation.sequence+=*p_mass_it;
        annotation.sequence += tmp.str();
        annotation.chunks.push_back(-1);

      }
    }

    return 0;
  }

  //compute all possible combinations of edge annotations for a single ID result (e set of edge lengths)
  UInt IdEval::getAllCombinations_(const std::vector<StringVec>& input, std::vector<IdAnnot>&output)
  {
    //number of Combinations
    UInt count = 0;

    //length of the input vector
    UInt size = input.size();

    //vector of iterators
    std::vector<StringVec::const_iterator> iters;

    //iterator on vector of iterators
    std::vector<StringVec::const_iterator>::iterator res_iter;

    //temporary AASequence
    String tmp_sequence;

    for (UInt i = 0; i < size; ++i)
    {
      iters.push_back(input[i].begin());
    }

    std::vector<String> temp(size);
    UInt total_combs = 1;
    for (Size i = 0; i < input.size(); ++i)
    {
      total_combs *= input[i].size();
    }

    //terminate when the last iterator has reached the end of the last vector
    while (count < total_combs)
    {
      std::cout << "in while  " << count << std::endl;
      ++count;

      UInt i = 0;
      while (i < size)
      {
        if (++iters[i] == input[i].end())
        {
          iters[i] = input[i].begin();
          ++i;
        }
        else
          break;
      }

      //output.push_back(IdAnnot());
      res_iter = iters.begin();

      tmp_sequence.clear();

      //now add the current combination to the output vector
      while (res_iter != iters.end())
      {
        tmp_sequence += **res_iter++;
        //output.back().push_back(**res_iter++);
      }

      std::cout << "in combination:  " << tmp_sequence << std::endl;
      output.push_back(IdAnnot(AASequence(tmp_sequence)));
    }
    return count;
  }

  //this function operates in three steps as described below
  UInt IdEval::generateAllAnnotations(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UInt> &result_masses,
      std::vector<IdAnnot> &all_candidates)
  {
    std::vector<std::vector<AASequence> > annotations;

    std::vector<std::vector<String> > all_permutations;
    std::vector<String> permutations;

    std::vector<IdSetup::AASeqVec>::iterator all_annot_it;
    IdSetup::AASeqVec::iterator single_annot_it;

    //first step, for each of the edge masses a vector of candidate annotations is returned and stored in annotations
    for (UInt i = 0; i < result_masses.size(); ++i)
    {
      //UInt index_value=(UInt)floor(result_masses[i]/precision +0.5);
      UInt index_value = result_masses[i];
      if (annot_map.find(index_value) != annot_map.end())
      {
        std::cerr<<"vectorSIZE: "<<annot_map.find(index_value)->second.size()<<std::endl;
        if(annot_map.find(index_value)->second[0].size()==1)
        {
          std::cerr<<"In if: "<<annot_map.find(index_value)->second[0]<<std::endl;
          annotations.push_back(IdSetup::AASeqVec(1,annot_map.find(index_value)->second[0]));
        }
        else
        {
          annotations.push_back(annot_map.find(index_value)->second);
        }
      }
      else
      {
        std::stringstream tmp;
        tmp << "[" << result_masses[i] * precision << "]";
        annotations.push_back(IdSetup::AASeqVec(1, AASequence(tmp.str())));
      }
    }
    //TODO this is now a hack to prevent resolution for last edge in tryptic ids
    //annotations.back().resize(1);

    //for each of the candidate annotations generate all possible permutations
    for (all_annot_it = annotations.begin(); all_annot_it != annotations.end(); ++all_annot_it)
    {
      for (single_annot_it = all_annot_it->begin(); single_annot_it != all_annot_it->end(); ++single_annot_it)
      {
        //compute all permutations for the given edge annotation of a single edge
        //getPermutations(*single_annot_it, permutations);
        String annot_string = single_annot_it->toString();
        if (annot_string[0] == '[')
        {
          permutations.push_back(annot_string);
        }
        else
        {
          //first the String needs to be sorted
          stable_sort(annot_string.begin(), annot_string.end());
          //UInt count=0;
          permutations.push_back(annot_string);
          //then all permutations are appended to the vector
          while (next_permutation(annot_string.begin(), annot_string.end()))
          {
            permutations.push_back(annot_string);
            //++count;
          }
        }
      }

      all_permutations.push_back(permutations);
      std::cout << "number of perms: " << permutations.size() << std::endl;

      permutations.clear();
      //}
    }

    return getAllCombinations_(all_permutations, all_candidates);
  }

  //TODO: modify this function to allow for smarter checking in case of gaps. actually they are also identified as matches but
  //Actually we call this function for a sample annotation so there are no gaps anyway.

  //calculate the quality of the prediction compared to the hopefully true annotation
  UInt IdEval::predictionQuality(const IdAnnot & cand_annot, const IdAnnot & true_annot, UInt &correct_AA, DoubleReal & correct_percentage)
  {

    const AASequence & cand_seq = cand_annot.sequence;
    const AASequence & true_seq = true_annot.sequence;

    correct_AA = 0;
    Size i = 0;
    Size j = 0;

    const DoubleReal cand_n_term_weight = cand_annot.n_gap;
    const DoubleReal true_n_term_weight = true_annot.n_gap;

    while (i < cand_seq.size() && j < true_seq.size())
    //while(cand_it!=cand_annot.sequence.end()-1 && true_it!=true_annot.sequence.end()-1)
    {
      if (fabs(cand_n_term_weight + cand_seq.getPrefix(i).getMonoWeight(Residue::Internal) - (true_n_term_weight + true_seq.getPrefix(j).getMonoWeight(Residue::Internal))) <= 2.5)
      {
        if(fabs (cand_seq[i].getMonoWeight()- true_seq[j].getMonoWeight())<0.1)//comparing residues not ok since K==Q and I==L
//        if (fabs(cand_n_term_weight + cand_seq.getPrefix(i + 1).getMonoWeight(Residue::Internal) - (true_n_term_weight + true_seq.getPrefix(j + 1).getMonoWeight(Residue::Internal))) <= 0.5)
        {
          ++correct_AA;
        }

        ++i;
        ++j;
      }

      else if (cand_n_term_weight + cand_seq.getPrefix(i).getMonoWeight(Residue::Internal) > true_n_term_weight + true_seq.getPrefix(j).getMonoWeight(Residue::Internal))
      {
        ++j;
      }

      else
      {
        ++i;
      }
    }

    correct_percentage = 100 * ((DoubleReal) correct_AA / true_seq.size());

    return 0;
  }

  //calculate the quality of the prediction compared to the hopefully true annotation
  UInt IdEval::predictionQuality(const std::vector<UInt> &edge_masses, const DoubleReal &precision, const AASequence & true_seq, UInt &correct_AA, DoubleReal & correct_percentage)
  {

    correct_AA = 0;
    UInt i = 0;
    UInt j = 0;

    DoubleReal predicted_prefix_mass = 0;//edge_masses[0]*precision;
    while (i < edge_masses.size() && j <= true_seq.size())
    //while(cand_it!=cand_annot.sequence.end()-1 && true_it!=true_annot.sequence.end()-1)
    {
      std::cout << predicted_prefix_mass << "  :!!: " << true_seq.getPrefix(j).getMonoWeight(Residue::Internal) << std::endl;
      DoubleReal predicted_prefix_mass_next = predicted_prefix_mass + edge_masses[i] * precision;
      if (fabs(predicted_prefix_mass - true_seq.getPrefix(j).getMonoWeight(Residue::Internal)) <= 0.5)
      {
        std::cout << predicted_prefix_mass_next << "  :++!!++: " << true_seq.getPrefix(j + 1).getMonoWeight(Residue::Internal) << std::endl;
        if (fabs(predicted_prefix_mass_next - true_seq.getPrefix(j + 1).getMonoWeight(Residue::Internal)) <= 0.5)
        {
          ++correct_AA;
        }
        predicted_prefix_mass = predicted_prefix_mass_next;
        ++i;
        ++j;
      }
      else if (predicted_prefix_mass > true_seq.getPrefix(j).getMonoWeight(Residue::Internal))
      {
        ++j;
      }
      else
      {
        predicted_prefix_mass = predicted_prefix_mass_next;
        ++i;
      }
    }
    correct_percentage = 100 * ((DoubleReal) correct_AA / true_seq.size());

    return 0;
  }


  //calculate the quality of the prediction compared to the hopefully true annotation
  UInt IdEval::possiblePredictionQuality(const IdSetup::AASeqVecMap &annot_map, const std::vector<UInt> &edge_masses, const DoubleReal &precision, const AASequence & true_seq)
  {
    std::cout<<"CHECK true_seq: "<<true_seq.toString()<<std::endl;
    std::vector<IdAnnot>candidates;
    std::vector<IdAnnot>::iterator cand_it;
    UInt c_gap, n_gap=0;

    //set the cgap to the complete weight
    c_gap=std::accumulate(edge_masses.begin(), edge_masses.end(), 0);

    UInt possible_correct_residues =0;

    for (UIntVec::const_iterator mass_it = edge_masses.begin(); mass_it != edge_masses.end(); ++mass_it)
    {
      c_gap-=*mass_it;

      //all combinations of permutations
      candidates.clear();
      generateAllAnnotations(annot_map, precision, UIntVec(1,*mass_it), candidates);
      UInt local_max_residues =0;

      for(cand_it=candidates.begin(); cand_it !=candidates.end(); ++cand_it)
      {
        IdAnnot tmp_annot(*cand_it);
        tmp_annot.c_gap=c_gap*precision;
        tmp_annot.n_gap=n_gap*precision;
        DoubleReal correct_per;
        UInt correct_residues;
        predictionQuality(tmp_annot, true_seq, correct_residues, correct_per);
        std::cout<<"comparing: "<<tmp_annot<<"  "<<true_seq<<std::endl;
        local_max_residues=std::max(local_max_residues, correct_residues);
      }
      possible_correct_residues+=local_max_residues;
      n_gap+=*mass_it;
    }

    return possible_correct_residues;
  }


  DoubleReal IdEval::rescoreGap(const AASequence & annot, const DoubleReal start_mz, const DoubleReal end_mz, const PeakSpectrum &spec, BayesScoring &scrfnc,
      const DoubleReal delta)
  {
    //check whether this annotation can possibly fill the gap, if not return score MINUS_INFINITY
    if (fabs(annot.getMonoWeight(Residue::Internal) - (end_mz - start_mz)) > delta)
    {
      return LagrangeProblem::MINUS_INF;
    }

    DoubleReal score = 0.0;
    std::vector<IdSetup::interpretation> interprets;

    for (Size i = 0; i < annot.size(); ++i)
    {
      score += scrfnc.getLogScore(spec, start_mz + annot.getPrefix(i).getMonoWeight(Residue::Internal), interprets);
    }
    //return score;
    return interprets.size() / annot.size();
  }

  void IdEval::brute_rescoring(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UIntVec> &result_masses,
      std::multimap<DoubleReal, IdAnnot> &ranked_candidates, const PeakSpectrum & orig_spec)
  {
    std::vector<IdAnnot> all_candidates;
    //run over all predictions from path search and generate all combinations and permutations of all ambiguities and perform similarity scoring of their theoretical spectra
    for (std::vector<UIntVec>::const_iterator cand_it = result_masses.begin(); cand_it != result_masses.end(); ++cand_it)
    {
      //all combinations of permutations are appended to all_candidates
      generateAllAnnotations(annot_map, precision, *cand_it, all_candidates);
    }

    TheoreticalSpectrumGenerator t_gen;
    Param params;
    params.setValue("add_isotopes", 1);
    params.setValue("add_losses", 1);
    t_gen.setParameters(params);

    //remove duplicate entries
    std::set<IdAnnot> annot_set;
    for (std::vector<IdAnnot>::iterator annot_it = all_candidates.begin(); annot_it != all_candidates.end(); ++annot_it)
    {
      if(annot_it->sequence.hasSuffix(AASequence("K")) || annot_it->sequence.hasSuffix(AASequence("R")))//TODO this is a hack to prevent non tryptic predictions
        annot_set.insert(*annot_it);
    }

    all_candidates.clear();
    all_candidates.assign(annot_set.begin(), annot_set.end());

    //perform similarity scoring
    for (std::vector<IdAnnot>::iterator annot_it = all_candidates.begin(); annot_it != all_candidates.end(); ++annot_it)
    {
      RichPeakSpectrum tmp;
      t_gen.getSpectrum(tmp, annot_it->sequence, 1);
      std::cout << "Sequence: " << annot_it->sequence << std::endl;
      t_gen.addPeaks(tmp, annot_it->sequence, Residue::AIon);
      t_gen.addPeaks(tmp, annot_it->sequence, Residue::BIon, 2);
      t_gen.addPeaks(tmp, annot_it->sequence, Residue::YIon, 2);

      //transform the RichPeakSpectrum into a PeakSpectrum
      PeakSpectrum tmp_peak_spec;
      for (RichPeakSpectrum::iterator spec_it = tmp.begin(); spec_it != tmp.end(); ++spec_it)
      {
        Peak1D tmp_peak;
        tmp_peak.setIntensity(spec_it->getIntensity());
        tmp_peak.setMZ(spec_it->getMZ());
        tmp_peak_spec.push_back(tmp_peak);
      }

      //DoubleReal score = ZhangSimilarityScore()(tmp_peak_spec, orig_spec);
      DoubleReal score = SpectrumCheapDPCorr()(tmp_peak_spec, orig_spec);
      //DoubleReal score = SpectrumAlignmentScore()(tmp_peak_spec, orig_spec);
      //DoubleReal score = SteinScottImproveScore()(tmp_peak_spec, orig_spec);
      //DoubleReal score = BinnedSharedPeakCount()(BinnedSpectrum(0.5, 1, tmp_peak_spec), BinnedSpectrum(0.5, 1, orig_spec));
      ranked_candidates.insert(std::pair<DoubleReal, IdAnnot>(score, *annot_it));
    }
  }

  void IdEval::smart_rescoring(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UIntVec> &result_masses,
        std::multimap<DoubleReal, IdAnnot> &ranked_candidates, const PeakSpectrum & orig_spec)
  {
    std::vector<IdAnnot>all_filtered_candidates,  all_candidates;
    std::vector<std::vector<IdAnnot> > all_edge_annots;
    std::vector<StringVec>selected_edge_annots;
    std::multimap<DoubleReal, String> scr_map;

    TheoreticalSpectrumGenerator t_gen;
    Param params;
    params.setValue("add_isotopes", 1);
    params.setValue("add_losses", 1);
    t_gen.setParameters(params);

    //run over all predictions
    for (std::vector<UIntVec>::const_iterator cand_it = result_masses.begin(); cand_it != result_masses.end(); ++cand_it)
    {
      selected_edge_annots.assign(cand_it->size(), StringVec());
      for(UIntVec::const_iterator edge_it= cand_it->begin(); edge_it!=cand_it->end(); ++edge_it)
      {
        all_candidates.clear();
        //all combinations of permutations are appended to all_candidates
        generateAllAnnotations(annot_map, precision, std::vector<UInt>(1,*edge_it), all_candidates);
        all_edge_annots.push_back(all_candidates);
      }

      //TODO  this is a hack to prevent resolution of last edge in tryptic peptides
      all_edge_annots.back().erase(all_edge_annots.back().begin(), all_edge_annots.back().begin()+2);

      for(Size i=0; i<all_edge_annots.size();++i)
      {
        for(Size j=0; j<all_edge_annots[i].size(); ++j)
        {
          String sequence;
          for(Size kk=0; kk<all_edge_annots.size(); ++kk)
          {
            if(kk!=i)
            {
              sequence+=all_edge_annots[kk][0].sequence.toString();
            }
            else
            {
              sequence+=all_edge_annots[kk][j].sequence.toString();
            }
          }

          RichPeakSpectrum tmp;
          t_gen.getSpectrum(tmp, sequence);
          t_gen.addPeaks(tmp, sequence, Residue::AIon);
          t_gen.addPeaks(tmp, sequence, Residue::BIon, 2);
          t_gen.addPeaks(tmp, sequence, Residue::YIon, 2);

          //transform the RichPeakSpectrum into a PeakSpectrum
          PeakSpectrum tmp_peak_spec;
          for (RichPeakSpectrum::iterator spec_it = tmp.begin(); spec_it != tmp.end(); ++spec_it)
          {
            Peak1D tmp_peak;
            tmp_peak.setIntensity(spec_it->getIntensity());
            tmp_peak.setMZ(spec_it->getMZ());
            tmp_peak_spec.push_back(tmp_peak);
          }

          DoubleReal score = SpectrumAlignmentScore()(tmp_peak_spec, orig_spec);
          scr_map.insert(std::pair<DoubleReal, String>(score, all_edge_annots[i][j].sequence.toString()));
        }
        //extract the best n edge annotations to be used in a full combinatorial trial
        for(std::multimap<DoubleReal, String>::reverse_iterator r_it = scr_map.rbegin(); r_it!=scr_map.rend(); ++r_it)
        {
          selected_edge_annots[i].push_back(r_it->second);
          if(selected_edge_annots[i].size()==3)
          {
            break;
          }
        }

        scr_map.clear();
        std::cout<<"i: "<<i<<"  sel_edge_size:  "<<selected_edge_annots[i].size()<<std::endl;
      }

      getAllCombinations_(selected_edge_annots, all_filtered_candidates);
      //std::cout << "TESTSequence: " << all_filtered_candidates[1].sequence << std::endl;
      all_edge_annots.clear();
    }//endof running over all predicted paths

    //now in a second run score all filtered candidates
    for (std::vector<IdAnnot>::iterator filtered_it = all_filtered_candidates.begin(); filtered_it != all_filtered_candidates.end(); ++filtered_it)
    {
      RichPeakSpectrum tmp;
      std::cout << "Sequence: " << filtered_it->sequence << std::endl;
      t_gen.getSpectrum(tmp, (*filtered_it).sequence, 1);
      t_gen.addPeaks(tmp, (*filtered_it).sequence, Residue::AIon);
      t_gen.addPeaks(tmp, (*filtered_it).sequence, Residue::BIon, 2);
      t_gen.addPeaks(tmp, (*filtered_it).sequence, Residue::YIon, 2);

      //transform the RichPeakSpectrum into a PeakSpectrum
      PeakSpectrum tmp_peak_spec;
      for (RichPeakSpectrum::iterator spec_it = tmp.begin(); spec_it != tmp.end(); ++spec_it)
      {
        Peak1D tmp_peak;
        tmp_peak.setIntensity(spec_it->getIntensity());
        tmp_peak.setMZ(spec_it->getMZ());
        tmp_peak_spec.push_back(tmp_peak);
      }

      //DoubleReal score = ZhangSimilarityScore()(tmp_peak_spec, orig_spec);
      DoubleReal score = SpectrumAlignmentScore()(tmp_peak_spec, orig_spec);
      //DoubleReal score = SteinScottImproveScore()(tmp_peak_spec, orig_spec);
      //DoubleReal score = BinnedSharedPeakCount()(BinnedSpectrum(0.5, 1, tmp_peak_spec), BinnedSpectrum(0.5, 1, orig_spec));
      ranked_candidates.insert(std::pair<DoubleReal, IdAnnot>(score, *filtered_it));
    }
  }


  void IdEval::smart_rescoring2(const IdSetup::AASeqVecMap &annot_map, const DoubleReal &precision, const std::vector<UIntVec> &result_masses,
          std::set<ScoreAnnotPair > &ranked_candidates, const PeakSpectrum & orig_spec, const BayesScoring &scr_func)
    {
      std::vector<IdAnnot>all_candidates;
      std::vector<std::vector<IdAnnot> > all_edge_annots;

      //map to each edge (i,j) the map of possible annotations with their score as key
      std::map<IndexPair, std::multimap<DoubleReal, AASequence> > edge_map;

      //run over all predictions
      for (std::vector<UIntVec>::const_iterator cand_it = result_masses.begin(); cand_it != result_masses.end(); ++cand_it)
      {
        UInt edge_start=0;
        for(UIntVec::const_iterator edge_it= cand_it->begin(); edge_it!=cand_it->end(); ++edge_it)
        {
          //if this edge has already been evaluated skip it
          if(edge_map.find(std::pair<UInt, UInt>(edge_start, *edge_it)) ==edge_map.end())
          {
            std::multimap<DoubleReal, AASequence>single_edge_map;
            //selected_edge_annots.assign(cand_it->size(), StringVec());
            all_candidates.clear();
            //all combinations of permutations are appended to all_candidates
            generateAllAnnotations(annot_map, precision, std::vector<UInt>(1,*edge_it), all_candidates);
            //all_edge_annots.push_back(all_candidates);

            //TODO  this is a hack to prevent resolution of last edge in tryptic peptides
            //all_edge_annots.back().erase(all_edge_annots.back().begin(), all_edge_annots.back().begin()+2);
//            if (edge_it+1 == cand_it->end())
//            {
//              all_candidates.erase(all_candidates.begin(), all_candidates.begin()+2);
//            }

            for(Size i=0; i<all_candidates.size();++i)
            {
              std::cout<<"Cand Sequence:  "<<all_candidates[i].sequence<<std::endl;
              AASequence sequence = all_candidates[i].sequence;
              DoubleReal support=0;
              //support+=rescoreGap(sequence, edge_start*precision, (edge_start+*edge_it)*precision, orig_spec, scr_func, 0.5);
              for(Size pos=1; pos<=sequence.size(); ++pos)
              {
                DoubleReal start_mass=(precision*edge_start) + sequence.getPrefix(pos).getMonoWeight(Residue::Internal);
                support+=scr_func.getSupportIonCount(orig_spec, start_mass);
              }
              single_edge_map.insert(std::pair<DoubleReal, AASequence>(support, sequence));
            }
            edge_map.insert(std::pair<IndexPair, std::multimap<DoubleReal, AASequence> >(IndexPair(edge_start, *edge_it), single_edge_map));
          }
          edge_start+=*edge_it;
        }

        edge_start=0;

        std::vector<std::multimap<DoubleReal, AASequence>::const_reverse_iterator>begin_iters, end_iters, iters;

        //now generate the combinations for one candidate (for each edge only the best n candidates are considered)
        for(UIntVec::const_iterator edge_it= cand_it->begin(); edge_it!=cand_it->end(); ++edge_it)
        {
          begin_iters.push_back(edge_map[IndexPair(edge_start, *edge_it)].rbegin());
          std::multimap<DoubleReal, AASequence>::const_reverse_iterator end_iter=begin_iters.back();
          int count=0;
          Size N=5;
          while(count++!=N && end_iter!=edge_map[IndexPair(edge_start, *edge_it)].rend())
          {
            ++end_iter;
          }
          end_iters.push_back(end_iter);
          edge_start+=*edge_it;
        }

        iters=begin_iters;

        //terminate when the last iterator has reached the end of the last vector
        while (iters.back() != end_iters.back())
        {
          AASequence tmp_sequence;
          DoubleReal score=0;
          UInt i = 0;
          while (i < begin_iters.size())
          {
            tmp_sequence+=iters[i]->second;
            score+=iters[i]->first;
            if (++iters[i] == end_iters[i])
            {
              if(i!=begin_iters.size()-1)
              {
                iters[i] = begin_iters[i];
              }
              ++i;
            }
            else
            {
              while (++i < begin_iters.size())
              {
                tmp_sequence+=iters[i]->second;
                score+=iters[i]->first;
              }
              break;
            }
          }
          std::cout << "in combination:  " << tmp_sequence << "   --   "<<score<<std::endl;
          ranked_candidates.insert(ScoreAnnotPair(score, IdAnnot(AASequence(tmp_sequence))));
        }        

//
//          AASequence tmp_sequence;
//          Size iter_index=0;
//          DoubleReal score=0;
//
//          //now add the current combination to the output vector
//          while (iter_index<iters.size())
//          {
//            tmp_sequence += iters[iter_index]->second;
//            score+=iters[iter_index]->first;
//            ++iter_index;
//          }
//          std::cout << "in combination:  " << tmp_sequence << std::endl;
//          ranked_candidates.insert(ScoreAnnotPair(score, IdAnnot(AASequence(tmp_sequence))));
//        }
      }
    }  

  void IdEval::getCIDSpectrum(PeakSpectrum& spec, const AASequence& sequence, Size charge, DoubleReal prefix, DoubleReal suffix)
  {
    static  DoubleReal min_mz =100;
    static  DoubleReal max_mz =2000;
    static Size max_isotope=2;
    static DoubleReal h2o_mass = EmpiricalFormula("H2O").getMonoWeight();
    static DoubleReal nh3_mass = EmpiricalFormula("NH3").getMonoWeight();
    static DoubleReal co_mass = EmpiricalFormula("CO").getMonoWeight();
    Peak1D p;
    DoubleReal b_pos(0 + prefix);
    DoubleReal y_pos(h2o_mass + suffix);
    bool b_H2O_loss(false), b_NH3_loss(false), y_NH3_loss(false);

    for (Size i = 0; i != sequence.size() - 1; ++i)
    {

      char aa = sequence[i].getOneLetterCode().c_str()[0];
      char aa2 = sequence[sequence.size()-i-1].getOneLetterCode().c_str()[0];
      b_pos = sequence.getPrefix(i+1).getMonoWeight(Residue::Internal);
      y_pos = sequence.getSuffix(i+1).getMonoWeight(Residue::Internal);

      for (Size z = 1; z <= charge && z < 3; ++z)
      {
        // b-ions
        if (b_pos >= min_mz && b_pos <= max_mz)
        {
          for (Size j = 0; j != max_isotope; ++j)
          {
            if (z == 1 /*|| b_pos > MIN_DOUBLE_MZ*/)
            {
              p.setPosition((b_pos + (DoubleReal)z * Constants::PROTON_MASS_U + (DoubleReal)j + Constants::NEUTRON_MASS_U)/(DoubleReal)z);
              IsotopeDistribution dist=sequence.getPrefix(i+1).getFormula().getIsotopeDistribution(max_isotope);
              //p.setIntensity(isotope_distributions_[(Size)b_pos][j] * 0.8 / (z*z));
              p.setIntensity(dist.getContainer()[j].second * 0.8 / (z*z));
              spec.push_back(p);
            }
          }
        }

        // b-ion losses
        if (b_pos - h2o_mass > min_mz && b_pos - h2o_mass < max_mz)
        {
          if (b_H2O_loss || aa == 'S' || aa == 'T' || aa == 'E' || aa == 'D')
          {
            b_H2O_loss = true;
            p.setPosition((b_pos + z * Constants::PROTON_MASS_U - h2o_mass)/z);
            p.setIntensity(0.02 / (DoubleReal)(z*z));
            if (z == 1/* || b_pos > MIN_DOUBLE_MZ*/)
            {
              spec.push_back(p);
            }
          }
          if (b_NH3_loss || aa == 'Q' || aa == 'N' || aa == 'R' || aa == 'K')
          {
            b_NH3_loss = true;
            p.setPosition((b_pos + z * Constants::PROTON_MASS_U - nh3_mass)/z);
            p.setIntensity(0.02 / (DoubleReal)(z*z));

            if (z == 1/* || b_pos > MIN_DOUBLE_MZ*/)
            {
              spec.push_back(p);
            }
          }
        }

        // a-ions only for charge 1
        if (z == 1)
        {
          if (b_pos - co_mass > min_mz && b_pos - co_mass < max_mz)
          {
            // a-ions
            p.setPosition((b_pos + z * Constants::PROTON_MASS_U - co_mass)/(DoubleReal)z);
            p.setIntensity(0.1f);
            spec.push_back(p);
          }
        }

        if (y_pos > min_mz && y_pos < max_mz)
        {
          // y-ions
          for (Size j = 0; j != max_isotope; ++j)
          {
            if (z == 1/* || y_pos > MIN_DOUBLE_MZ*/)
            {
              p.setPosition((y_pos + (DoubleReal)z * Constants::PROTON_MASS_U + (DoubleReal)j * Constants::NEUTRON_MASS_U)/(DoubleReal)z);
              IsotopeDistribution dist=sequence.getSuffix(i+1).getFormula().getIsotopeDistribution(max_isotope);
              p.setIntensity(dist.getContainer()[j].second / (DoubleReal)(z*z));
              spec.push_back(p);
            }
          }

          // H2O loss
          p.setPosition((y_pos + z * Constants::PROTON_MASS_U - h2o_mass)/(DoubleReal)z);
          p.setIntensity(0.1 / (DoubleReal)(z*z));
          if (aa2 == 'Q') // pyroglutamic acid formation
          {
            p.setIntensity(0.5f);
          }
          if (z == 1/* || y_pos > MIN_DOUBLE_MZ*/)
          {
             spec.push_back(p);
          }

          // NH3 loss
          if (y_NH3_loss || aa2 == 'Q' || aa2 == 'N' || aa2 == 'R' || aa2 == 'K')
          {
            y_NH3_loss = true;
            p.setPosition((y_pos + z * Constants::PROTON_MASS_U - nh3_mass)/(DoubleReal)z);
            p.setIntensity(0.1 / (DoubleReal)(z*z));

            if (z == 1 /*|| y_pos > MIN_DOUBLE_MZ*/)
            {
              spec.push_back(p);
            }
          }
        }
      }
    }

    // if Q1 abundant loss of water -> pyroglutamic acid formation
    //if (sequence[0] == 'Q' && prefix == 0 && suffix == 0)
    {
      /*
      for (PeakSpectrum::Iterator it = spec.begin(); it != spec.end(); ++it)
      {
        it->setIntensity(it->getIntensity() * 0.5);
      }*/

      /*
      for (Size j = 0; j != max_isotope; ++j)
      {
        p.setPosition((precursor_weight + charge - 1 + j)/(DoubleReal)charge);
        p.setIntensity(isotope_distributions_[(Int)p.getPosition()[0]][j] * 0.1);
        spec.push_back(p);
      }
      */
    }


    spec.sortByPosition();

    return;
  }

  DoubleReal IdEval::scorePSM(const PeakSpectrum & spec, const AASequence & sequence, const std::set<IdSetup::ion_type>scored_types, /*const PriorTable & priors,*/ DoubleReal delta)
  {
    const static Size number_sectors = 3;
    static std::vector<std::vector<DoubleReal> >priors(number_sectors);
    if(priors.empty())
    {
      DoubleReal prior;
      for(Size sector=0; sector<number_sectors; ++sector)
      {
        for(Size type=0; type<=IdSetup::LAST_TYPE; ++type)
        {
          if(type==IdSetup::BIon || type==IdSetup::YIon)
            prior=1;
          else if(type==IdSetup::BIon2 || type==IdSetup::YIon2)
            prior=0.5;
          else if(type==IdSetup::CIon || type==IdSetup::ZIon)
            prior=0.2;
          else //neutral losses
            prior=0.5;

          priors[sector].push_back(prior);
        }
      }
    }

    static const DoubleReal delta_for_isotopes=0.1;
    static const DoubleReal weight_for_isotopes=0.5;
    static const DoubleReal secondary_peak_weight = 0.8; //a penalty factor is a peak is a secondary one
    static const DoubleReal missing_peak_penalty = 2.0; //the penalty factor for the missing peak. penalty is -1 * missing_peak_penalty * prior

    //for each of the supported ion types we check for witness peaks and score them
    DoubleReal parent_mass=sequence.getMonoWeight(Residue::Full,1);
    DoubleReal total_score=0;
    for(Size pos=1; pos<sequence.size(); ++pos)
    {
      DoubleReal prm=sequence.getPrefix(pos).getMonoWeight(Residue::Internal);
      Size sector=std::min(number_sectors-1, (Size)floor(number_sectors * prm/parent_mass));

      std::set<IdSetup::ion_type>::const_iterator type_it;
      for(type_it=scored_types.begin(); type_it!=scored_types.end(); ++type_it)
      {
        DoubleReal score=0;
        DoubleReal offset_pos=IdSetup::getFragmentMass(*type_it, prm, parent_mass);
        Size peak_index = spec.findNearest(offset_pos);
        //check whether peak is in allowed range
        DoubleReal mass_error = fabs(spec[peak_index].getMZ()-offset_pos);

        if(mass_error<delta)
        {
          bool is_primary=false, is_secondary=false;
          DoubleReal parent_intensity=0.0, child_intensity=0.0;
          //check whether peak is lone, primary or secondary peak
          Size index=peak_index;
          Size charge=1;
          if(*type_it==IdSetup::BIon2 || *type_it==IdSetup::YIon2)
          {
            charge=2;
          }
          DoubleReal target_pos=offset_pos - Constants::PROTON_MASS_U/charge;
          //check left side of the peak
          while(index>=0)
          {
            if(fabs(spec[index].getMZ()-target_pos)<delta_for_isotopes)
            {
              if(spec[index].getIntensity()>spec[peak_index].getIntensity())
              {
                is_secondary=true;
                parent_intensity=spec[index].getIntensity();
              }
              break;
            }
            else if(spec[index].getMZ()<target_pos-delta_for_isotopes)
            {
              break;
            }
            --index;
          }
          //check right side of the peak
          index=peak_index;
          target_pos=offset_pos + Constants::PROTON_MASS_U/charge;
          while(index<spec.size())
          {
            if(fabs(spec[index].getMZ()-target_pos)<delta_for_isotopes)
            {
              if(spec[index].getIntensity()<spec[peak_index].getIntensity())
              {
                is_primary=true;
                child_intensity=spec[index].getIntensity();
              }
              break;
            }
            else if(spec[index].getMZ()>target_pos+delta_for_isotopes)
            {
              break;
            }
            ++index;
          }

          //sum up the score respecting the peak status (lone, primary, secondary), the mass_error and the prior
          score=spec[peak_index].getIntensity()*priors[sector][*type_it];
          score*= ((delta-mass_error)/delta);
          if(is_primary)
          {
            score+=weight_for_isotopes*child_intensity;
          }
          else if(is_secondary) //penalty if its a secondary peak
          {
            score=secondary_peak_weight*score;
          }
        }

        //if no peak was found penalize it
        else
        {
          score=(-1)*missing_peak_penalty*priors[sector][*type_it];
        }
        total_score+=score;
      }
    }
    return total_score;
  }


}//end of namespace

