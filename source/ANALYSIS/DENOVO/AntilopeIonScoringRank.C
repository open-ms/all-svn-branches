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

#include <numeric>
#include <cmath>

#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringRank.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>

#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{

  //default constructor
  RankScoringFunction::RankScoringFunction() :
    DefaultParamHandler("RankScoringFunction"),
    log_odds_b_(0),
    log_odds_y_(0)
  {
    defaults_.setValue("number_of_sectors", 3, "The number of sectors");
    defaults_.setMaxInt("number_of_sectors", 5);
    defaults_.setMinInt("number_of_sectors", 1);
    defaults_.setValue("delta", 0.5, "allowed mass deviation from PRM for each peak");
    defaults_.setValue("number_of_ranks", 100, "number of considered ranks");
    defaults_.setValue("scoring_file", "../../data/scoring/rankScores.dat", "file to read/store rank scores");
    defaultsToParam_();
  }


  //custom constructor
//  RankScoringFunction::RankScoringFunction() :
//    DefaultParamHandler("RankScoringFunction")
//  {
//    defaults_.setValue("number_of_sectors", 3, "The number of sectors");
//    defaults_.setMaxInt("number_of_sectors", 5);
//    defaults_.setMinInt("number_of_sectors", 1);
//    defaults_.setValue("delta", 0.5, "allowed mass deviation from PRM for each peak");
//    defaults_.setValue("number_of_ranks", 100, "number of considered ranks");
//    defaults_.setValue("scoring_file", "../data/scoring/rankScores.bin", "file to read/store rank scores");
//    defaultsToParam_();
//
//    loadModel();
//
//  }

  //copy constructor
  RankScoringFunction::RankScoringFunction(RankScoringFunction & in) :
    DefaultParamHandler(in),
    log_odds_b_(in.log_odds_b_),
    log_odds_y_(in.log_odds_y_)
  {
    updateMembers_();
  }

  //assigment operator
  RankScoringFunction& RankScoringFunction::operator=(const RankScoringFunction& in)
  {
    if (this != &in)
    {
      DefaultParamHandler::operator =(in);
      log_odds_b_ = in.log_odds_b_;
      log_odds_y_ = in.log_odds_y_;
      updateMembers_();
    }
    return *this;
  }

  //update the member variables delta and precision
  void RankScoringFunction::updateMembers_()
  {
    delta_ = param_.getValue("delta");
    number_of_sectors_ = param_.getValue("number_of_sectors");
    max_rank_ = param_.getValue("number_of_ranks");
  }

  //the input Map has to be annotated with the true sequences
  void RankScoringFunction::generateRankProbs(const PeakMap &spec_map)
  {
    std::vector<UIntVec>b_counts(number_of_sectors_, UIntVec(max_rank_,0));
    std::vector<UIntVec>y_counts(b_counts);
    std::vector<UIntVec>noise_counts(b_counts);

    DMatrix b_probs(number_of_sectors_, DVector(max_rank_,0.0));
    DMatrix y_probs(number_of_sectors_, DVector(max_rank_,0.0));

    DVector b_ion_masses, y_ion_masses;

    for (PeakMap::const_iterator spec_it = spec_map.begin(); spec_it != spec_map.end(); ++spec_it)
    //for (RichPeakMap::const_iterator spec_it = specs.begin(); spec_it != specs.begin()+1; ++spec_it)
    {
      PeakSpectrum local_spec(*spec_it);
      IdSetup::filterSpectrumInspect(local_spec);
      //std::cout<<local_spec[0].getMZ()<< " !!!!! "<<local_spec[1].getMZ()<<std::endl;
      //std::cout<<local_spec[0].getIntensity()<< " !!!!! "<<local_spec[1].getIntensity()<<std::endl;
      AASequence annotation = spec_it->getPeptideIdentifications()[0].getHits()[0].getSequence();
      if (!annotation.isValid())
      {
        std::cout<<"jump out"<<std::endl;
        continue;
      }
      b_ion_masses.clear();
      y_ion_masses.clear();
      b_ion_masses.reserve(annotation.size());
      y_ion_masses.reserve(annotation.size());

      for(Size substr_len=1; substr_len<annotation.size();++substr_len)
      {
        b_ion_masses.push_back(annotation.getPrefix(substr_len).getMonoWeight(Residue::BIon, 1));
        y_ion_masses.push_back(annotation.getSuffix(substr_len).getMonoWeight(Residue::YIon, 1));
      }

      const DoubleReal charge = spec_it->getPrecursors()[0].getCharge();
      const DoubleReal precursor_mass = (spec_it->getPrecursors()[0].getPosition()[0]*charge) - (charge-1);

      //set intensity for each peak to its rank
      local_spec.sortByIntensity(true);
      for(Size i=0; i<local_spec.size();++i)
      {
        local_spec[i].setIntensity(i);
      }
      local_spec.sortByPosition();

      //now start simply counting
      PeakSpectrum::iterator peak_it = local_spec.begin();
      Size nr_y_ion = 0, nr_b_ion = 0;
      while (peak_it != local_spec.end() && nr_b_ion < b_ion_masses.size() && nr_y_ion < y_ion_masses.size())
      {
        Size rank = (Size)peak_it->getIntensity();
        if(rank >= max_rank_)
        {
         ++peak_it;
         continue;
        }
        Size sector = std::min(number_of_sectors_-1, (Size)(floor((peak_it->getMZ() / precursor_mass) * number_of_sectors_)));
        //std::cout<<peak_it->getMZ()<<"  "<<peak_it->getIntensity()<<"  "<<b_ion_masses[nr_b_ion]<<"  "<<y_ion_masses[nr_y_ion]<<std::endl;

        if (fabs(peak_it->getMZ() - b_ion_masses[nr_b_ion]) <= delta_)
        {
          ++b_counts[sector][rank];
          ++peak_it;
        }
        else if (fabs(peak_it->getMZ() - y_ion_masses[nr_y_ion]) <= delta_)
        {
          ++y_counts[sector][rank];
          ++peak_it;
        }
        else if (peak_it->getMZ() > y_ion_masses[nr_y_ion])
        {
          ++nr_y_ion;
        }
        else if (peak_it->getMZ() > b_ion_masses[nr_b_ion])
        {
          ++nr_b_ion;
        }
        else
        {
          ++noise_counts[sector][rank];
          ++peak_it;
        }
      }
    }

    //test output
    for (Size sector = 0; sector < number_of_sectors_; ++sector)
    {
      std::cout << "Sector " << sector << std::endl;
      for (Size rank = 0; rank < b_counts[sector].size(); ++rank)
      {
        UInt total=b_counts[sector][rank] + y_counts[sector][rank] + noise_counts[sector][rank];
        std::cout << "b-ion  " << rank << "  ::  " << (DoubleReal) b_counts[sector][rank]/total << std::endl;
        b_probs[sector][rank]=(DoubleReal) b_counts[sector][rank]/total;
      }
      for (Size rank = 0; rank < y_counts[sector].size(); ++rank)
      {
        UInt total=b_counts[sector][rank] + y_counts[sector][rank] + noise_counts[sector][rank];
        std::cout << "y-ion  " << rank << "  ::  " << (DoubleReal) y_counts[sector][rank]/total << std::endl;
        y_probs[sector][rank]=(DoubleReal) y_counts[sector][rank]/total;
      }
    }

    log_odds_b_.assign(number_of_sectors_, DVector(max_rank_,0.0));
    log_odds_y_.assign(number_of_sectors_, DVector(max_rank_,0.0));

    for(Size i=0; i<number_of_sectors_;++i)
    {
      UInt total_b = std::accumulate(b_counts[i].begin(),b_counts[i].end(),0);
      UInt total_y = std::accumulate(y_counts[i].begin(),y_counts[i].end(),0);
      UInt total_noise = std::accumulate(noise_counts[i].begin(),noise_counts[i].end(),0);
      DoubleReal background_probs_b = (DoubleReal)total_b/(total_b + total_y + total_noise);
      DoubleReal background_probs_y = (DoubleReal)total_y/(total_b + total_y + total_noise);

      for(Size rank=0; rank<max_rank_;++rank)
      {
        if(b_probs[i][rank]!=0)
        {
          log_odds_b_[i][rank]=log((b_probs[i][rank]/(1-b_probs[i][rank]))/(background_probs_b/(1-background_probs_b)));
        }
        else
        {
          log_odds_b_[i][rank]=-(std::numeric_limits<float>::max());
        }
        if(y_probs[i][rank]!=0)
        {
          log_odds_y_[i][rank]=log((y_probs[i][rank]/(1-y_probs[i][rank]))/(background_probs_y/(1-background_probs_y)));
          //log_odds_y_[i][rank]=log(y_probs[i][rank]/background_probs_y);
        }
        else
        {          
          log_odds_y_[i][rank]=-(std::numeric_limits<float>::max());
        }
      }
    }
  }

  void RankScoringFunction::writeModel()
  {
    String file_name = param_.getValue("scoring_file");
    TextFile  file;
    //store the number of sectors
    file.push_back("<NumberSectors>");
    file.push_back(number_of_sectors_);
    file.push_back("</NumberSectors>");
    std::vector<Int>tmp_sectors(1,number_of_sectors_);
    std::streampos pos= IdSetup::writeVectorToBinFile<Int>(file_name, pos, tmp_sectors);

    //store the maximum rank
    file.push_back("<NumberRanks>");
    file.push_back(max_rank_);
    file.push_back("</NumberRanks>");

    for(Size sector=0; sector<number_of_sectors_; ++sector)
    {
      file.push_back("<Sector " + String(sector) + ">");
      std::cout<<"writing sector: "<<sector<<std::endl;
      file.push_back("<Prob. B-ion>");
      for(Size index=0; index<log_odds_b_[sector].size(); ++index)
      {
        file.push_back(log_odds_b_[sector][index]);
      }
      file.push_back("</Prob. B-ion>");
      file.push_back("<Prob. Y-ion>");
      for(Size index=0; index<log_odds_b_[sector].size(); ++index)
      {
        file.push_back(log_odds_y_[sector][index]);
      }
      file.push_back("</Prob. Y-ion>");
      file.push_back("</Sector " + String(sector) + ">");
    }
    file.store(file_name);
  }

  void RankScoringFunction::loadModel()
  {
    String file_name = param_.getValue("scoring_file");
    TextFile infile(file_name);
    log_odds_b_.clear();
    log_odds_y_.clear();

    TextFile::iterator left_marker = infile.search("<NumberSectors>");
    TextFile::iterator right_marker = infile.search("</NumberSectors>");
    if(left_marker!=infile.end() && right_marker!=infile.end())
    {
      ++left_marker;
      number_of_sectors_= left_marker->toInt();
    }
    left_marker = infile.search(right_marker, "<NumberRanks>");
    right_marker = infile.search(right_marker, "</NumberRanks>");
    if(left_marker!=infile.end() && right_marker!=infile.end())
    {
      ++left_marker;
      max_rank_= left_marker->toInt();
    }
    log_odds_b_.assign(number_of_sectors_,DVector(max_rank_));
    log_odds_y_.assign(number_of_sectors_,DVector(max_rank_));

    for(Size sector=0; sector<number_of_sectors_; ++sector)
    {
      left_marker = infile.search(right_marker, "<Sector " + String(sector) + ">");
      right_marker = infile.search(right_marker, "</Sector " + String(sector) + ">");
      if(left_marker!=infile.end() && right_marker!=infile.end())
      {
        left_marker = infile.search(left_marker, "<Prob. B-ion>");
        right_marker = infile.search(left_marker, "</Prob. B-ion>");
        if(left_marker!=infile.end() && right_marker!=infile.end())
        {
          Size index=0;
          while(++left_marker!=right_marker)
          {
            log_odds_b_[sector][index]=left_marker->toDouble();
            ++index;
          }
        }
        left_marker = infile.search(right_marker, "<Prob. Y-ion>");
        right_marker = infile.search(right_marker, "</Prob. Y-ion>");
        if(left_marker!=infile.end() && right_marker!=infile.end())
        {
          Size index=0;
          while(++left_marker!=right_marker)
          {
            log_odds_y_[sector][index]=left_marker->toDouble();
            ++index;
          }
        }
      }
    }
    infile.store(file_name);
  }

  DoubleReal RankScoringFunction::getRankScores(DVector &b_scores, DVector &y_scores, const PeakSpectrum &spec)
  {
    //set intensity for each peak to its rank
    PeakSpectrum local_spec(spec);

    b_scores.clear();
    y_scores.clear();
    b_scores.reserve(local_spec.size());
    y_scores.reserve(local_spec.size());

    const DoubleReal charge = local_spec.getPrecursors()[0].getCharge();
    const DoubleReal precursor_mass = (local_spec.getPrecursors()[0].getPosition()[0]*charge) - (charge-1);

    local_spec.sortByIntensity(true);
    for(Size i=0; i<local_spec.size();++i)
    {
      local_spec[i].setIntensity(std::min(i,max_rank_-1));
    }
    local_spec.sortByPosition();

    for(PeakSpectrum::ConstIterator peak_it = local_spec.begin(); peak_it!=local_spec.end();++peak_it)
    {
      Size sector = std::min(number_of_sectors_-1, (Size)(floor((peak_it->getMZ() / precursor_mass) * number_of_sectors_)));
      b_scores.push_back(log_odds_b_[sector][(Size)peak_it->getIntensity()]);
      std::cout<<"b_scores "<<sector<<"  "<<b_scores.back()<<"  "<<(Size)peak_it->getIntensity()<<std::endl;
      y_scores.push_back(log_odds_y_[sector][(Size)peak_it->getIntensity()]);
      std::cout<<"y_scores "<<sector<<"  "<<y_scores.back()<<std::endl;
    }
  }


}//namespace
