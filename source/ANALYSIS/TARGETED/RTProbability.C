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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/RTProbability.h>
#include <OpenMS/SIMULATION/RTSimulation.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#include <fstream>
#include <limits>
#include <gsl/gsl_cdf.h>
#include <OpenMS/SIMULATION/SimTypes.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define RTPROB_DEBUG
//#undef  RTPROB_DEBUG

namespace OpenMS 
{
  RTProbability::RTProbability()
    : DefaultParamHandler("RTProbability"),
      sigma_(0.),
      mu_(0.)

  {
    defaults_.setValue("number_of_bins", 40, "Number of bins used for the fitting, if sparse datasets are used, this number should be smaller", StringList::create("advanced"));
    defaults_.setValue("rt_settings:min_rt", 960., "Minimal RT in the experiment (in seconds)");
    defaults_.setMinFloat("rt_settings:min_rt",0.);
    defaults_.setValue("rt_settings:max_rt", 3840., "Maximal RT in the experiment (in seconds)");
    defaults_.setMinFloat("rt_settings:min_rt",1.);
    defaults_.setValue("rt_settings:rt_step_size", 30., "Time between two consecutive spectra (in seconds)");
    defaults_.setMinFloat("rt_settings:min_rt",1.);
    defaults_.setValue("gauss_amplitude", 0.08, "Initial Gauss amplitude");
    defaults_.setMinFloat("gauss_amplitude",0.);
    defaults_.setValue("gauss_mean", -1., "Initial Gauss mean");
    defaults_.setValue("gauss_sigma", -3., "Initial Gauss sigma");
    
  
#ifdef RTPROB_DEBUG
    defaults_.setValue("rev_filename", "rev", "bla", StringList::create("advanced"));
    defaults_.setValue("fwd_filename", "fwd", "bla", StringList::create("advanced"));
#endif

    defaultsToParam_();
  }

  RTProbability::RTProbability(const RTProbability& rhs)
    : DefaultParamHandler(rhs),
      sigma_(rhs.sigma_),
      mu_(rhs.mu_)
  {

  }

  RTProbability::~RTProbability()
  {

  }

  DoubleReal RTProbability::getRTProbability(DoubleReal min_obs_rt,DoubleReal max_obs_rt,DoubleReal theo_rt)
  {
    // first adapt gaussian RT error distribution to a normal distribution with \mu = 0
    Int theo_scan = getScanNumber_(theo_rt);
    Int obs_scan_begin = getScanNumber_(min_obs_rt)-1;
    Int obs_scan_end = getScanNumber_(max_obs_rt)+1;

    if(obs_scan_begin == -1 || obs_scan_end == -1)
      {
        std::cout << "Probably an error occured during RTProb-calc: scan = -1: "
                  << obs_scan_begin << " "<< obs_scan_end << std::endl;
        return 0.;
      }
    
    obs_scan_begin -= mu_;
    obs_scan_end -= mu_;
    
    DoubleReal x1 = theo_scan - obs_scan_end;
    DoubleReal x2 = theo_scan - obs_scan_begin;
    
    DoubleReal prob;
    if(x2 > x1) prob = gsl_cdf_gaussian_P(x2,sigma_) - gsl_cdf_gaussian_P(x1,sigma_);
    else prob = gsl_cdf_gaussian_P(x1,sigma_) -  gsl_cdf_gaussian_P(x2,sigma_);
    if((prob < 0.) || (obs_scan_begin == obs_scan_end))
      {
        std::cerr << min_obs_rt << " "<< obs_scan_begin << " " << max_obs_rt << " "<< obs_scan_end << " "
                  << theo_rt << " " << theo_scan << " " << mu_ << " "<< x1 << " "<<x2 << " "<< prob<<std::endl;
      }
    return prob;
  }


  void RTProbability::learnGaussian(FeatureMap<>& features,String rt_model_path,bool higherScoreBetter,DoubleReal min_score)
  {
    Size number_of_bins(param_.getValue("number_of_bins"));
    std::cout << number_of_bins <<" bins, min_score "<< min_score<<std::endl;
    std::vector<AASequence> peptide_sequences;
    std::vector<std::pair<Int,Int> > measured_rts;
    std::vector<DoubleReal> predicted_rts;
    // filter for save peptide ids
    for(Size f = 0; f < features.size();++f)
      {
        std::vector<PeptideIdentification> & pep_ids = features[f].getPeptideIdentifications();
        DoubleReal max_score = higherScoreBetter ? 0. : 1.;
        AASequence max_seq = "";
        for(Size pi = 0; pi < pep_ids.size(); ++pi)
          {
            if(!higherScoreBetter && pep_ids[pi].isHigherScoreBetter())
              {
                std::cout << "Attention: Different id score types mixed up!"<<std::endl;
                break;
              }
            for(Size ph = 0; ph < pep_ids[pi].getHits().size();++ph)
              {
                //  std::cout << pep_ids[pi].getHits()[ph].getScore() << "\t";
                if((higherScoreBetter && pep_ids[pi].getHits()[ph].getScore() > max_score)||
                   (!higherScoreBetter && pep_ids[pi].getHits()[ph].getScore() < max_score))
                  {
                    // std::cout << pep_ids[pi].getHits()[ph].getScore() << "\t"
                    //           << pep_ids[pi].getHits()[ph].getSequence() << "\n";
                    max_score = pep_ids[pi].getHits()[ph].getScore();
                    max_seq = pep_ids[pi].getHits()[ph].getSequence();
                  }
              }
          }
        //     std::cout << "pep_ids.size() "<<pep_ids.size()<<"max_score "<< max_score << std::endl;
        if(max_seq != "" &&
           ( (higherScoreBetter && (max_score > min_score)) || (!higherScoreBetter && (max_score < min_score)) ))
          {
            peptide_sequences.push_back(max_seq);
            DoubleReal rt_begin = features[f].getConvexHull().getBoundingBox().minPosition()[0];
            DoubleReal rt_end =   features[f].getConvexHull().getBoundingBox().maxPosition()[0];
            Int scan_begin = getScanNumber_(rt_begin);
            Int scan_end = getScanNumber_(rt_end);
            if(scan_begin != -1 && scan_end != -1) measured_rts.push_back(std::make_pair(scan_begin,scan_end));
//             std::cout << "min_rt "<<min_rt<<" rt_step_size "<< rt_step_size <<" max_rt "<<max_rt
//                       << " rt "<<features[f].getRT()<<" scan "<<scan <<std::endl;
//             std::cout << "max_score "<<max_score << std::endl;
          }
      }
    // make rt predictions
    SimRandomNumberGenerator rnd_gen;
//     rnd_gen.biological_rng = gsl_rng_alloc(gsl_rng_mt19937);
//     gsl_rng_set(rnd_gen.biological_rng, time(0));
//     gsl_rng_set(rnd_gen.technical_rng, time(0));
    RTSimulation rt_sim(rnd_gen);
    Param param;
    param.setValue("HPLC:model_file",rt_model_path);
    rt_sim.setParameters(param);
    rt_sim.wrapSVM(peptide_sequences,predicted_rts);
    mu_ = 0.;
    sigma_ = 0.;
    // calculate mean
    std::vector<Int> predicted_scans;
    std::vector<Int> diffs;
    Size invalid_scans = 0;
    for(Size s = 0; s < predicted_rts.size();++s)
      {
        std::cout << predicted_rts[s] << " "<< measured_rts[s].first << " "<<measured_rts[s].second << std::endl;
        Int scan = getScanNumber_(predicted_rts[s]);
        if(scan == -1) ++invalid_scans;
        predicted_scans.push_back(scan);
      }

    std::ofstream out("rt_fehler");
    for(Size s = 0; s < predicted_scans.size();++s)
      {
        //				std::cout << predicted_rts[s] << " "<< measured_rts[s] << std::endl;
        if(predicted_scans[s] != -1)
          {
            if(predicted_scans[s] < measured_rts[s].first) // smaller than feature's min_rt
              {
                mu_ +=  (predicted_scans[s]-measured_rts[s].first);
                diffs.push_back(predicted_scans[s] - measured_rts[s].first);
                out << predicted_scans[s] - measured_rts[s].first<<std::endl;
              }
            else if(predicted_scans[s] < measured_rts[s].second) // between feature's min_rt and max_rt
              {
                //                mu_ +=  0;
                diffs.push_back(0);
                out << 0<<std::endl;
              }
            else // bigger than feature's max_rt
              {
                mu_ +=  (predicted_scans[s] - measured_rts[s].second);
                diffs.push_back(predicted_scans[s] - measured_rts[s].second);
                out << predicted_scans[s] - measured_rts[s].second <<std::endl;
              }
          }
      }
    out.close();
    mu_ /= (DoubleReal)(predicted_scans.size()-invalid_scans);
		
    std::cout << "rt_diff_mean "<<mu_ << std::endl;
    // calculate std
    for(Size s = 0; s < predicted_scans.size();++s)
      {
        if(predicted_scans[s] != -1)
          {
            DoubleReal diff;
            if(predicted_scans[s] < measured_rts[s].first) // smaller than feature's min_rt
              {
                diff =  (predicted_scans[s]-measured_rts[s].first) - mu_;
              }
            else if(predicted_scans[s] < measured_rts[s].second) // between feature's min_rt and max_rt
              {
                diff = 0 - mu_;
              }
            else // bigger than feature's max_rt
              {
                diff =  (predicted_scans[s] - measured_rts[s].second) - mu_;
              }
            sigma_ += diff*diff;
          }
      }
    sigma_ /= (DoubleReal)(predicted_scans.size()-invalid_scans-1);
    sigma_ = sqrt(sigma_);
    std::cout << "rt_diff_std "<< sigma_ <<std::endl;

    // TODO: compute Histogramm and use GaussFitter
    // get the range of the diffs
    Int max(std::numeric_limits<Int>::min()), min(std::numeric_limits<Int>::max());
    for (std::vector<Int>::const_iterator it = diffs.begin(); it != diffs.end(); ++it)
      {
        if (*it > max)
          {
            max = *it;
          }
        if (*it < min)
          {
            min = *it;
          }
      }
    // make the binning
    Int diff = max - min;
#ifdef RTPROB_DEBUG
    std::cout << "min "<< min << " max "<< max << " diff " << diff << std::endl;
#endif
    Size number_of_diffs(diffs.size());
    std::vector<DoubleReal> binned_diffs(number_of_diffs);
    for (std::vector<Int>::const_iterator it = diffs.begin(); it != diffs.end(); ++it)
      {
        Size bin = (Size)((DoubleReal)(*it - min) / diff * (DoubleReal)(number_of_bins - 1));
        binned_diffs[bin] += 1.;
      }


    std::vector<DPosition<2> > diff_data(number_of_bins);;
    // normalize to \sum = 1 and store in diff_data
    for (Size i=0; i < number_of_bins; ++i)
      {
        binned_diffs[i] /= (DoubleReal) number_of_diffs * 3.;
        DPosition<2> pos;
        pos.setX(min + diff*(DoubleReal)i / (DoubleReal)number_of_bins);
#ifdef RTPROB_DEBUG
        std::cout << min << " + "<<(DoubleReal)i / (DoubleReal)number_of_bins << " * "<< i << " = "
                  << pos.getX()<<std::endl;
#endif
        pos.setY(binned_diffs[i]);
        diff_data.push_back(pos);
      }

    Math::GaussFitter gf;
    Math::GaussFitter::GaussFitResult result_1st;
    result_1st.A = param_.getValue("gauss_amplitude");//1./(sigma_*sqrt(2.*Constants::PI));//0.028; //gauss_A; //0.06;
    result_1st.x0 = mu_;//param_.getValue("gauss_mean");// -1.;//mu_;//0.008;//gauss_x0; //0.7;
    result_1st.sigma = sigma_;//param_.getValue("gauss_sigma"); //3.;//sigma_;//0.25;//gauss_sigma; //0.5;
    gf.setInitialParameters(result_1st);
#ifdef RTPROB_DEBUG
    std::cerr << "Initial Gauss guess: A=" << result_1st.A << ", mu=" << result_1st.x0 << ", sigma=" << result_1st.sigma << std::endl;
#endif
    Math::GaussFitter::GaussFitResult result_gauss;
    try
      {
        //result_gauss = gf.fit(diff_data);
        if (gf.getGnuplotFormula() == "")
          {
            result_gauss.A = result_1st.A;
            result_gauss.x0 = result_1st.x0;
            result_gauss.sigma = result_1st.sigma;
          }

      }
    catch(...)
      {
        std::cerr <<"fit failed"<<std::endl;
        // fit failed?
        if (gf.getGnuplotFormula() == "")
          {
            result_gauss.A = result_1st.A;
            result_gauss.x0 = result_1st.x0;
            result_gauss.sigma = result_1st.sigma;
          }
      }
    
    sigma_ = result_gauss.sigma; 
    mu_ = result_gauss.x0;
    
#ifdef RTPROB_DEBUG
    std::cerr << gf.getGnuplotFormula() << std::endl;
    String fwd_filename = "gaussian_rt_error";//param_.getValue("fwd_filename");
    if (gf.getGnuplotFormula() == "")
      {
        String formula("f(x)=" + String(result_1st.A) + " * exp(-(x - " + String(result_1st.x0) + ") ** 2 / 2 / (" + String(result_1st.sigma) + ") ** 2)");
        generateDistributionImage_(binned_diffs, formula, fwd_filename,min,diff); 
      }
    else
      {
        generateDistributionImage_(binned_diffs, gf.getGnuplotFormula(), fwd_filename,min,diff);
      }

    std::ofstream os2("rt_range_vs_pred_rt");
    
    
#endif
  }
  

  
  Int RTProbability::getScanNumber_(DoubleReal rt)
  {
    DoubleReal min_rt = param_.getValue("rt_settings:min_rt");
    DoubleReal max_rt = param_.getValue("rt_settings:max_rt");
    DoubleReal rt_step_size = param_.getValue("rt_settings:rt_step_size");
    
    if(rt > max_rt || rt < min_rt) return -1;

    Int scan = (Int)floor((rt - min_rt) / rt_step_size);
    return scan;
  }


  void RTProbability::generateDistributionImage_(const std::vector<DoubleReal>& ids,
                                                 const String& formula, const String& filename,
                                                 DoubleReal min, DoubleReal diff)
  {
    Size number_of_bins(param_.getValue("number_of_bins"));
		
    // write distribution to file
    std::ofstream o((filename + "_dist_tmp.dat").c_str());
    for (Size i = 0; i < number_of_bins; ++i)
      {
        o << min + diff*(DoubleReal)i / (DoubleReal)number_of_bins << " " << ids[i] << std::endl;
      }
    o.close();

    std::ofstream os((filename + "_gnuplot.gpl").c_str());
    os << "set terminal png" << std::endl;
    os << "set output '" << filename << "_distribution.png'" << std::endl;
    os << formula <<std::endl;
    os << "plot f(x), '" << filename << "_dist_tmp.dat' w boxes" << std::endl;
    os.close();

    system(("gnuplot " + filename + "_gnuplot.gpl").c_str());

    return;
  }
	

	void RTProbability::setGaussianParameters(DoubleReal mu, DoubleReal sigma)
	{
		sigma_ = sigma;
		mu_ = mu;
	}


}//namespace
