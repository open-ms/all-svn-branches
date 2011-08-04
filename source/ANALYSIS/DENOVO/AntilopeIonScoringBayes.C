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


//STL includes
#include <map>
#include <list>
#include <bitset>
#include <string>

#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringBayes.h>

//OpenMS includes
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>


#define PEPNOVO 1
#define DEBUG 1


namespace OpenMS
{
 //-----------------------------------------------------------------------------
//----------------- Definitions for BayesScoringSmile class--------------
//-----------------------------------------------------------------------------


//default constructor
BayesScoring::BayesScoring() :  
  number_intensity_levels_(0),
  intensity_limits_(0),
  delta_(0.5)
{
}

//custom constructor
BayesScoring::BayesScoring(String parameter_file) :  
  number_intensity_levels_(0),
  intensity_limits_(0),
  delta_(0.5)
{
  loadFile(parameter_file);
}

//copy constructor
BayesScoring::BayesScoring(const BayesScoring & in) :
  number_intensity_levels_(in.number_intensity_levels_),
  intensity_limits_(in.intensity_limits_),  
  precision_(in.precision_),
  delta_(in.delta_),
  number_of_sectors_(in.number_of_sectors_),
  networks_(in.networks_)
{  
}

//assigment operator
BayesScoring& BayesScoring::operator=(const BayesScoring& in)
{
  if (this != &in)
  {
    number_intensity_levels_=in.number_intensity_levels_;
    intensity_limits_=in.intensity_limits_;
    precision_= in.precision_;
    delta_ = in.delta_;
    number_of_sectors_ = in.number_of_sectors_,
    networks_= in.networks_;
  }
  return *this;
}

//load the statistics parameter file
UInt BayesScoring::loadFile(String bayes_info_file)
{
  //delete all previously loaded models
  networks_.clear();

  TextFile info_file;
  info_file.load(bayes_info_file);

  TextFile::iterator left_marker = info_file.search("<NetworkList>");
  TextFile::iterator right_marker = info_file.search("</NetworkList>");  
  if(left_marker>=right_marker)
  {
    std::cerr<<"Parse error"<<std::endl;
    return 1;
  }
  left_marker = info_file.search(left_marker, "<NumSectors>");
  right_marker = info_file.search(left_marker, "</NumSectors>");
  if(++left_marker>=right_marker)
  {
    std::cerr<<"Parse error"<<std::endl;
    return 1;
  }
  number_of_sectors_=(*left_marker).toInt();

  left_marker = info_file.search(left_marker, "<NumIntenLevels>");
  right_marker = info_file.search(left_marker, "</NumIntenLevels>");
  if(++left_marker>=right_marker)
  {
    std::cerr<<"Parse error"<<std::endl;
    return 1;
  }
  number_intensity_levels_=(*left_marker).toInt();

  std::vector<Size>sectors(number_of_sectors_);
  std::vector<String>filenames(number_of_sectors_);

  left_marker = info_file.search(left_marker, "<Network>");
  right_marker = info_file.search(left_marker, "</Network>");

  for(Size sec=0; sec<number_of_sectors_; ++sec)
  {
    if(left_marker>=right_marker)
    {
      std::cerr<<"Parse error"<<std::endl;
      return 1;
    }

    left_marker = info_file.search(left_marker, "<Sector>");
    right_marker = info_file.search(left_marker, "</Sector>");
    if(++left_marker<right_marker)
    {
      sectors[sec]=(*left_marker).toInt();
    }
    else
    {
      std::cerr<<"Parse error"<<std::endl;
      return 1;
    }

    left_marker = info_file.search(left_marker, "<NetworkFile>");
    right_marker = info_file.search(left_marker, "</NetworkFile>");
    if(++left_marker<right_marker)
    {
      filenames[sec] = *left_marker;
    }
    else
    {
      std::cerr<<"Parse error"<<std::endl;
      return 1;
    }
    left_marker = info_file.search(left_marker, "<Network>");
    right_marker = info_file.search(left_marker, "</Network>");
  }

  networks_.assign(number_of_sectors_, BayesianNetwork());
  for(Size i=0; i<number_of_sectors_; ++i)
  {
    networks_[sectors[i]].load(filenames[i]);
  }
  selected_types_.clear();
  selected_types_.assign(number_of_sectors_, std::set<IdSetup::ion_type>());
  std::set<String>tmp_types;
  for(Size sector=0; sector<number_of_sectors_; ++sector)
  {
    tmp_types.clear();
    networks_[sector].getIdentifiers(tmp_types);
    std::set<String>::const_iterator t_it;
    for(t_it=tmp_types.begin(); t_it!=tmp_types.end(); ++t_it)
    {
      IdSetup::ion_type tmp_type = IdSetup::str_to_type(*t_it);
      if(tmp_type!=IdSetup::Invalid)
      {
        selected_types_[sector].insert(tmp_type);
      }
    }
  }

  return 0;
}


//return the log score of a given prefix mass given the MS/MS Spectrum
DoubleReal BayesScoring::getLogScore(const PeakSpectrum& spectrum, const DoubleReal prefix_mass,
    std::vector<IdSetup::interpretation> &interprets)
{
  std::cout<<"enter getlogscore"<<std::endl; 
  interprets.clear();

  //the score to be returned
  DoubleReal log_score = 0;

  DoubleReal parent_mass = spectrum.getPrecursors()[0].getPosition()[0];
  DoubleReal charge = spectrum.getPrecursors()[0].getCharge();
  parent_mass = (parent_mass * charge) - (charge - 1); //TODO PROTONMASS

  //compute the sector of the prefix mass
  UInt sector = std::min(number_of_sectors_-1, (Size)floor(number_of_sectors_ * prefix_mass/parent_mass));

  //boolean vector showing whether peak at expected offset for a each ion type is abundant or absent
  std::vector<UInt> peak_levels(IdSetup::LAST_TYPE+1, 0);

  //given the mass check for all ion types for abundant peaks
  std::set<IdSetup::ion_type>::const_iterator ion_it;
  for (ion_it = selected_types_[sector].begin(); ion_it!= selected_types_[sector].end(); ++ion_it)
  {
    ion_type type_i = (ion_type) *ion_it;

    DoubleReal mz_pos = 0;

    mz_pos = IdSetup::getFragmentMass(type_i, prefix_mass, parent_mass);

    //find nearest peak to that mz_pos TODO: Maybe search for highest peak in allowed intervall as closest peak could be low level noise?
    UInt nearest_peak_id = spectrum.findNearest(mz_pos);

    //store the intensity level of the peak found
    if (fabs(mz_pos - spectrum[nearest_peak_id].getMZ()) < delta_)
    {
      interprets.push_back(IdSetup::interpretation(nearest_peak_id, type_i));
      peak_levels[type_i] = (UInt)spectrum[nearest_peak_id].getIntensity();

      //std::cout<<"mzpos "<<mz_pos<<"nearest peak id: "<<nearest_peak_id<<"  intens: "
      //<<spectrum[nearest_peak_id].getIntensity()<<"  mass: "<<spectrum[nearest_peak_id].getMZ()<<std::endl;
    }
  }
  //Set observed values in the network
  //given the detected ions compute the score
  for(ion_it = selected_types_[sector].begin(); ion_it!= selected_types_[sector].end(); ++ion_it)
  {
    networks_[sector].setOutComeIndex(IdSetup::type_to_str(*ion_it), peak_levels[*ion_it]);
  }
  networks_[sector].setOutComeIndex("class", 0);
  log_score = log(networks_[sector].getLikelihood());

  networks_[sector].setOutComeIndex("class", 1);
  log_score -= log(networks_[sector].getLikelihood());

  return log_score;
}

UInt BayesScoring::normalizeIntensity(PeakSpectrum &S, Size number_intensity_levels)
{
  const UInt max_rank=number_intensity_levels-1;
  DoubleReal parent_mass = S.getPrecursors()[0].getPosition()[0];
  DoubleReal charge = S.getPrecursors()[0].getCharge();
  parent_mass = (parent_mass * charge) - (charge - 1);

  UInt bin_size = parent_mass/100;

  S.sortByIntensity(true);

  PeakSpectrum::iterator fwit = S.begin();

  UInt count = 0;
  DoubleReal prev_inten = -10;
  for(fwit=S.begin(); fwit!=S.end(); ++fwit, ++count)
  {
    if(prev_inten == fwit->getIntensity())
    {
      fwit->setIntensity((fwit-1)->getIntensity());
    }
    else
    {
      //normalize
      prev_inten=fwit->getIntensity();
      fwit->setIntensity(std::max(0, (Int)(max_rank-(count/bin_size))));
    }
    std::cout<<"original rank: "<<count<<"  discretized:  "<<fwit->getIntensity()<<std::endl;
  }
  //bring spectrum back into normal order
  S.sortByPosition();

  return 0;
}

//normalize peak intensities as in Sherenga paper
UInt BayesScoring::normalizeIntensity(PeakSpectrum &S) const
{  
  return normalizeIntensity(S, number_intensity_levels_);
}


//return the log score of a given prefix mass given the MS/MS Spectrum
DoubleReal BayesScoring::getSupportIonCount(const PeakSpectrum& spectrum, const DoubleReal prefix_mass) const
{
  DoubleReal supporters=0.0;  

  DoubleReal parent_mass = spectrum.getPrecursors()[0].getPosition()[0];
  DoubleReal charge = spectrum.getPrecursors()[0].getCharge();
  parent_mass = (parent_mass * charge) - (charge - 1);  

  //given the mass check for all ion types for abundant peaks
  std::set<ion_type>::const_iterator ion_it;
  for (ion_it = selected_types_[0].begin(); ion_it!= selected_types_[0].end(); ++ion_it)
  {
    ion_type type_i = (ion_type) *ion_it;

    DoubleReal mz_pos = IdSetup::getFragmentMass(type_i, prefix_mass, parent_mass);

    //find nearest peak to that mz_pos
    UInt nearest_peak_id = spectrum.findNearest(mz_pos);

    if(fabs(mz_pos - spectrum[nearest_peak_id].getMZ()) <= delta_ && spectrum[nearest_peak_id].getIntensity() !=0)
    {
      //count the monoisotopic peak      
      supporters+=spectrum[nearest_peak_id].getIntensity();

      //search for isotopic peak (+1, +0.5 for doubly charged)
      DoubleReal isotope_offset=1;
      if(type_i==IdSetup::BIon2 || type_i==IdSetup::YIon2)
      {
        isotope_offset=0.5;
      }
      Size iso_peak_id = spectrum.findNearest(spectrum[nearest_peak_id].getMZ()+isotope_offset);
      if(fabs((spectrum[iso_peak_id].getMZ()-spectrum[nearest_peak_id].getMZ())-isotope_offset) <=0.1)
      {
        if(spectrum[iso_peak_id].getIntensity()<spectrum[nearest_peak_id].getIntensity())
        {          
          supporters+=spectrum[iso_peak_id].getIntensity();
        }
      }            
    }
  }
  return supporters;
}




}//end of namespace



