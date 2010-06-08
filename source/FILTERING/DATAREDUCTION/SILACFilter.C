// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <cmath>

namespace OpenMS
{



SILACFilter::SILACFilter(DoubleReal mass_separation_light_heavy, Int charge_) {
	silac_type=DOUBLE;
	envelope_distance_light_medium=0;
	charge=charge_;
	envelope_distance_light_heavy=mass_separation_light_heavy/(DoubleReal)charge;
	isotope_distance=1.0/(DoubleReal) charge;
}

SILACFilter::SILACFilter(DoubleReal mass_separation_light_medium,DoubleReal mass_separation_light_heavy,Int charge_) {
	silac_type=TRIPLE;
	charge=charge_;
	envelope_distance_light_medium=mass_separation_light_medium/(DoubleReal)charge;
	envelope_distance_light_heavy=mass_separation_light_heavy/(DoubleReal)charge;
	isotope_distance=1.0/(DoubleReal) charge;
}

SILACFilter::~SILACFilter() {
}

bool SILACFilter::blacklisted(DoubleReal value)
{
	std::set<DoubleReal,doubleCmp>::iterator value_pos=blacklist.find(value);
	return value_pos!=blacklist.end();
}


void SILACFilter::blockValue(DoubleReal value)
{
	blacklist.insert(value);
}


void SILACFilter::reset()
{
	blacklist.clear();
}

bool SILACFilter::isFeature(DoubleReal act_rt,DoubleReal act_mz)
{
	if (!blacklisted(act_mz) /*&& !blacklisted(act_mz+isotope_distance) && !blacklisted(act_mz+2*isotope_distance) && !blacklisted(act_mz+envelope_distance_light_heavy) && !blacklisted(act_mz+envelope_distance_light_heavy+isotope_distance) && !blacklisted(act_mz+envelope_distance_light_heavy+2*isotope_distance)*/)
	{
		std::pair<bool,std::vector<DoubleReal> > val=getIntensities(act_mz,SILACFiltering::acc_spl,SILACFiltering::spline_spl);
		if (!val.first)
			return false;
		std::vector<DoubleReal> intensities_spl =val.second;
		std::vector<DoubleReal> intensities_lin = getIntensities(act_mz,SILACFiltering::acc_lin,SILACFiltering::spline_lin).second;
		bool condDouble1=false;
		bool condDouble2=false;
		bool condDouble3=false;
		bool condTriple1=false;
		bool condTriple2=false;
		bool condTriple3=false;
		if (silac_type==DOUBLE)
		{
			condDouble1 = (intensities_lin[1] >= SILACFiltering::intensity_cutoff) && (intensities_lin[2] >= SILACFiltering::intensity_cutoff) && (intensities_lin[3] >= SILACFiltering::intensity_cutoff) && (intensities_lin[5] >= SILACFiltering::intensity_cutoff) && (intensities_lin[6] >= SILACFiltering::intensity_cutoff) && (intensities_lin[7] >= SILACFiltering::intensity_cutoff); // all six intensities peak simultaneously
			condDouble2 = (intensities_spl[1] >= intensities_spl[2]) && (intensities_spl[2] >= intensities_spl[3]) && (intensities_spl[5] >= intensities_spl[6]) && (intensities_spl[6] >= intensities_spl[7]); // isotopic peaks within one envelop decrease
			condDouble3 = (intensities_spl[0]<=intensities_spl[1] && intensities_spl[4]<= intensities_spl[5]) || (intensities_spl[0]>=intensities_spl[1] && intensities_spl[4]>= intensities_spl[5]); //(std::abs((intensities_lin[1]/intensities_lin[2])-(intensities_lin[5]/intensities_lin[6])) < 0.5 /*&& std::abs((intensities_spl[2]/intensities_spl[3])-(intensities_spl[6]/intensities_spl[7])) < 1*/);
		}
		if (silac_type==TRIPLE)
		{
			condTriple1 = (intensities_lin[1] >= SILACFiltering::intensity_cutoff) && (intensities_lin[2] >= SILACFiltering::intensity_cutoff) && (intensities_lin[3] >= SILACFiltering::intensity_cutoff) && (intensities_lin[5] >= SILACFiltering::intensity_cutoff) && (intensities_lin[6] >= SILACFiltering::intensity_cutoff) && (intensities_lin[7] >= SILACFiltering::intensity_cutoff) && (intensities_lin[9] >= SILACFiltering::intensity_cutoff) && (intensities_lin[10] >= SILACFiltering::intensity_cutoff) && (intensities_lin[11] >= SILACFiltering::intensity_cutoff); // all nine intensities peak simultaneously
			condTriple2 = (intensities_lin[0] >= intensities_lin[1]) && (intensities_lin[1] >= intensities_lin[2]) && (intensities_lin[3] >= intensities_lin[4]) && (intensities_lin[4] >= intensities_lin[5]) && (intensities_lin[6] >= intensities_lin[7]) && (intensities_lin[7] >= intensities_lin[8]); // isotopic peaks within one envelop decrease
			condTriple3 = ((intensities_spl[0]<=intensities_spl[1] && intensities_spl[4]<= intensities_spl[5] && intensities_spl[8]<= intensities_spl[9]) || (intensities_spl[0]>=intensities_spl[1] && intensities_spl[4]>= intensities_spl[5] && intensities_spl[8]>= intensities_spl[9]));
		}
		if (condTriple1 && condTriple2 && condTriple3 && silac_type==TRIPLE)
		{
			DataPoint next_element;
			next_element.feature_id=SILACFiltering::feature_id;
			next_element.rt=act_rt;
			next_element.mz=act_mz;
			next_element.silac_type=TRIPLE;
			next_element.charge=charge;
			next_element.envelope_distance_light_heavy=envelope_distance_light_heavy;
			next_element.intensities.push_back(intensities_spl[1]);
			next_element.intensities.push_back(intensities_spl[2]);
			next_element.intensities.push_back(intensities_spl[3]);
			next_element.intensities.push_back(intensities_spl[5]);
			next_element.intensities.push_back(intensities_spl[6]);
			next_element.intensities.push_back(intensities_spl[7]);
			next_element.intensities.push_back(intensities_spl[9]);
			next_element.intensities.push_back(intensities_spl[10]);
			next_element.intensities.push_back(intensities_spl[11]);
			elements.push_back(next_element);
			peak_values.clear();
			peak_values.push_back(act_mz);
			peak_values.push_back(act_mz+isotope_distance);
			peak_values.push_back(act_mz+2*isotope_distance);
			peak_values.push_back(act_mz+envelope_distance_light_medium);
			peak_values.push_back(act_mz+envelope_distance_light_medium+isotope_distance);
			peak_values.push_back(act_mz+envelope_distance_light_medium+2*isotope_distance);
			peak_values.push_back(act_mz+envelope_distance_light_heavy);
			peak_values.push_back(act_mz+envelope_distance_light_heavy+isotope_distance);
			peak_values.push_back(act_mz+envelope_distance_light_heavy+2*isotope_distance);
			return true;
		}
		if (condDouble1 && condDouble2 && condDouble3 && silac_type==DOUBLE)
		{
//			if (SILACFiltering::feature_id==1512)
//			{
//				std::vout << intensities_spl[1] << " " << intensities_spl[2] << " " << intensities_spl[3] << " "<< intensities_spl[] << " "
//			}
			DataPoint next_element;
			next_element.feature_id=SILACFiltering::feature_id;
			next_element.rt=act_rt;
			next_element.mz=act_mz;
			next_element.silac_type=DOUBLE;
			next_element.charge=charge;
			next_element.envelope_distance_light_heavy=envelope_distance_light_heavy;
			next_element.intensities.push_back(intensities_spl[1]);
			next_element.intensities.push_back(intensities_spl[2]);
			next_element.intensities.push_back(intensities_spl[3]);
			next_element.intensities.push_back(intensities_spl[5]);
			next_element.intensities.push_back(intensities_spl[6]);
			next_element.intensities.push_back(intensities_spl[7]);
			elements.push_back(next_element);
			peak_values.clear();
			peak_values.push_back(act_mz);
			peak_values.push_back(act_mz+isotope_distance);
			peak_values.push_back(act_mz+2*isotope_distance);
			peak_values.push_back(act_mz+envelope_distance_light_heavy);
			peak_values.push_back(act_mz+envelope_distance_light_heavy+isotope_distance);
			peak_values.push_back(act_mz+envelope_distance_light_heavy+2*isotope_distance);
			return true;
		}
	}
//	else
//	{
//		for (std::set<DoubleReal>::iterator it=blacklist.begin();it!=blacklist.end();++it)
//		{
//			std::cout << *it << " ";
//		}
//		std::cout << std::endl << act_mz << " | " << *blacklist.find(act_mz) << std::endl;
//	}
	return false;
}


std::pair<bool,std::vector<DoubleReal> > SILACFilter::getIntensities(DoubleReal act_mz, gsl_interp_accel* acc,gsl_spline* spline)
{
	std::vector<DoubleReal> splines;
	double isMaximum=true;
	splines.push_back(gsl_spline_eval (spline, act_mz-isotope_distance, acc));
	DoubleReal value1=gsl_spline_eval (spline, act_mz-0.01, acc);
	DoubleReal value=gsl_spline_eval (spline, act_mz, acc);
	DoubleReal value2=gsl_spline_eval (spline, act_mz+0.01, acc);
	if (value1>value || value2>value)
		isMaximum=false;
	splines.push_back(value);
	value1=gsl_spline_eval (spline, act_mz+isotope_distance-0.01, acc);
	value=gsl_spline_eval (spline, act_mz+isotope_distance, acc);
	value2=gsl_spline_eval (spline, act_mz+isotope_distance+0.01, acc);
	if (value1>value || value2>value)
		isMaximum=false;
	splines.push_back(value);
	value1=gsl_spline_eval (spline, act_mz+2*isotope_distance-0.01, acc);
	value=gsl_spline_eval (spline, act_mz+2*isotope_distance, acc);
	value2=gsl_spline_eval (spline, act_mz+2*isotope_distance+0.01, acc);
	if (value1>value || value2>value)
		isMaximum=false;
	splines.push_back(value);
	if (silac_type==TRIPLE)
	{
		splines.push_back(gsl_spline_eval (spline, act_mz-isotope_distance, acc));
		splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_medium, acc));
		splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_medium+isotope_distance, acc));
		splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_medium+0.01, acc));
	}
	splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy-isotope_distance, acc));
	value1=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy-0.01, acc);
	value=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy, acc);
	value2=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+0.01, acc);
	if (value1>value || value2>value)
		isMaximum=false;
	splines.push_back(value);
	value1=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+isotope_distance-0.01, acc);
	value=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+isotope_distance, acc);
	value2=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+isotope_distance+0.01, acc);
	if (value1>value || value2>value)
		isMaximum=false;
	splines.push_back(value);
	value1=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+2*isotope_distance-0.01, acc);
	value=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+2*isotope_distance, acc);
	value2=gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+2*isotope_distance+0.01, acc);
	if (value1>value || value2>value)
		isMaximum=false;
	splines.push_back(value);
	return std::make_pair(isMaximum,splines);
}

bool SILACFilter::doubleCmp::operator()(DoubleReal a, DoubleReal b) const
{
	if ( fabs(a-b) < 0.022/*2*SILACFiltering::mz_stepwidth*/)
	{
		return false;
	}
	else
	{
		return a<b;
	}
}

Int SILACFilter::getSILACType()
{
	return silac_type;
}

std::vector<DoubleReal> SILACFilter::getPeakValues()
{
	return peak_values;
}

DoubleReal SILACFilter::getEnvelopeDistanceLightMedium()
{
	return envelope_distance_light_medium;
}
DoubleReal SILACFilter::getEnvelopeDistanceLightHeavy()
{
	return envelope_distance_light_heavy;
}

DoubleReal SILACFilter::getIsotopeDistance()
{
	return isotope_distance;
}

std::vector<DataPoint> SILACFilter::getElements()
{
	return elements;
}

}
