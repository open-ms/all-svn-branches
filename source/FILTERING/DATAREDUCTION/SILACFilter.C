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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <cmath>

namespace OpenMS
{



SILACFilter::SILACFilter(DoubleReal mass_separation_light_heavy, Int charge_,DoubleReal model_deviation_) {
	silac_type=DOUBLE;
	envelope_distance_light_medium=0;
	charge=charge_;
	envelope_distance_light_heavy=mass_separation_light_heavy/(DoubleReal)charge;
	isotope_distance=1.000495/(DoubleReal) charge;
	model_deviation=model_deviation_;
}

SILACFilter::SILACFilter(DoubleReal mass_separation_light_medium,DoubleReal mass_separation_light_heavy,Int charge_,DoubleReal model_deviation_) {
	silac_type=TRIPLE;
	charge=charge_;
	envelope_distance_light_medium=mass_separation_light_medium/(DoubleReal)charge;
	envelope_distance_light_heavy=mass_separation_light_heavy/(DoubleReal)charge;
	isotope_distance=1.000495/(DoubleReal) charge;
	model_deviation=model_deviation_;
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

DoubleReal mean(std::vector<DoubleReal>& values)
{
	DoubleReal mean= 0.0;
	for ( Size i = 0; i < values.size(); i++ )
	{
		mean+=values[i];
	}
	return mean/values.size();
}

DoubleReal standardDeviation(std::vector<DoubleReal>& values)
{
	DoubleReal sd = 0.0;
	DoubleReal mv=mean(values);

	for( Size i = 0; i < values.size(); i++ )
	{
		sd += pow((values[i] - mv),2);

	}
	sd /= values.size();
	return std::sqrt(sd);
}

bool SILACFilter::isFeature(DoubleReal act_rt,DoubleReal act_mz)
{
	if (!blacklisted(act_mz) &&  !blacklisted(act_mz+isotope_distance) && !blacklisted(act_mz+2*isotope_distance) && !blacklisted(act_mz+envelope_distance_light_heavy) && !blacklisted(act_mz+envelope_distance_light_heavy+isotope_distance) && !blacklisted(act_mz+envelope_distance_light_heavy+2*isotope_distance))
	{
		std::vector<DoubleReal> intensities_spl = getIntensities(act_mz,SILACFiltering::acc_spl,SILACFiltering::spline_spl);
		std::vector<DoubleReal> intensities_lin = getIntensities(act_mz,SILACFiltering::acc_lin,SILACFiltering::spline_lin);
		bool condDouble1=false;
		bool condDouble3=false;
		bool condTriple1=false;
		bool condTriple2=false;
		bool condTriple3=false;
		if (silac_type==DOUBLE)
		{
			condDouble1 = ((intensities_lin[1] >= SILACFiltering::intensity_cutoff) && (intensities_lin[2] >= SILACFiltering::intensity_cutoff) && (intensities_lin[3] >= SILACFiltering::intensity_cutoff) && (intensities_lin[5] >= SILACFiltering::intensity_cutoff) && (intensities_lin[6] >= SILACFiltering::intensity_cutoff)) || ((intensities_lin[1] >= SILACFiltering::intensity_cutoff) && (intensities_lin[2] >= SILACFiltering::intensity_cutoff) && (intensities_lin[5] >= SILACFiltering::intensity_cutoff) && (intensities_lin[6] >= SILACFiltering::intensity_cutoff) && (intensities_lin[7] >= SILACFiltering::intensity_cutoff)); // all six intensities peak simultaneously
			condDouble3 = (intensities_spl[0]<=intensities_spl[1] && intensities_spl[4]<= intensities_spl[5]) || (intensities_spl[0]>=intensities_spl[1] && intensities_spl[4]>= intensities_spl[5]);
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
		if (condDouble1/* && condDouble2 */&& condDouble3 && silac_type==DOUBLE)
		{
			DoubleReal light_correlation1=getPeakCorrelation(act_mz);
			DoubleReal light_correlation2=0.85;
			if (intensities_lin[3]>0.0)
				light_correlation2=getPeakCorrelation(act_mz+isotope_distance);
			DoubleReal heavy_correlation1=getPeakCorrelation(act_mz+envelope_distance_light_heavy);
			DoubleReal heavy_correlation2=0.85;
			if (intensities_lin[7]>0.0)
				heavy_correlation2=getPeakCorrelation(act_mz+envelope_distance_light_heavy+isotope_distance);
			DoubleReal correlation=light_correlation1+light_correlation2+heavy_correlation1+heavy_correlation2;
			correlation/=4;

			if (correlation < 0.85)
				return false;

			ExtendedIsotopeModel model_light;
			Param tmp;
			tmp.setValue( "isotope:monoisotopic_mz", act_mz );
			tmp.setValue( "interpolation_step", 0.01 );
			tmp.setValue("charge",charge);
			tmp.setValue("isotope:stdev",getPeakWidth(act_mz)/2);
			model_light.setParameters( tmp );
			std::vector<Peak1D> model_data;
			model_light.getSamples(model_data);
			if (intensities_lin[3]==0.0)
				intensities_spl[3]=(1-model_deviation)*(intensities_spl[2]/model_light.getIntensity(act_mz+isotope_distance));

			DoubleReal quality_l1=(intensities_spl[1]*model_light.getIntensity(act_mz+isotope_distance))/(intensities_spl[2]*model_light.getIntensity(act_mz));
			DoubleReal quality_l2=(intensities_spl[2]*model_light.getIntensity(act_mz+2*isotope_distance))/(intensities_spl[3]*model_light.getIntensity(act_mz+isotope_distance));

			ExtendedIsotopeModel model_heavy;
			Param tmp1;
			tmp1.setValue( "isotope:monoisotopic_mz", act_mz+envelope_distance_light_heavy );
			tmp1.setValue( "interpolation_step", 0.01 );
			tmp1.setValue("charge",charge);
			tmp1.setValue("isotope:stdev",getPeakWidth(act_mz+envelope_distance_light_heavy)/2);
			model_heavy.setParameters( tmp1 );

			if (intensities_lin[7]==0.0)
				intensities_spl[7]=(1-model_deviation)*(intensities_spl[6]/model_heavy.getIntensity(act_mz+envelope_distance_light_heavy+isotope_distance));

			DoubleReal quality_h1=(intensities_spl[5]*model_heavy.getIntensity(act_mz+envelope_distance_light_heavy+isotope_distance))/(intensities_spl[6]*model_heavy.getIntensity(act_mz+envelope_distance_light_heavy));
			DoubleReal quality_h2=(intensities_spl[6]*model_heavy.getIntensity(act_mz+envelope_distance_light_heavy+2*isotope_distance))/(intensities_spl[7]*model_heavy.getIntensity(act_mz+envelope_distance_light_heavy+isotope_distance));

			//False negative debug output
//			if (act_rt < 2627.2 && act_rt > 2627.0 && act_mz > 426.0 && act_mz < 430.0)
//			{
//				std::cout << std::endl << "mz: " << act_mz << " rt: " << act_rt << std::endl;
//				std::cout << "light structure/model: " << intensities_spl[1] << "/" << model_light.getIntensity(act_mz)  << " " << intensities_spl[2] << "/" << model_light.getIntensity(act_mz+isotope_distance) << " " << intensities_spl[3] << "/" << model_light.getIntensity(act_mz+2*isotope_distance) << std::endl;
//				std::cout << "heavy structure/model: " << intensities_spl[5] << "/" << model_heavy.getIntensity(act_mz+envelope_distance_light_heavy)  << " " << intensities_spl[6] << "/" << model_heavy.getIntensity(act_mz+isotope_distance+envelope_distance_light_heavy) << " " << intensities_spl[7] << "/" << model_heavy.getIntensity(act_mz+2*isotope_distance+envelope_distance_light_heavy) << std::endl;
//				std::cout << "intensity values (light value vs. model value):" << std::endl;
//				for (DoubleReal position=act_mz-0.1;position<=act_mz+2*isotope_distance+0.1;position+=0.01)
//				{
//					DoubleReal intensity=gsl_spline_eval (SILACFiltering::spline_spl, position, SILACFiltering::acc_spl);
//					std::cout << intensity << "\t" << model_light.getIntensity(position) << std::endl;
//				}
//				std::cout << std::endl;
//			}

			if (std::abs(quality_l1-1) > model_deviation || std::abs(quality_l2-1) > model_deviation)
				return false;

			if (std::abs(quality_h1-1) > model_deviation ||std::abs(quality_h2-1) > model_deviation)
				return false;

			//False positive debug output
//			if (SILACFiltering::feature_id == 140)
//			{
//				std::cout << std::endl << "mz: " << act_mz << " rt: " << act_rt << "\t" << act_mz+envelope_distance_light_heavy << std::endl;
//				std::cout << "light structure/model: " << "("  << intensities_spl[0] << ")" << intensities_spl[1] << "/" << model_light.getIntensity(act_mz)  << " " << intensities_spl[2] << "/" << model_light.getIntensity(act_mz+isotope_distance) << " " << intensities_spl[3] << "/" << model_light.getIntensity(act_mz+2*isotope_distance) << std::endl;
//				std::cout << "heavy structure/model: "<< "("  << intensities_spl[4] << ")" << intensities_spl[5] << "/" << model_heavy.getIntensity(act_mz+envelope_distance_light_heavy)  << " " << intensities_spl[6] << "/" << model_heavy.getIntensity(act_mz+isotope_distance+envelope_distance_light_heavy) << " " << intensities_spl[7] << "/" << model_heavy.getIntensity(act_mz+2*isotope_distance+envelope_distance_light_heavy) << std::endl;
//				std::cout << "intensity values (light value vs. model value):" << std::endl;
//				for (DoubleReal position=act_mz+envelope_distance_light_heavy-0.1;position<=act_mz+envelope_distance_light_heavy+2*isotope_distance+0.1;position+=0.0001)
//				{
//					DoubleReal intensity=gsl_spline_eval (SILACFiltering::spline_spl, position, SILACFiltering::acc_spl);
//					std::cout << intensity << "\t" << model_light.getIntensity(position) << std::endl;
//				}
//				std::cout << std::endl;
//			}

			DoubleReal quality= correlation;
			DataPoint next_element;
			next_element.feature_id=SILACFiltering::feature_id;
			next_element.rt=act_rt;
			next_element.mz=act_mz;
			next_element.silac_type=DOUBLE;
			next_element.charge=charge;
			next_element.envelope_distance_light_heavy=envelope_distance_light_heavy;
			next_element.quality=quality;
			next_element.intensities.clear();
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
	return false;
}

DoubleReal SILACFilter::getPeakCorrelation(DoubleReal act_mz)
{
	std::vector<DoubleReal> first_values;
	std::vector<DoubleReal> second_values;
	for (DoubleReal pos=act_mz-2*SILACFiltering::mz_stepwidth;pos<=act_mz+2*SILACFiltering::mz_stepwidth;pos+=SILACFiltering::mz_stepwidth)
	{
		DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl);
		DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+isotope_distance, SILACFiltering::acc_spl);
		first_values.push_back(intensity1);
		second_values.push_back(intensity2);
	}
	return Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
}


std::vector<DoubleReal> SILACFilter::getIntensities(DoubleReal act_mz, gsl_interp_accel* acc,gsl_spline* spline)
{
	std::vector<DoubleReal> splines;
	splines.push_back(gsl_spline_eval (spline, act_mz-isotope_distance, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+isotope_distance, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+2*isotope_distance, acc));
	if (silac_type==TRIPLE)
	{
		splines.push_back(gsl_spline_eval (spline, act_mz-isotope_distance, acc));
		splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_medium, acc));
		splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_medium+isotope_distance, acc));
		splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_medium+0.01, acc));
	}
	splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy-isotope_distance, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+isotope_distance, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+envelope_distance_light_heavy+2*isotope_distance, acc));
	return splines;
}

bool SILACFilter::doubleCmp::operator()(DoubleReal a, DoubleReal b) const
{
	DoubleReal peak_width=SILACFilter::getPeakWidth((a+b)/2);
	if ( fabs(a-b) < peak_width/*2*SILACFiltering::mz_stepwidth*/)
	{
		return false;
	}
	else
	{
		return a<b;
	}
}

DoubleReal SILACFilter::getPeakWidth(DoubleReal mz)
{
	return 11*(1.828e-7*pow(mz,1.504));
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
