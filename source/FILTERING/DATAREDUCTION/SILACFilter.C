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



SILACFilter::SILACFilter(std::vector<DoubleReal> mass_separations, Int charge_,DoubleReal model_deviation_) {
	silac_type=mass_separations.size();
	charge=charge_;
	for (std::vector<DoubleReal>::iterator it=mass_separations.begin();it!=mass_separations.end();++it)
	{
		envelope_distances.insert(*it/(DoubleReal)charge);
	}
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
	if (blacklisted(act_mz) ||  blacklisted(act_mz+isotope_distance) || blacklisted(act_mz+2*isotope_distance))
		return false;
	std::vector<DoubleReal> intensities;
	std::vector<DoubleReal> intensities_spl_light = getIntensities(act_mz,SILACFiltering::acc_spl,SILACFiltering::spline_spl);
	std::vector<DoubleReal> intensities_lin_light = getIntensities(act_mz,SILACFiltering::acc_lin,SILACFiltering::spline_lin);
	intensities.push_back(intensities_spl_light[1]);
	intensities.push_back(intensities_spl_light[2]);
	intensities.push_back(intensities_spl_light[3]);
	std::vector<DoubleReal> act_peak_positions;
	act_peak_positions.push_back(act_mz);
	act_peak_positions.push_back(act_mz+isotope_distance);
	act_peak_positions.push_back(act_mz+2*isotope_distance);
	for (std::set<DoubleReal>::iterator envelope_iterator=envelope_distances.begin();envelope_iterator!=envelope_distances.end();++envelope_iterator)
	{
		DoubleReal envelope_distance=*envelope_iterator;
		 if ( blacklisted(act_mz+envelope_distance) ||  blacklisted(act_mz+envelope_distance+isotope_distance) || blacklisted(act_mz+envelope_distance+2*isotope_distance))
				 return false;
		std::vector<DoubleReal> intensities_spl_heavy = getIntensities(act_mz+envelope_distance,SILACFiltering::acc_spl,SILACFiltering::spline_spl);
		std::vector<DoubleReal> intensities_lin_heavy = getIntensities(act_mz+envelope_distance,SILACFiltering::acc_lin,SILACFiltering::spline_lin);

		std::vector<DoubleReal> intensities_spl(intensities_spl_light.begin(),intensities_spl_light.end());
		intensities_spl.insert(intensities_spl.end(),intensities_spl_heavy.begin(),intensities_spl_heavy.end());
		std::vector<DoubleReal> intensities_lin(intensities_lin_light.begin(),intensities_lin_light.end());
		intensities_lin.insert(intensities_lin.end(),intensities_lin_heavy.begin(),intensities_lin_heavy.end());

		if (checkArea(act_mz, envelope_distance,intensities_spl, intensities_lin))
		{
			intensities.push_back(intensities_spl_heavy[1]);
			intensities.push_back(intensities_spl_heavy[2]);
			intensities.push_back(intensities_spl_heavy[3]);
			act_peak_positions.push_back(act_mz+envelope_distance);
			act_peak_positions.push_back(act_mz+envelope_distance+isotope_distance);
			act_peak_positions.push_back(act_mz+envelope_distance+2*isotope_distance);
		}
		else
		{
			return false;
		}
	}
	DataPoint next_element;
	next_element.feature_id=SILACFiltering::feature_id;
	next_element.rt=act_rt;
	next_element.mz=act_mz;
	next_element.silac_type=envelope_distances.size();
	next_element.charge=charge;
	next_element.envelope_distance_light_heavy=*envelope_distances.rbegin();
	next_element.intensities.insert(next_element.intensities.end(),intensities.begin(),intensities.end());
	elements.push_back(next_element);
	peak_positions.clear();
	peak_positions.insert(peak_positions.begin(),act_peak_positions.begin(),act_peak_positions.end());
	return true;

}

bool SILACFilter::checkArea(DoubleReal act_mz, DoubleReal envelope_distance,std::vector<DoubleReal>& intensities_spl, std::vector<DoubleReal>& intensities_lin)
{
		bool condition1 = ((intensities_lin[1] >= SILACFiltering::intensity_cutoff) && (intensities_lin[2] >= SILACFiltering::intensity_cutoff) && (intensities_lin[3] >= SILACFiltering::intensity_cutoff) && (intensities_lin[5] >= SILACFiltering::intensity_cutoff) && (intensities_lin[6] >= SILACFiltering::intensity_cutoff)  && (intensities_lin[7] >= SILACFiltering::intensity_cutoff)); // all six intensities peak simultaneously
		bool condition2 = (intensities_spl[0]<=intensities_spl[1] && intensities_spl[4]<= intensities_spl[5]) || (intensities_spl[0]>=intensities_spl[1] && intensities_spl[4]>= intensities_spl[5]);
		if (condition1 && condition2)
		{
			DoubleReal light_correlation1=getPeakCorrelation(act_mz);
			DoubleReal light_correlation2=getPeakCorrelation(act_mz+isotope_distance);

			DoubleReal heavy_correlation1=getPeakCorrelation(act_mz+envelope_distance);
			DoubleReal heavy_correlation2=getPeakCorrelation(act_mz+envelope_distance+isotope_distance);

			DoubleReal correlation=light_correlation1+light_correlation2+heavy_correlation1+heavy_correlation2;
			correlation/=4;
//			std::cout << light_correlation1 << " " << light_correlation2 << " " << heavy_correlation1 << " " << heavy_correlation1 << std::endl;

			if (light_correlation1 < 0.9 || light_correlation2 < 0.9 || heavy_correlation1 < 0.9 || heavy_correlation2 < 0.9)
				return false;

			ExtendedIsotopeModel model_light;
			Param param_light;
			param_light.setValue( "isotope:monoisotopic_mz", act_mz );
			param_light.setValue( "interpolation_step", 0.01 );
			param_light.setValue("charge",charge);
			param_light.setValue("isotope:stdev",getPeakWidth(act_mz)/2);
			model_light.setParameters( param_light );
			std::vector<Peak1D> model_data;
			model_light.getSamples(model_data);

			DoubleReal quality_l1=(intensities_spl[1]*model_light.getIntensity(act_mz+isotope_distance))/(intensities_spl[2]*model_light.getIntensity(act_mz));
			DoubleReal quality_l2=(intensities_spl[2]*model_light.getIntensity(act_mz+2*isotope_distance))/(intensities_spl[3]*model_light.getIntensity(act_mz+isotope_distance));

			ExtendedIsotopeModel model_heavy;
			Param param_heavy;
			param_heavy.setValue( "isotope:monoisotopic_mz", act_mz+envelope_distance);
			param_heavy.setValue( "interpolation_step", 0.01 );
			param_heavy.setValue("charge",charge);
			param_heavy.setValue("isotope:stdev",getPeakWidth(act_mz+envelope_distance)/2);
			model_heavy.setParameters( param_heavy );

			DoubleReal quality_h1=(intensities_spl[5]*model_heavy.getIntensity(act_mz+envelope_distance+isotope_distance))/(intensities_spl[6]*model_heavy.getIntensity(act_mz+envelope_distance));
			DoubleReal quality_h2=(intensities_spl[6]*model_heavy.getIntensity(act_mz+envelope_distance+2*isotope_distance))/(intensities_spl[7]*model_heavy.getIntensity(act_mz+envelope_distance+isotope_distance));

			//False negative debug output
//			if (act_rt < 1665.55 && act_rt > 1665.45 && act_mz > 627.0 && act_mz < 632.0)
//			{
//				std::cout << std::endl << "mz: " << act_mz << std::endl;
//				std::cout << "light structure/model: " << intensities_spl[1] << "/" << model_light.getIntensity(act_mz)  << " " << intensities_spl[2] << "/" << model_light.getIntensity(act_mz+isotope_distance) << " " << intensities_spl[3] << "/" << model_light.getIntensity(act_mz+2*isotope_distance) << std::endl;
//				std::cout << "heavy structure/model: " << intensities_spl[5] << "/" << model_heavy.getIntensity(act_mz+envelope_distance)  << " " << intensities_spl[6] << "/" << model_heavy.getIntensity(act_mz+isotope_distance+envelope_distance) << " " << intensities_spl[7] << "/" << model_heavy.getIntensity(act_mz+2*isotope_distance+envelope_distance) << std::endl;
//				std::cout << "intensity values (light value vs. model value):" << std::endl;
//				std::cout << "Correlation: " << light_correlation1  << " " << light_correlation2  << " " << heavy_correlation1  << " " << heavy_correlation2 << std::endl;
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

//			False positive debug output
//			if (SILACFiltering::feature_id == 305)
//			{
//				std::cout << std::endl << "mz: " << act_mz <<  "\t" << act_mz+envelope_distance << std::endl;
//				std::cout << "light structure/model: " << "("  << intensities_spl[0] << ")" << intensities_spl[1] << "/" << model_light.getIntensity(act_mz)  << " " << intensities_spl[2] << "/" << model_light.getIntensity(act_mz+isotope_distance) << " " << intensities_spl[3] << "/" << model_light.getIntensity(act_mz+2*isotope_distance) << std::endl;
//				std::cout << "heavy structure/model: "<< "("  << intensities_spl[4] << ")" << intensities_spl[5] << "/" << model_heavy.getIntensity(act_mz+envelope_distance)  << " " << intensities_spl[6] << "/" << model_heavy.getIntensity(act_mz+isotope_distance+envelope_distance) << " " << intensities_spl[7] << "/" << model_heavy.getIntensity(act_mz+2*isotope_distance+envelope_distance) << std::endl;
//				std::cout << "intensity values (light value vs. model value):" << std::endl;
//				for (DoubleReal position=act_mz+envelope_distance-0.1;position<=act_mz+envelope_distance+2*isotope_distance+0.1;position+=0.0001)
//				{
//					DoubleReal intensity=gsl_spline_eval (SILACFiltering::spline_spl, position, SILACFiltering::acc_spl);
//					std::cout << intensity << "\t" << model_light.getIntensity(position) << std::endl;
//				}
//				std::cout << std::endl;
//			}
			return true;
		}
	return false;

}



DoubleReal SILACFilter::getPeakCorrelation(DoubleReal act_mz)
{
	std::vector<DoubleReal> first_values;
	std::vector<DoubleReal> second_values;
//	if (SILACFiltering::feature_id==1401)
//	std::cout << SILACFiltering::feature_id << " " << act_mz << std::endl;
		for (DoubleReal pos=act_mz-0.04;pos<=act_mz+0.04;pos+=0.01)
		{
			DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl);
			DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+isotope_distance, SILACFiltering::acc_spl);
			first_values.push_back(intensity1);
			second_values.push_back(intensity2);
//	if (SILACFiltering::feature_id==1401)
//	std::cout << intensity1 << "\t" << intensity2 << std::endl;
		}
//	if (SILACFiltering::feature_id==1401)
//	std::cout << std::endl << std::endl;
//		for (Int i=0;i<first_values.size();++i)
//		{
//			std::cout << first_values[i] << " ";
//		}
//		std::cout << std::endl;
//		for (Int i=0;i<first_values.size();++i)
//				{
//					std::cout << second_values[i] << " ";
//				}
//		std::cout << std::endl;
	return Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
}


std::vector<DoubleReal> SILACFilter::getIntensities(DoubleReal act_mz, gsl_interp_accel* acc,gsl_spline* spline)
{
	std::vector<DoubleReal> splines;
	splines.push_back(gsl_spline_eval (spline, act_mz-isotope_distance, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+isotope_distance, acc));
	splines.push_back(gsl_spline_eval (spline, act_mz+2*isotope_distance, acc));
	return splines;
}

bool SILACFilter::doubleCmp::operator()(DoubleReal a, DoubleReal b) const
{
	DoubleReal peak_width=getPeakWidth((a+b)/2);
	if ( fabs(a-b) < peak_width)
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
	return 3*(1.828e-7*pow(mz,1.504));
}

Int SILACFilter::getSILACType()
{
	return envelope_distances.size();
}

std::vector<DoubleReal> SILACFilter::getPeakPositions()
{
	return peak_positions;
}

std::set<DoubleReal> SILACFilter::getEnvelopeDistances()
{
	return envelope_distances;
}

DoubleReal SILACFilter::getIsotopeDistance()
{
	return isotope_distance;
}

std::vector<DataPoint> SILACFilter::getElements()
{
	return elements;
}

Int SILACFilter::getCharge()
{
	return charge;
}

}
