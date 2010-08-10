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

bool debug=false;

SILACFilter::SILACFilter(std::set<DoubleReal> mass_separations, Int charge_,DoubleReal model_deviation_) {
	silac_type=mass_separations.size();
	charge=charge_;
	for (std::set<DoubleReal>::iterator it=mass_separations.begin();it!=mass_separations.end();++it)
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

bool SILACFilter::isFeature(DoubleReal act_rt,DoubleReal act_mz)
{
	if (act_rt > 1601 && act_rt < 1602)
		debug=true;
	else
		debug=false;
	if (blacklisted(act_mz) ||  blacklisted(act_mz+isotope_distance) || blacklisted(act_mz+2*isotope_distance))
		return false;
	std::vector<DoubleReal> intensities;
	std::vector<DoubleReal> act_peak_positions;
	std::vector<DoubleReal> intensities_spl_light = getIntensities(act_mz,SILACFiltering::acc_spl,SILACFiltering::spline_spl);
	std::vector<DoubleReal> intensities_lin_light = getIntensities(act_mz,SILACFiltering::acc_lin,SILACFiltering::spline_lin);

	if (checkArea(act_mz, intensities_spl_light, intensities_lin_light))
	{
		intensities.push_back(intensities_spl_light[1]);
		intensities.push_back(intensities_spl_light[2]);
		intensities.push_back(intensities_spl_light[3]);
		act_peak_positions.push_back(act_mz);
		act_peak_positions.push_back(act_mz+isotope_distance);
		act_peak_positions.push_back(act_mz+2*isotope_distance);
	}
	else
	{
		return false;
	}
//	debug=false;
	for (std::set<DoubleReal>::iterator envelope_iterator=envelope_distances.begin();envelope_iterator!=envelope_distances.end();++envelope_iterator)
	{
		DoubleReal envelope_distance=*envelope_iterator;
		 if ( blacklisted(act_mz+envelope_distance) ||  blacklisted(act_mz+envelope_distance+isotope_distance) || blacklisted(act_mz+envelope_distance+2*isotope_distance))
				 return false;
		std::vector<DoubleReal> intensities_spl_envelope = getIntensities(act_mz+envelope_distance,SILACFiltering::acc_spl,SILACFiltering::spline_spl);
		std::vector<DoubleReal> intensities_lin_envelope = getIntensities(act_mz+envelope_distance,SILACFiltering::acc_lin,SILACFiltering::spline_lin);

		std::pair<DoubleReal,DoubleReal> deviation=getPeakCorrelation(act_mz,envelope_distance,0.05);
//		std::cout.precision(10);
//		if (SILACFiltering::feature_id == 229)
//		{
//			std::cout << act_mz << "\t" << deviation.first << "\t" << deviation.second << std::endl;
//		}

//		if ((gsl_spline_eval (SILACFiltering::spline_spl, act_mz-isotope_distance, SILACFiltering::acc_spl) < gsl_spline_eval (SILACFiltering::spline_spl, act_mz, SILACFiltering::acc_spl) && gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_distance-isotope_distance, SILACFiltering::acc_spl) > gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_distance, SILACFiltering::acc_spl)) || (gsl_spline_eval (SILACFiltering::spline_spl, act_mz-isotope_distance, SILACFiltering::acc_spl) > gsl_spline_eval (SILACFiltering::spline_spl, act_mz, SILACFiltering::acc_spl) && gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_distance-isotope_distance, SILACFiltering::acc_spl) < gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_distance, SILACFiltering::acc_spl)))
//		{
//			return false;
//		}

		if (checkArea(act_mz+envelope_distance+deviation.second, intensities_spl_envelope, intensities_lin_envelope))
		{
			intensities.push_back(intensities_spl_envelope[1]);
			intensities.push_back(intensities_spl_envelope[2]);
			intensities.push_back(intensities_spl_envelope[3]);
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
	//TODO
	next_element.silac_type=DataPoint::DOUBLE;
	elements.push_back(next_element);
	peak_positions.clear();
	peak_positions.insert(peak_positions.begin(),act_peak_positions.begin(),act_peak_positions.end());
	return true;

}

bool SILACFilter::checkArea(DoubleReal act_mz, std::vector<DoubleReal>& intensities_spl, std::vector<DoubleReal>& intensities_lin)
{
//		bool condition1 = (); // all six intensities peak simultaneously
		//bool condition2 = (intensities_spl[0]<=intensities_spl[1] && intensities_spl[4]<= intensities_spl[5]) || (intensities_spl[0]>=intensities_spl[1] && intensities_spl[4]>= intensities_spl[5]);

	if ((intensities_lin[1] < SILACFiltering::intensity_cutoff) || (intensities_lin[2] < SILACFiltering::intensity_cutoff) || (intensities_lin[3] < SILACFiltering::intensity_cutoff))
	{
		return false;
	}

	std::pair<DoubleReal,DoubleReal> correlation1=getPeakCorrelation(act_mz,isotope_distance,0.01);
	if (correlation1.first < 0.999)
		return false;
	std::pair<DoubleReal,DoubleReal> correlation2=getPeakCorrelation(act_mz+isotope_distance+correlation1.second,isotope_distance,0.01);
	if (correlation2.first < 0.999)
		return false;

	ExtendedIsotopeModel model;
	Param param;
	param.setValue( "isotope:monoisotopic_mz", act_mz );
	param.setValue( "interpolation_step", 0.01 );
	param.setValue("charge",charge);
	param.setValue("isotope:stdev",getPeakWidth(act_mz)/2);
	model.setParameters( param );
	std::vector<Peak1D> model_data;
	model.getSamples(model_data);

	DoubleReal quality1=(gsl_spline_eval (SILACFiltering::spline_spl, act_mz, SILACFiltering::acc_spl)*model.getIntensity(act_mz+isotope_distance+correlation1.second))/(gsl_spline_eval (SILACFiltering::spline_spl, act_mz+isotope_distance+correlation1.second, SILACFiltering::acc_spl)*model.getIntensity(act_mz));
	DoubleReal quality2=(gsl_spline_eval (SILACFiltering::spline_spl, act_mz+isotope_distance+correlation1.second, SILACFiltering::acc_spl)*model.getIntensity(act_mz+2*isotope_distance+correlation1.second+correlation2.second))/(gsl_spline_eval (SILACFiltering::spline_spl, act_mz+2*isotope_distance+correlation1.second+correlation2.second, SILACFiltering::acc_spl)*model.getIntensity(act_mz+isotope_distance+correlation1.second));

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

	if (std::abs(log(quality1)) > model_deviation || std::abs(log(quality2)) > model_deviation)
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

	std::cout.precision(10);
	DoubleReal area_size=0.3*getPeakWidth(act_mz);
	//			if (SILACFiltering::feature_id == 254)
	//			{
	////				for (DoubleReal pos=act_mz-0.003;pos<=act_mz+0.003;pos+=0.0003)
	////				{
	////					std::cout << pos << "\t" << gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl) << "\t" << gsl_spline_eval (SILACFiltering::spline_spl, pos+isotope_distance+correlation1.second, SILACFiltering::acc_spl)<< "\t" << getPeakCorrelation(pos).first<< "\t" << getPeakCorrelation(pos).second << std::endl;
	////				}
	//				std::cout << act_mz << "\t" << correlation1.first << "\t" << correlation1.second << std::endl;
	//				std::vector<DoubleReal> first_values;
	//				std::vector<DoubleReal> second_values;
	//				for (DoubleReal pos=act_mz-area_size;pos<=act_mz+area_size;pos+=area_size/10)
	//				{
	//					DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl);
	//					DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+correlation1.second+isotope_distance, SILACFiltering::acc_spl);
	//					first_values.push_back(intensity1);
	//					second_values.push_back(intensity2);
	//					std::cout << intensity1 << "\t" << intensity2 << std::endl;
	//				}
	//				DoubleReal act_correlation=Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
	//				std::cout << act_correlation << std::endl;
	//			}
	return true;
}



std::pair<DoubleReal,DoubleReal> SILACFilter::getPeakCorrelation(DoubleReal act_mz,DoubleReal distance,DoubleReal deviation)
{

//	if (SILACFiltering::feature_id==140)
//	std::cout << SILACFiltering::feature_id << " " << act_mz << std::endl;
	DoubleReal best_correlation=0.0;
	DoubleReal best_correlation_deviation=0.0;
	DoubleReal offset=0.0;
	DoubleReal area_size=0.25*getPeakWidth(act_mz);
//	std::cout << act_mz << area_size << std::endl;
	for (int i=0;i<2;++i)
	{
		for (DoubleReal act_deviation=offset+deviation*(-1);act_deviation<=offset+deviation;act_deviation+=deviation/10)
				{
					std::vector<DoubleReal> first_values;
					std::vector<DoubleReal> second_values;
					for (DoubleReal pos=act_mz-area_size;pos<=act_mz+area_size;pos+=area_size/5)
					{
						DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl);
						DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+act_deviation+distance, SILACFiltering::acc_spl);
						first_values.push_back(intensity1);
						second_values.push_back(intensity2);
					}
					DoubleReal act_correlation=Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
					if (act_correlation > best_correlation)
					{
						best_correlation=act_correlation;
						best_correlation_deviation=act_deviation;
					}
				}
		deviation/=20;
		offset=best_correlation_deviation;
	}

//		std::cout << best_correlation_deviation << " ";
	return std::make_pair(best_correlation,best_correlation_deviation);
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
	return 11*(1.889e-7*pow(mz,1.5));
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
