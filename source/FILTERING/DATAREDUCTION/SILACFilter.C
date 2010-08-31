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
#include <OpenMS/CONCEPT/Constants.h>
 #include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_randist.h>
#include <cmath>

namespace OpenMS
{
bool debug=false;
DoubleReal pos1;

SILACFilter::SILACFilter(std::set<DoubleReal> mass_separations, Int charge_,DoubleReal model_deviation_) {
	silac_type=mass_separations.size();
	charge=charge_;
	envelope_distances.insert(0.0);
	for (std::set<DoubleReal>::iterator it=mass_separations.begin();it!=mass_separations.end();++it)
	{
		envelope_distances.insert(*it/(DoubleReal)charge);
	}
	isotope_distance=1.000495/(DoubleReal) charge;
	model_deviation=model_deviation_;
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
}

SILACFilter::SILACFilter(Int charge_,DoubleReal model_deviation_) {
	silac_type=0;
	charge=charge_;
	envelope_distances.insert(0.0);
	isotope_distance=1.000495/(DoubleReal) charge;
	model_deviation=model_deviation_;
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
}

SILACFilter::SILACFilter(Int charge_) {
	silac_type=0;
	charge=charge_;
	envelope_distances.insert(0.0);
	isotope_distance=1.000495/(DoubleReal) charge;
	model_deviation=-1;
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
	blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
}

SILACFilter::~SILACFilter() {
}

bool SILACFilter::blacklisted(DoubleReal value)
{
	bool contained=false;
	for (std::list<std::set<DoubleReal,doubleCmp> >::iterator it=blacklist_lifetime.begin();it!=blacklist_lifetime.end();++it)
	{
		std::set<DoubleReal,doubleCmp>::iterator value_pos=it->find(value);
//		std::cout << (value_pos!=blacklist.end()) << std::endl;
		contained=contained || value_pos!=it->end();
	}
	return contained;
}


void  SILACFilter::blockValue(DoubleReal value)
{
	blacklist_lifetime.back().insert(value);
}


void SILACFilter::reset()
{
	blacklist_lifetime.erase(blacklist_lifetime.begin());

//	blacklist.clear();
//	 for_each( iterators.begin(), iterators.end(), blacklist.erase );
	 blacklist_lifetime.push_back(std::set<DoubleReal,doubleCmp>());
}

bool SILACFilter::isFeature(DoubleReal act_rt,DoubleReal act_mz)
{
	if (act_rt < 1797 && act_rt>1796)
		debug=true;
	else
		debug=false;
	if (blacklisted(act_mz) /*||  blacklisted(act_mz+isotope_distance) || blacklisted(act_mz+2*isotope_distance)*/)
		return false;


	if (gsl_spline_eval (SILACFiltering::spline_lin, act_mz, SILACFiltering::acc_lin) < SILACFiltering::intensity_cutoff)
	{
		return false;
	}

	bool monoisotopic_smaller=true;


	std::vector<DoubleReal> light_intensities;
	std::vector<DoubleReal> intensities;
	std::vector<DoubleReal> act_peak_positions;
	DoubleReal deviation=0.0;
	bool missing_peak=false;


	for (std::set<DoubleReal>::iterator envelope_iterator=envelope_distances.begin();envelope_iterator!=envelope_distances.end();++envelope_iterator)
	{
		DoubleReal envelope_distance=*envelope_iterator;
		DoubleReal tolerance=getPeakWidth(act_mz+envelope_distance);

		DoubleReal max_previous_intensity=0.0;
		DoubleReal max_heavy_intensity=0.0;

		for (DoubleReal pos=act_mz+envelope_distance-0.5*tolerance;pos<=act_mz+envelope_distance+0.5*tolerance;pos+=tolerance/10)
		{
			DoubleReal act_intensity=gsl_spline_eval (SILACFiltering::spline_lin, pos, SILACFiltering::acc_lin);
			DoubleReal previous_intensity=gsl_spline_eval (SILACFiltering::spline_lin,pos-isotope_distance, SILACFiltering::acc_lin);
			if (act_intensity > max_heavy_intensity)
				max_heavy_intensity=act_intensity;
			if (previous_intensity > max_previous_intensity)
				max_previous_intensity=previous_intensity;

		}
		if (envelope_iterator==envelope_distances.begin() && max_previous_intensity>=gsl_spline_eval (SILACFiltering::spline_lin, act_mz, SILACFiltering::acc_lin))
		{
			monoisotopic_smaller=false;
		}
		if (max_heavy_intensity < SILACFiltering::intensity_cutoff)
		{
			return false;
		}

		std::vector<DoubleReal> expected_positions;
		std::vector<DoubleReal> exact_positions;
		std::vector<DoubleReal> data;

		computeCorrelation(act_mz,envelope_distance,tolerance,data);
		DoubleReal exact_position=computeExactPosition(act_mz,envelope_distance,tolerance,gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance, SILACFiltering::acc_lin),data);
		deviation=std::abs(exact_position-envelope_distance);
//		deviation=exact_position;

		if (exact_position < 0.0)
			return false;

		if (envelope_iterator!=envelope_distances.begin() && max_previous_intensity >= gsl_spline_eval (SILACFiltering::spline_lin, act_mz+exact_position, SILACFiltering::acc_lin) && monoisotopic_smaller)
		{
			return false;
		}
		exact_positions.push_back(exact_position);
		data.clear();

		computeCorrelation(act_mz,envelope_distance+isotope_distance,tolerance,data);
		exact_position=computeExactPosition(act_mz,envelope_distance+isotope_distance,tolerance,gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance+isotope_distance, SILACFiltering::acc_lin),data);

		if (exact_position < 0.0)
			return false;

		exact_positions.push_back(exact_position);
		data.clear();

		computeCorrelation(act_mz,envelope_distance+2*isotope_distance,tolerance,data);
		exact_position=computeExactPosition(act_mz,envelope_distance+2*isotope_distance,tolerance,gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance+2*isotope_distance, SILACFiltering::acc_lin),data);

		double act_missing_peak=false;
		if (exact_position < 0.0 && missing_peak)
			return false;
		else if (exact_position < 0.0 && !missing_peak)
		{
			act_missing_peak=true;
		}
		exact_positions.push_back(exact_position);

		Size peak_number=3;
		do
		{
			DoubleReal further_intensity=gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance+peak_number*isotope_distance, SILACFiltering::acc_lin);
			if (further_intensity > gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance+(peak_number-1)*isotope_distance, SILACFiltering::acc_lin))
				break;
			computeCorrelation(act_mz,envelope_distance+peak_number*isotope_distance,tolerance,data);
			exact_position=computeExactPosition(act_mz,envelope_distance+peak_number*isotope_distance,tolerance,further_intensity,data);
			if (exact_position > 0.0)
			{
				exact_positions.push_back(exact_position);
				++peak_number;
			}
		}while(exact_position > 0.0);


		bool correct_ratios=true;
		if (envelope_iterator!=envelope_distances.begin())
			correct_ratios=checkRatios(act_mz,light_intensities,exact_positions);
		else
			light_intensities.insert(light_intensities.end(),exact_positions.begin(),exact_positions.end());

		if (checkArea(act_mz, exact_positions, intensities, missing_peak) && correct_ratios)
		{
			for (Size i=0;i<exact_positions.size();++i)
			{
				exact_positions[i]+=act_mz;
			}
			act_peak_positions.insert(act_peak_positions.end(),exact_positions.begin(),exact_positions.end());
		}
		else
		{
			return false;
		}
		missing_peak=act_missing_peak;
	}
	DataPoint next_element;
	next_element.feature_id=SILACFiltering::feature_id;
	next_element.rt=act_rt;
	next_element.mz=act_mz;
	next_element.silac_type=envelope_distances.size();
	next_element.charge=charge;
	next_element.envelope_distance_light_heavy=*envelope_distances.rbegin();
	//TODO !!
	next_element.intensities.insert(next_element.intensities.end(),intensities.begin(),intensities.end());
	//TODO
	next_element.silac_type=DataPoint::DOUBLE;
	elements.push_back(next_element);
	peak_positions.clear();
	peak_positions.insert(peak_positions.begin(),act_peak_positions.begin(),act_peak_positions.end());

	peak.setMZ(act_mz);
	peak.setIntensity(deviation*10);
	return true;

}

bool SILACFilter::checkArea(DoubleReal act_mz, const std::vector<DoubleReal>& exact_positions, std::vector<DoubleReal>& intensities,bool missing_peak)
{
	std::vector<DoubleReal> act_intensities;
	for (std::vector<DoubleReal>::const_iterator position_it=exact_positions.begin();position_it!=exact_positions.end();++position_it)
	{
		act_intensities.push_back(gsl_spline_eval (SILACFiltering::spline_spl, act_mz+*position_it, SILACFiltering::acc_spl));
	}

	DoubleReal area_width=getPeakWidth(act_mz);
	for (Size i=1;i< 3;++i)
	{
		std::vector<DoubleReal> first_values;
		std::vector<DoubleReal> second_values;
		DoubleReal max1=0.0;
		DoubleReal max2=0.0;
		for (DoubleReal pos=act_mz-0.3*area_width;pos<=act_mz+0.3*area_width;pos+=0.06*area_width)
		{
			DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl);
			DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+exact_positions[i], SILACFiltering::acc_spl);
			if (intensity1>max1)
				max1=intensity1;
			if (intensity2>max2)
				max2=intensity2;
		}
		DoubleReal scaling_factor=max1/max2;
		for (DoubleReal pos=act_mz-0.5*area_width;pos<=act_mz+0.5*area_width;pos+=0.05*area_width)
		{
			DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::acc_spl);
			DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+exact_positions[i], SILACFiltering::acc_spl);
			first_values.push_back(intensity1);
			second_values.push_back(intensity2*scaling_factor);
		}
		DoubleReal act_correlation=Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
//		if (debug && act_mz > 616.322 &&  act_mz < 616.352 && exact_positions[0]<0.01 && i==1)
//		{
//			std::cout << std::setprecision(10);
//			std::cout << act_mz << "\t" << act_intensities[0] << "\t" << act_intensities[1] << "\t"<< std::setprecision(3) << act_correlation << std::endl;
//		}

		if ((act_correlation < 0.95 && i<2) || (i==2 && missing_peak && act_correlation < 0.95))
			return false;
	}

	ExtendedIsotopeModel model;
	Param param;
	param.setValue( "isotope:monoisotopic_mz", act_mz+exact_positions[0] );
	param.setValue("charge",charge);
	param.setValue("isotope:stdev",getPeakWidth(act_mz+exact_positions[0])*0.42466);
	model.setParameters( param );
	std::vector<Peak1D> model_data;
	model.getSamples(model_data);

	for (Size i=0;i< 2;++i)
	{
		DoubleReal quality=(act_intensities[i]*model.getIntensity(act_mz+exact_positions[i+1]))/(act_intensities[i+1]*model.getIntensity(act_mz+exact_positions[i]));
		if ((std::abs(log(quality)) > model_deviation && i!=1) || (std::abs(log(quality)) > model_deviation && missing_peak && i==1))
			return false;
	}

	intensities.insert(intensities.end(),act_intensities.begin(),act_intensities.end());
	return true;
}

void SILACFilter::computeCorrelation(DoubleReal act_mz,DoubleReal offset,DoubleReal tolerance,std::vector<DoubleReal>& data)
{
	if (offset > 0.0)
	{
		DoubleReal stepwidth=tolerance/28;
		//n must be a power of two for gsl fast fourier transformation; take next higher size for n, which is a power of two and fill the rest with zeros
		Size vector_size = pow(2,(ceil(log2((3*tolerance)/stepwidth))));
		data.clear();
		data.resize(vector_size,0);

			//Create a vector containing interpolated values with a spacing of "stepwidth"
			Size i=0;
			DoubleReal starting_offset=std::min(act_mz-SILACFiltering::mz_min,tolerance);

			for (DoubleReal x=act_mz-starting_offset;x<=act_mz+tolerance;x+=stepwidth)
			{
				data[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x, SILACFiltering::acc_lin);
				++i;
			}

			std::vector<DoubleReal> gauss_fitted_data(vector_size,0);
			i=0;
			for (DoubleReal x=act_mz-starting_offset-tolerance;x<=act_mz+tolerance;x+=stepwidth)
			{
				gauss_fitted_data[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x+offset, SILACFiltering::acc_lin)/**gsl_ran_gaussian_pdf (x-act_mz, (tolerance/3)*0.75)*gauss_normalization_factor*/;
				++i;
			}
//			if (SILACFiltering::feature_id>96 && SILACFiltering::feature_id<118 &&offset>3.8 && offset<4.2)
//							std::cout << "------------------" << std::endl;

			//Fourier transformation of the values
			gsl_fft_real_radix2_transform (&(*data.begin()), 1, vector_size);
			gsl_fft_real_radix2_transform (&(*gauss_fitted_data.begin()), 1, vector_size);

			//Multiply the fourier transformed complex values with the complex conjugate
			//Have a look at the GNU Scientific Library reference manual for a description of the data structure
			data[0]=data[0]*gauss_fitted_data[0];
			data[vector_size/2]=data[vector_size/2]*gauss_fitted_data[vector_size/2];
			for (i = 1; i <= vector_size/2; ++i)
			{
				data[i]=data[i]*gauss_fitted_data[i]+data[vector_size-i]*gauss_fitted_data[vector_size-i];
				data[vector_size-i]=0.0;
			}

			//Compute inverse fourier transformation
			gsl_fft_halfcomplex_radix2_inverse (&(*data.begin()), 1, vector_size);

			data.resize(vector_size/2);
	}
}


DoubleReal SILACFilter::computeExactPosition(DoubleReal act_mz,DoubleReal expected_position,DoubleReal tolerance,DoubleReal cutoff,std::vector<DoubleReal> data)
{
	if (expected_position>0)
	{
		DoubleReal stepwidth=tolerance/28;
//		std::cout << stepwidth << std::endl;

		std::vector<DoubleReal> x(data.size(),0);
		std::vector<DoubleReal> y(data.size(),0);
		DoubleReal shift=expected_position-tolerance;
		for (Size i=0;i<data.size();++i)
		{
			x[i]=shift;
			y[i]=data[i];
			shift+=stepwidth;
		}

		gsl_interp_accel* acc_correlation = gsl_interp_accel_alloc();
		gsl_spline* spline_correlation = gsl_spline_alloc(gsl_interp_cspline, x.size());
		gsl_spline_init(spline_correlation, &(*x.begin()), &(*y.begin()), x.size());

		DoubleReal exact_position=0.0;
		DoubleReal max_intensity=0.0;

		for (DoubleReal act_position=0.0;act_position<=tolerance; act_position+=stepwidth)
		{
			DoubleReal last_intensity=gsl_spline_eval (spline_correlation, expected_position+act_position-stepwidth, acc_correlation);
			DoubleReal act_intensity=gsl_spline_eval (spline_correlation, expected_position+act_position, acc_correlation);
			DoubleReal next_intensity=gsl_spline_eval (spline_correlation, expected_position+act_position+stepwidth, acc_correlation);
			if (act_intensity > last_intensity && act_intensity > next_intensity && act_intensity > 1000 && gsl_spline_eval (SILACFiltering::spline_lin, act_mz+expected_position+act_position, SILACFiltering::acc_lin) > SILACFiltering::intensity_cutoff)
			{
				gsl_spline_free(spline_correlation);
				gsl_interp_accel_free(acc_correlation);
				return expected_position+act_position;
			}
			last_intensity=gsl_spline_eval(spline_correlation, expected_position-act_position-stepwidth, acc_correlation);
			act_intensity=gsl_spline_eval(spline_correlation, expected_position-act_position, acc_correlation);
			next_intensity=gsl_spline_eval(spline_correlation, expected_position-act_position+stepwidth, acc_correlation);
			if (act_intensity > last_intensity && act_intensity > next_intensity && act_intensity > 1000 && gsl_spline_eval (SILACFiltering::spline_lin, act_mz+expected_position-act_position, SILACFiltering::acc_lin) > SILACFiltering::intensity_cutoff)
			{
				gsl_spline_free(spline_correlation);
				gsl_interp_accel_free(acc_correlation);
				return expected_position-act_position;
			}
			//			if (SILACFiltering::feature_id==2889)
			//																		{
			//																			std::cout << act_position << "\t" << act_intensity << std::endl;
			//																		}
		}
		gsl_spline_free(spline_correlation);
		gsl_interp_accel_free(acc_correlation);
		return -1;
	}
	else
	{
		return 0.0;
	}
}

bool SILACFilter::checkRatios(DoubleReal act_mz,const std::vector<DoubleReal>& light_positions, const std::vector<DoubleReal>& envelope_positions)
{
	std::vector<DoubleReal> light_intensities(3,0);
	light_intensities[0]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+light_positions[0], SILACFiltering::acc_spl);
	light_intensities[1]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+light_positions[1], SILACFiltering::acc_spl);
	light_intensities[2]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+light_positions[2], SILACFiltering::acc_spl);

	std::vector<DoubleReal> envelope_intensities(3,0);
	envelope_intensities[0]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_positions[0], SILACFiltering::acc_spl);
	envelope_intensities[1]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_positions[1], SILACFiltering::acc_spl);
	envelope_intensities[2]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+envelope_positions[2], SILACFiltering::acc_spl);

	DoubleReal ratio1=log((light_intensities[0]/light_intensities[1])/(envelope_intensities[0]/envelope_intensities[1]));
	DoubleReal ratio2=log((light_intensities[1]/light_intensities[2])/(envelope_intensities[1]/envelope_intensities[2]));

//	std::cout << envelope_intensities[0] << "\t" << envelope_intensities[1] << "\t" << (envelope_intensities[0]/envelope_intensities[1]) << std::endl;

	peak_ratio1.setMZ(act_mz);
	peak_ratio1.setIntensity(std::abs(ratio1));
	peak_ratio2.setMZ(act_mz);
	peak_ratio2.setIntensity(std::abs(ratio2));

//	return (std::abs(ratio1) < 1.5 && std::abs(ratio2) < 1);
	return true;
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
	DoubleReal peak_width=3*getPeakWidth((a+b)/2);
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
	return 5*(1.889e-7*pow(mz,1.5));
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
