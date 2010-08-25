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
	if (act_rt < 2001)
		debug=true;
	else
		debug=false;
	if (blacklisted(act_mz) ||  blacklisted(act_mz+isotope_distance) || blacklisted(act_mz+2*isotope_distance))
		return false;

	std::vector<DoubleReal> intensities;
	std::vector<DoubleReal> act_peak_positions;
	DoubleReal deviation=0.0;


	for (std::set<DoubleReal>::iterator envelope_iterator=envelope_distances.begin();envelope_iterator!=envelope_distances.end();++envelope_iterator)
	{
		DoubleReal envelope_distance=*envelope_iterator;
		DoubleReal tolerance=getPeakWidth(act_mz+envelope_distance);
		std::vector<DoubleReal> expected_positions;
		std::vector<DoubleReal> exact_positions;
		std::vector<DoubleReal> data;

		computeCorrelation(act_mz,envelope_distance,tolerance,data);
		DoubleReal exact_position=computeExactPosition(envelope_distance,tolerance,gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance, SILACFiltering::acc_lin),data);
		deviation=std::abs(exact_position-envelope_distance);
//		deviation=exact_position;
		if (exact_position < 0.0)
			return false;
		else
			exact_positions.push_back(exact_position);
		data.clear();

		computeCorrelation(act_mz,envelope_distance+isotope_distance,tolerance,data);
		exact_position=computeExactPosition(envelope_distance+isotope_distance,tolerance,gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance+isotope_distance, SILACFiltering::acc_lin),data);

		if (exact_position < 0.0)
			return false;
		else
			exact_positions.push_back(exact_position);
		data.clear();

		computeCorrelation(act_mz,envelope_distance+2*isotope_distance,tolerance,data);
		exact_position=computeExactPosition(envelope_distance+2*isotope_distance,tolerance,gsl_spline_eval (SILACFiltering::spline_lin, act_mz+envelope_distance+2*isotope_distance, SILACFiltering::acc_lin),data);
		if (exact_position < 0.0)
			return false;
		else
			exact_positions.push_back(exact_position);

		if (checkArea(act_mz, exact_positions, intensities))
		{
			act_peak_positions.insert(act_peak_positions.end(),exact_positions.begin(),exact_positions.end());
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
	peak_positions.insert(peak_positions.begin(),++act_peak_positions.begin(),act_peak_positions.end());

	peak.setMZ(act_mz);
	peak.setIntensity(deviation*10);
	return true;

}

bool SILACFilter::checkArea(DoubleReal act_mz, const std::vector<DoubleReal>& exact_positions, std::vector<DoubleReal>& intensities)
{
	std::vector<DoubleReal> act_intensities(3,0);
	std::vector<DoubleReal> mono_intensities;

	act_intensities[0]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+exact_positions[0], SILACFiltering::acc_spl);
	act_intensities[1]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+exact_positions[1], SILACFiltering::acc_spl);
	act_intensities[2]=gsl_spline_eval (SILACFiltering::spline_spl, act_mz+exact_positions[2], SILACFiltering::acc_spl);

//	if (debug)
//		{
//			std::cout.precision(10);
//			std::cout << act_mz << " " << act_intensities[0] << " " << act_intensities[1] << " " << act_intensities[2] << std::endl;
//				for (std::vector<DoubleReal>::const_iterator position_it=exact_positions.begin();position_it!=exact_positions.end();++position_it)
//				{
//					std::vector<DoubleReal> first_values;
//					std::vector<DoubleReal> second_values;
//					for (DoubleReal pos=act_mz-0.02;pos<=act_mz+0.02;pos+=0.001)
//					{
//						DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_lin, pos, SILACFiltering::acc_lin);
//						DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_lin, pos+*position_it, SILACFiltering::acc_lin);
//						first_values.push_back(intensity1);
//						second_values.push_back(intensity2);
//					}
//					DoubleReal act_correlation=Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
//
//					std::cout << *position_it << "[" << act_correlation << "] ";
//				}
//				std::cout << std::endl;
//		}

	if ((act_intensities[0] < SILACFiltering::intensity_cutoff))
	{
		return false;
	}

	DoubleReal area_width=getPeakWidth(act_mz);
	for (std::vector<DoubleReal>::const_iterator position_it=exact_positions.begin();position_it!=exact_positions.end();++position_it)
	{
		std::vector<DoubleReal> first_values;
		std::vector<DoubleReal> second_values;
		for (DoubleReal pos=act_mz-0.5*area_width;pos<=act_mz+0.5*area_width;pos+=0.001)
		{
			DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_lin, pos, SILACFiltering::acc_lin);
			DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_lin, pos+*position_it, SILACFiltering::acc_lin);
			first_values.push_back(intensity1);
			second_values.push_back(intensity2);
		}
		DoubleReal act_correlation=Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
		if (act_correlation < 0.99)
			return false;
		intensities.push_back(gsl_spline_eval (SILACFiltering::spline_spl, act_mz+*position_it, SILACFiltering::acc_spl));
	}

//	if (SILACFiltering::feature_id==420)
//										{
//											std::cout << act_mz << std::endl;
//										}


	ExtendedIsotopeModel model;
	Param param;
	param.setValue( "isotope:monoisotopic_mz", act_mz+exact_positions[0] );
	param.setValue( "interpolation_step", 0.01 );
	param.setValue("charge",charge);
	param.setValue("isotope:stdev",getPeakWidth(act_mz+exact_positions[0])/2);
	model.setParameters( param );
	std::vector<Peak1D> model_data;
	model.getSamples(model_data);

	DoubleReal quality1=(act_intensities[0]*model.getIntensity(act_mz+exact_positions[1]))/(act_intensities[1]*model.getIntensity(act_mz+exact_positions[0]));
	DoubleReal quality2=(act_intensities[1]*model.getIntensity(act_mz+exact_positions[2]))/(act_intensities[2]*model.getIntensity(act_mz+exact_positions[1]));

	if (std::abs(log(quality1)) > model_deviation || std::abs(log(quality2)) > model_deviation)
		return false;

	//	False positive debug output
//	if (SILACFiltering::feature_id == 275)
//	{
//	std::cout <<
//
//	}
//
//	std::cout.precision(10);
	return true;
}

void SILACFilter::computeCorrelation(DoubleReal act_mz,DoubleReal offset,DoubleReal tolerance,std::vector<DoubleReal>& data)
{
	if (offset > 0.0)
	{
		DoubleReal stepwidth=0.001;
		//n must be a power of two for gsl fast fourier transformation; take next higher size for n, which is a power of two and fill the rest with zeros
		Size vector_size = pow(2,(ceil(log2((3*tolerance)/stepwidth))));
		data.clear();
		data.resize(vector_size,0);

			//Create a vector containing interpolated values with a spacing of "stepwidth"
			Size i=0;
			DoubleReal starting_offset=std::min(act_mz-SILACFiltering::mz_min,tolerance);

//			if (SILACFiltering::feature_id>96 && SILACFiltering::feature_id<118 &&offset>3.8 && offset<4.2)
//				{
//					std::cout << act_mz << std::endl;
//					std::cout << tolerance << std::endl;
//					std::cout << offset << std::endl;
//					std::cout << "id: " << SILACFiltering::feature_id << std::endl;
//				}



			for (DoubleReal x=act_mz-starting_offset;x<=act_mz+tolerance;x+=stepwidth)
			{
//				if (SILACFiltering::feature_id>96 && SILACFiltering::feature_id<118 &&offset>3.8 && offset<4.2)
//					std::cout << x << "\t" << gsl_spline_eval (SILACFiltering::spline_lin, x, SILACFiltering::acc_lin) << std::endl;
				data[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x, SILACFiltering::acc_lin);
				++i;
			}

//			if (SILACFiltering::feature_id>96 && SILACFiltering::feature_id<118)
//							std::cout << "\n\n" << std::endl;



//			DoubleReal gauss_normalization_factor=(tolerance/3)*0.75*sqrt(2*Constants::PI);
			std::vector<DoubleReal> gauss_fitted_data(vector_size,0);
			i=0;
			for (DoubleReal x=act_mz-starting_offset-tolerance;x<=act_mz+tolerance;x+=stepwidth)
			{
//				if (SILACFiltering::feature_id>96 && SILACFiltering::feature_id<118 &&offset>3.8 && offset<4.2)
//									std::cout << x+offset << "\t" << gsl_spline_eval (SILACFiltering::spline_lin, x+offset, SILACFiltering::acc_lin) << std::endl;
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


DoubleReal SILACFilter::computeExactPosition(DoubleReal expected_position,DoubleReal tolerance,DoubleReal cutoff,std::vector<DoubleReal> data)
{
	if (expected_position>0)
	{
		DoubleReal stepwidth=0.001;

		std::vector<DoubleReal> x(data.size(),0);
		std::vector<DoubleReal> y(data.size(),0);
		DoubleReal shift=expected_position-tolerance;
		for (Size i=0;i<data.size();++i)
		{
			x[i]=shift;
			y[i]=data[i];
			shift+=stepwidth;
		}

		gsl_interp_accel* acc_spl = gsl_interp_accel_alloc();
		gsl_spline* spline_spl = gsl_spline_alloc(gsl_interp_akima, x.size());
		gsl_spline_init(spline_spl, &(*x.begin()), &(*y.begin()), x.size());

		DoubleReal exact_position=0.0;
		DoubleReal max_intensity=0.0;

		for (DoubleReal act_position=expected_position-tolerance;act_position<=expected_position+tolerance; act_position+=stepwidth)
		{

			DoubleReal act_intensity=gsl_spline_eval (spline_spl, act_position, acc_spl);
//			if (SILACFiltering::feature_id==2889)
//																		{
//																			std::cout << act_position << "\t" << act_intensity << std::endl;
//																		}

			if (act_intensity > max_intensity)
			{
				exact_position=act_position;
				max_intensity=act_intensity;
			}
		}
		if (gsl_spline_eval (spline_spl, exact_position-stepwidth, acc_spl) > max_intensity || gsl_spline_eval (spline_spl, exact_position+stepwidth, acc_spl) > max_intensity || max_intensity < 1000 )
		{
			gsl_spline_free(spline_spl);
			gsl_interp_accel_free(acc_spl);
			return -1;
		}
		else
		{
//			if (SILACFiltering::feature_id==420)
//			{
//				std::cout << "expe_pos " << expected_position << std::endl;
//				std::cout << max_intensity << "\t" << exact_position << "\t(" << expected_position-tolerance <<"\t" << expected_position+tolerance<< ")\t" << cutoff << std::endl;
//			}

//			if (SILACFiltering::feature_id>2801 && SILACFiltering::feature_id<2806 &&expected_position>3.8 && expected_position<4.2)
//			{
//				for (DoubleReal act_position=expected_position-tolerance;act_position<=expected_position+tolerance; act_position+=stepwidth)
//				{
//
//					DoubleReal act_intensity=gsl_spline_eval (spline_spl, act_position, acc_spl);
//					std::cout << act_position << "\t" << act_intensity << std::endl;
//				}
//				std::cout << "pos: " << exact_position << "\n#############################################################" << std::endl;
//			}
//			if (SILACFiltering::feature_id>96 && SILACFiltering::feature_id<118)
//				std::cout << "\n\n";
			gsl_spline_free(spline_spl);
			gsl_interp_accel_free(acc_spl);
			return exact_position;
		}
	}
	else
	{
		return 0.0;
	}
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
	DoubleReal peak_width=2*getPeakWidth((a+b)/2);
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
