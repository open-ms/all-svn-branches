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
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
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
	SILACFilter::SILACFilter(std::set<DoubleReal> mass_separations, Int charge_,DoubleReal model_deviation_, Int isotopes_per_peptide_)
	{
		numberOfPeptides = mass_separations.size(); // number of labelled peptides +1 [e.g. for SILAC triplet =3]
		isotopes_per_peptide=isotopes_per_peptide_;
		charge=charge_;
		
		mz_peptide_separations.push_back(0.0);
		for (std::set<DoubleReal>::iterator it=mass_separations.begin(); it!=mass_separations.end(); ++it)
		{
			mz_peptide_separations.push_back(*it/(DoubleReal)charge);
		}
		
		
		
		
		
		
		
		silac_type=mass_separations.size(); // rausschmeissen??
		isotope_distance=1.000495/(DoubleReal) charge;
		model_deviation=model_deviation_;
	}

SILACFilter::~SILACFilter() {
}

bool SILACFilter::isPair(DoubleReal rt, DoubleReal mz)
{
	//---------------------------------------------------------------
	// BLUNT INTENSITY FILTER (Just check that intensity at current position is above the intensity cutoff.)
	//---------------------------------------------------------------
	if (gsl_spline_eval (SILACFiltering::spline_lin, mz, SILACFiltering::current_lin) < SILACFiltering::intensity_cutoff)
	{
		return false;
	}

	//---------------------------------------------------------------
	// EXACT m/z POSITIONS (Determine the actual positions of peaks. Say 4 Th is the theoretic shift. In the experimental data it will be 4.0029 Th.)
	//---------------------------------------------------------------
	std::vector<std::vector<DoubleReal> > exact_positions;
	for (Int peptide = 0; peptide <= numberOfPeptides; peptide++) // loop over labelled peptides [e.g. for SILAC triplets: 0=light 1=medium 2=heavy]
	{
		for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++) // loop over isotopic peaks within a peptide [0=mono-isotopic peak etc.]
		{
			std::vector<DoubleReal> tempVector;
			computeCorrelation(mz, mz_peptide_separations[peptide] + isotope*isotope_distance, 0.001, tempVector);
			DoubleReal deltaMZ = computeExactDistance(mz, mz_peptide_separations[peptide] + isotope*isotope_distance, 0.001, tempVector);
			std::cout << " m/z = " << mz << ",    delta m/z = " << mz_peptide_separations[peptide] + isotope*isotope_distance << ",    exact delta m/z = " << deltaMZ << std::endl;
		}
	}
	
	

	bool monoisotopic_smaller=true;
	std::vector<std::vector<DoubleReal> > intensities;
	peak_positions.clear();
	bool missing_peak=false;

	//Iterate over all mass shifts (including 0) and check each isotope pattern
	for (std::vector<DoubleReal>::iterator envelope_iterator=mz_peptide_separations.begin();envelope_iterator!=mz_peptide_separations.end();++envelope_iterator)
	{
		DoubleReal envelope_distance=*envelope_iterator;
		//Estimate the correlation tolerance
		DoubleReal tolerance=getPeakWidth(mz+envelope_distance);

		DoubleReal max_previous_intensity=0.0;
		DoubleReal max_heavy_intensity=0.0;

		//Check if the intensity of the monoisotopic peak of each isotope pattern is higher than its predecessor; use the tolerance area for the intensity maximum
		for (DoubleReal pos=mz+envelope_distance-0.5*tolerance;pos<=mz+envelope_distance+0.5*tolerance;pos+=tolerance/10)
		{
			DoubleReal act_intensity=gsl_spline_eval (SILACFiltering::spline_lin, pos, SILACFiltering::current_lin);
			DoubleReal previous_intensity=gsl_spline_eval (SILACFiltering::spline_lin,pos-isotope_distance, SILACFiltering::current_lin);
			if (act_intensity > max_heavy_intensity)
				max_heavy_intensity=act_intensity;
			if (previous_intensity > max_previous_intensity)
				max_previous_intensity=previous_intensity;
		}
		//Check if the light monoisotopic peak is smaller than its predecessor
		if (envelope_iterator==mz_peptide_separations.begin() && max_previous_intensity>=gsl_spline_eval (SILACFiltering::spline_lin, mz, SILACFiltering::current_lin))
		{
			monoisotopic_smaller=false;
		}
		//Check if the intensity of the current mass shift is higher than the threshold
		if (max_heavy_intensity < SILACFiltering::intensity_cutoff)
		{
			return false;
		}

		std::vector<DoubleReal> expected_positions;
		std::vector<DoubleReal> exact_positions;
		std::vector<DoubleReal> data;

		//Compute autocorrelation and exact positions for all peaks of the SILAC pattern
		//Monoisotopic peak
		computeCorrelation(mz,envelope_distance,tolerance,data);
		DoubleReal exact_position=computeExactDistance(mz,envelope_distance,tolerance,data);

		if (exact_position < 0.0)
			return false;

		//Check if the intensities of the predecessors of the light monoisotopic peak and the monoisotopic peak of the current mass shift are both either smaller or higher
		if (envelope_iterator!=mz_peptide_separations.begin() && max_previous_intensity >= gsl_spline_eval (SILACFiltering::spline_lin, mz+exact_position, SILACFiltering::current_lin) && monoisotopic_smaller)
		{
			return false;
		}
		exact_positions.push_back(exact_position);
		data.clear();

		//First isotope peak
		computeCorrelation(mz,envelope_distance+isotope_distance,tolerance,data);
		exact_position=computeExactDistance(mz,envelope_distance+isotope_distance,tolerance,data);

		if (exact_position < 0.0)
			return false;

		exact_positions.push_back(exact_position);
		data.clear();

		// second isotope peak
		computeCorrelation(mz, envelope_distance + 2 * isotope_distance, tolerance, data);
		exact_position=computeExactDistance(mz, envelope_distance + 2 * isotope_distance, tolerance, data);

		//Check for missing peaks. One second isotope peak may be missing.
		//If the current peak is missing and there is already a missing peak, return false
		double act_missing_peak=false;
		if (exact_position < 0.0 && missing_peak)
			return false;
		else if (exact_position < 0.0 && !missing_peak)
		{
			act_missing_peak=true;
		}
		exact_positions.push_back(exact_position);

		//Add the current m/z position the exact positions of the peaks, to get the position in the spectrum
		std::vector<DoubleReal> act_intensities;

		//Check the shape of the current isotope pattern
		if (checkPattern(mz, exact_positions, act_intensities, missing_peak))
		{
			for (Size i=0;i<exact_positions.size();++i)
			{
				exact_positions[i]+=mz;
			}
		}
		else
		{
			return false;
		}
		//store all peak positions for blacklisting
		peak_positions.insert(peak_positions.end(),exact_positions.begin(),exact_positions.end());
		missing_peak=act_missing_peak;
		intensities.push_back(act_intensities);
	}
	//After all filters are passed and all distances are exact, create a data point at the given position
	DataPoint next_element;
	next_element.feature_id=SILACFiltering::feature_id;
	next_element.rt=rt;
	next_element.mz=mz;
	next_element.charge=charge;
	next_element.intensities.insert(next_element.intensities.end(),intensities.begin(),intensities.end());
	next_element.mass_shifts.insert(next_element.mass_shifts.begin(),mz_peptide_separations.begin(),mz_peptide_separations.end());
	elements.push_back(next_element);
	return true;

}

bool SILACFilter::checkPattern(DoubleReal mz, const std::vector<DoubleReal>& exact_positions, std::vector<DoubleReal>& intensities,bool missing_peak)
{
	std::vector<DoubleReal> act_intensities;
	//Store all intensities of the current isotope pattern
	for (std::vector<DoubleReal>::const_iterator position_it=exact_positions.begin();position_it!=exact_positions.end();++position_it)
	{
		DoubleReal act_intensity=gsl_spline_eval (SILACFiltering::spline_spl, mz+*position_it, SILACFiltering::current_spl);
		if (act_intensity<=0)
			return false;
		act_intensities.push_back(act_intensity);
	}

	//Compute the pearson correlation of all peaks within the current isotope pattern in the area +/- 0.7*peak_width around the current position
	DoubleReal area_width=getPeakWidth(mz);
	for (Size i=1;i< 3;++i)
	{
		std::vector<DoubleReal> first_values;
		std::vector<DoubleReal> second_values;
		for (DoubleReal pos=mz-0.7*area_width;pos<=mz+0.7*area_width;pos+=0.14*area_width)
		{
			DoubleReal intensity1=gsl_spline_eval (SILACFiltering::spline_spl, pos, SILACFiltering::current_spl);
			DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+exact_positions[i], SILACFiltering::current_spl);
			first_values.push_back(intensity1);
			second_values.push_back(intensity2);
		}
		DoubleReal act_correlation=Math::pearsonCorrelationCoefficient(first_values.begin(), first_values.end(), second_values.begin(), second_values.end());
		//Do not take a missing peak into account
		if ((act_correlation < 0.95 && i<2) || (i==2 && missing_peak && act_correlation < 0.95))
			return false;
	}

	//Isotope model filtering
	ExtendedIsotopeModel model;
	Param param;
	param.setValue( "isotope:monoisotopic_mz", mz+exact_positions[0] );
	param.setValue("charge",charge);
	param.setValue("isotope:stdev",getPeakWidth(mz+exact_positions[0])*0.42466);
	model.setParameters( param );
	std::vector<Peak1D> model_data;
	model.getSamples(model_data);

	for (Size i=0;i< 2;++i)
	{
		DoubleReal quality=(act_intensities[i]*model.getIntensity(mz+exact_positions[i+1]))/(act_intensities[i+1]*model.getIntensity(mz+exact_positions[i]));
		if ((std::abs(log(quality)) > model_deviation && i!=1) || (std::abs(log(quality)) > model_deviation && missing_peak && i==1))
			return false;
	}
	intensities.insert(intensities.end(),act_intensities.begin(),act_intensities.end());
	return true;
}

void SILACFilter::computeCorrelation(DoubleReal mz,DoubleReal offset,DoubleReal tolerance,std::vector<DoubleReal>& data)
{
	if (offset > 0.0)
	{
		DoubleReal stepwidth=tolerance/18;
		//n must be a power of two for gsl fast fourier transformation; take next higher size for n, which is a power of two and fill the rest with zeros
		Size vector_size = pow(2,(ceil(log2((3*tolerance)/stepwidth))));
		data.clear();
		data.resize(vector_size,0);

			//Create a vector containing interpolated values with a spacing of "stepwidth"
			Size i=0;
			DoubleReal starting_offset=std::min(mz-SILACFiltering::mz_min,tolerance);

			for (DoubleReal x=mz-starting_offset;x<=mz+tolerance;x+=stepwidth)
			{
				data[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x, SILACFiltering::current_lin);
				++i;
			}

			std::vector<DoubleReal> gauss_fitted_data(vector_size,0);
			i=0;
			for (DoubleReal x=mz-starting_offset-tolerance;x<=mz+tolerance;x+=stepwidth)
			{
				gauss_fitted_data[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x+offset, SILACFiltering::current_lin)/**gsl_ran_gaussian_pdf (x-mz, (tolerance/3)*0.75)*gauss_normalization_factor*/;
				++i;
			}

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


DoubleReal SILACFilter::computeExactDistance(DoubleReal mz,DoubleReal expected_distance,DoubleReal tolerance,std::vector<DoubleReal> data)
{
	if (expected_distance>0)
	{
		DoubleReal stepwidth=tolerance/18;

		//Interpolate the current autocorrelation peak
		std::vector<DoubleReal> x(data.size(),0);
		std::vector<DoubleReal> y(data.size(),0);
		DoubleReal shift=expected_distance-tolerance;
		for (Size i=0;i<data.size();++i)
		{
			x[i]=shift;
			y[i]=data[i];
			shift+=stepwidth;
		}

		gsl_interp_accel* acc_correlation = gsl_interp_accel_alloc();
		gsl_spline* spline_correlation = gsl_spline_alloc(gsl_interp_cspline, x.size());
		gsl_spline_init(spline_correlation, &(*x.begin()), &(*y.begin()), x.size());

		//Search for first maximum in + and - direction and check preconditions
		for (DoubleReal act_position=0.0;act_position<=tolerance; act_position+=stepwidth)
		{
			DoubleReal last_intensity=gsl_spline_eval (spline_correlation, expected_distance+act_position-stepwidth, acc_correlation);
			DoubleReal act_intensity=gsl_spline_eval (spline_correlation, expected_distance+act_position, acc_correlation);
			DoubleReal next_intensity=gsl_spline_eval (spline_correlation, expected_distance+act_position+stepwidth, acc_correlation);
			if (act_intensity > last_intensity && act_intensity > next_intensity && act_intensity > 1000 && gsl_spline_eval (SILACFiltering::spline_lin, mz+expected_distance+act_position, SILACFiltering::current_lin) > SILACFiltering::intensity_cutoff)
			{
				gsl_spline_free(spline_correlation);
				gsl_interp_accel_free(acc_correlation);
				return expected_distance+act_position;
			}
			last_intensity=gsl_spline_eval(spline_correlation, expected_distance-act_position-stepwidth, acc_correlation);
			act_intensity=gsl_spline_eval(spline_correlation, expected_distance-act_position, acc_correlation);
			next_intensity=gsl_spline_eval(spline_correlation, expected_distance-act_position+stepwidth, acc_correlation);
			if (act_intensity > last_intensity && act_intensity > next_intensity && act_intensity > 1000 && gsl_spline_eval (SILACFiltering::spline_lin, mz+expected_distance-act_position, SILACFiltering::current_lin) > SILACFiltering::intensity_cutoff)
			{
				gsl_spline_free(spline_correlation);
				gsl_interp_accel_free(acc_correlation);
				return expected_distance-act_position;
			}
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

DoubleReal SILACFilter::computeActualMzShift(DoubleReal mz, DoubleReal expectedMzShift, DoubleReal maxMzDeviation)
{
	
}

bool SILACFilter::doubleCmp::operator()(DoubleReal a, DoubleReal b) const
{
	//If a m/z position is blacklisted, no SILAC pair may start within the area of 0.8*getPeakWidth(m/z)
	DoubleReal peak_width=0.8*getPeakWidth((a+b)/2);
	if (std::abs(a-b) < peak_width)
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
	return mz_peptide_separations.size();
}

std::vector<DoubleReal> SILACFilter::getPeakPositions()
{
	return peak_positions;
}

/*std::vector<DoubleReal> SILACFilter::getMzPeptideSeparations()
{
	return mz_peptide_separations;
}*/

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
	
Int SILACFilter::getIsotopesPerPeptide()
{
	return isotopes_per_peptide;
}
	
}
