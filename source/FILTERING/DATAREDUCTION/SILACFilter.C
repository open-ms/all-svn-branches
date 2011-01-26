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
// $Authors: Lars Nilse, Steffen Sass, Holger Plattfaut $
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
		numberOfPeptides = mass_separations.size();    // number of labelled peptides +1 [e.g. for SILAC triplet =3]
		isotopes_per_peptide=isotopes_per_peptide_;    // isotopic peaks per peptide
		charge=charge_;    // peptide charge
		model_deviation=model_deviation_;    // allowed deviation from averegine model
		isotope_distance=1.000495/(DoubleReal) charge;    // distance between isotopic peaks of a peptide [Th]
		
		// m/z shifts from mass shifts
		mz_peptide_separations.push_back(0.0);
		for (std::set<DoubleReal>::iterator it = mass_separations.begin(); it!=mass_separations.end(); ++it)
		{
			mz_peptide_separations.push_back(*it/(DoubleReal)charge);
		}
		
		silac_type=mass_separations.size(); // rausschmeissen??
	}

SILACFilter::~SILACFilter() {
}

bool SILACFilter::isSILACPattern(DoubleReal rt, DoubleReal mz)
{
	current_mz = mz;
	
	//---------------------------------------------------------------
	// BLUNT INTENSITY FILTER (Just check that intensity at current position is above the intensity cutoff.)
	//---------------------------------------------------------------
	if (gsl_spline_eval (SILACFiltering::spline_lin, mz, SILACFiltering::current_lin) < SILACFiltering::intensity_cutoff)
	{
		return false;
	}

	
	//---------------------------------------------------------------
	// EXACT m/z SHIFTS (Determine the actual shifts between peaks. Say 4 Th is the theoretic shift. In the experimental data it will be 4.0029 Th.)
	//---------------------------------------------------------------
	exact_shifts.clear();
	exact_intensities.clear();
	//std::cout << " m/z = " << mz << " RT = " << rt << std::endl;
	for (Int peptide = 0; peptide <= numberOfPeptides; peptide++) // loop over labelled peptides [e.g. for SILAC triplets: 0=light 1=medium 2=heavy]
	{
		std::vector<DoubleReal> exact_shifts_singlePeptide;
		std::vector<DoubleReal> exact_intensities_singlePeptide;
		for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++) // loop over isotopic peaks within a peptide [0=mono-isotopic peak etc.]
		{
			DoubleReal deltaMZ = computeActualMzShift(mz, mz_peptide_separations[peptide] + isotope*isotope_distance, 0.001);
			exact_shifts_singlePeptide.push_back( deltaMZ );
			exact_intensities_singlePeptide.push_back( gsl_spline_eval (SILACFiltering::spline_spl, mz + deltaMZ, SILACFiltering::current_spl) );
			//std::cout << "   " << mz_peptide_separations[peptide] + isotope*isotope_distance << "  (" << deltaMZ << ")";
		}
		exact_shifts.push_back(exact_shifts_singlePeptide);
		exact_intensities.push_back(exact_intensities_singlePeptide);
		//std::cout << std::endl;
	}
	//std::cout << std::endl << std::endl;
	
	
	//---------------------------------------------------------------
	// COMPLETE INTENSITY FILTER (Check that all of the intensities are above the cutoff.)
	//---------------------------------------------------------------
	for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
	{
		for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++)
		{
			if (exact_shifts[peptide][isotope] == -1) // no correlating signal was found
			{
				return false;
			}
			else
			{
				if (exact_intensities[peptide][isotope] < SILACFiltering::intensity_cutoff)
				{
					return false;    // If only one intensity is below the cutoff, return 'false'.
				}
			}
		}
	}
	
	
	//---------------------------------------------------------------
	// CORRELATION FILTER (Check that for each peptide all possible combinations of isotopic peaks correlate.)
	//---------------------------------------------------------------
	for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
	{
		for (Int isotope1 = 0; isotope1 < isotopes_per_peptide; isotope1++)
		{
			for (Int isotope2 = 0; isotope2 < isotopes_per_peptide; isotope2++)
			{
				if (isotope1 != isotope2)
				{
					std::vector<DoubleReal> intensities1;    // intensities in region around isotopic peak 1
					std::vector<DoubleReal> intensities2;    // intensities in region around isotopic peak 2
					DoubleReal mzWindow = 0.7 * getPeakWidth(mz);    // width of the window around m/z in which the correlation is calculated
					for (DoubleReal dmz = - mzWindow; dmz <= mzWindow; dmz += 0.2 * mzWindow)     // fill intensity vectors
					{
						DoubleReal intens1 = gsl_spline_eval(SILACFiltering::spline_spl, mz + exact_shifts[peptide][isotope1] + dmz, SILACFiltering::current_spl);
						DoubleReal intens2 = gsl_spline_eval(SILACFiltering::spline_spl, mz + exact_shifts[peptide][isotope2] + dmz, SILACFiltering::current_spl);
						intensities1.push_back( intens1 );
						intensities2.push_back( intens2 );
					}
					DoubleReal intensityCorrelation = Math::pearsonCorrelationCoefficient( intensities1.begin(), intensities1.end(), intensities2.begin(), intensities2.end());    // calculate Pearson correlation coefficient
					if ( intensityCorrelation < 0.95 ) return false;
				}
			}
			
		}
	}
	
	
	//---------------------------------------------------------------
	// AVERAGINE FILTER (Check if realtive ratios confirm with an averagine model of all peptides.)
	//---------------------------------------------------------------
	ExtendedIsotopeModel averagineModel;
	Param avaregineParam;
	avaregineParam.setValue( "isotope:monoisotopic_mz", mz );
	avaregineParam.setValue("charge", charge);
	avaregineParam.setValue("isotope:stdev", getPeakWidth(mz)*0.42466);    // Where does the 0.42466 come from?
	averagineModel.setParameters( avaregineParam );
	std::vector<Peak1D> model_data;
	averagineModel.getSamples(model_data);
	// ...
	
	
	//---------------------------------------------------------------
	// ALL FILTERS PASSED => CREATE DATAPOINT
	//---------------------------------------------------------------
	DataPoint newElement;    // Raw data point at this particular RT and m/z passed all filters. Store it for further clustering.
	newElement.feature_id = SILACFiltering::feature_id;
	newElement.rt = rt;
	newElement.mz = mz;
	newElement.charge = charge;
	newElement.intensities.insert(newElement.intensities.begin(), exact_intensities.begin(), exact_intensities.end());
	newElement.mass_shifts.insert(newElement.mass_shifts.begin(), mz_peptide_separations.begin(), mz_peptide_separations.end());
	elements.push_back(newElement);
	
	return true;
}

	
	
	
	
	
	
	
	
	
bool SILACFilter::checkPattern(DoubleReal mz, const std::vector<DoubleReal>& exact_shifts, std::vector<DoubleReal>& intensities,bool missing_peak)
{
	std::vector<DoubleReal> act_intensities;
	//Store all intensities of the current isotope pattern
	for (std::vector<DoubleReal>::const_iterator position_it=exact_shifts.begin();position_it!=exact_shifts.end();++position_it)
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
			DoubleReal intensity2=gsl_spline_eval (SILACFiltering::spline_spl, pos+exact_shifts[i], SILACFiltering::current_spl);
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
	param.setValue( "isotope:monoisotopic_mz", mz+exact_shifts[0] );
	param.setValue("charge",charge);
	param.setValue("isotope:stdev",getPeakWidth(mz+exact_shifts[0])*0.42466);
	model.setParameters( param );
	std::vector<Peak1D> model_data;
	model.getSamples(model_data);

	for (Size i=0;i< 2;++i)
	{
		DoubleReal quality=(act_intensities[i]*model.getIntensity(mz+exact_shifts[i+1]))/(act_intensities[i+1]*model.getIntensity(mz+exact_shifts[i]));
		if ((std::abs(log(quality)) > model_deviation && i!=1) || (std::abs(log(quality)) > model_deviation && missing_peak && i==1))
			return false;
	}
	intensities.insert(intensities.end(),act_intensities.begin(),act_intensities.end());
	return true;
}

	
	

	
	
	
	
	

DoubleReal SILACFilter::computeActualMzShift(DoubleReal mz, DoubleReal expectedMzShift, DoubleReal maxMzDeviation)
{
	if (expectedMzShift <= 0.0)
	{
		return 0;
	}
	else
	{
		std::vector<DoubleReal> tempVector;
		DoubleReal stepwidth = maxMzDeviation/18;
		//n must be a power of two for gsl fast fourier transformation; take next higher size for n, which is a power of two and fill the rest with zeros
		Size tempVector_size = pow(2,(ceil(log2((3*maxMzDeviation)/stepwidth))));
		
		tempVector.clear();
		tempVector.resize(tempVector_size);
		
		//Create a vector containing interpolated values with a spacing of 'stepwidth'
		Size i=0;
		DoubleReal starting_offset=std::min(mz - SILACFiltering::mz_min, maxMzDeviation);  // Why no mz_max?
		for (DoubleReal x = mz - starting_offset; x<=mz + maxMzDeviation; x+=stepwidth)
		{
			tempVector[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x, SILACFiltering::current_lin);
			++i;
		}
		
		std::vector<DoubleReal> gauss_fitted_data(tempVector_size,0);
		i=0;
		for (DoubleReal xx=mz-starting_offset-maxMzDeviation;xx<=mz+maxMzDeviation;xx+=stepwidth)
		{
			gauss_fitted_data[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, xx+expectedMzShift, SILACFiltering::current_lin);
			++i;
		}
		
		//Fourier transformation of the values
		gsl_fft_real_radix2_transform (&(*tempVector.begin()), 1, tempVector_size);
		gsl_fft_real_radix2_transform (&(*gauss_fitted_data.begin()), 1, tempVector_size);
		
		//Multiply the fourier transformed complex values with the complex conjugate
		//Have a look at the GNU Scientific Library reference manual for a description of the data structure.
		tempVector[0]=tempVector[0]*gauss_fitted_data[0];
		tempVector[tempVector_size/2]=tempVector[tempVector_size/2]*gauss_fitted_data[tempVector_size/2];
		for (i = 1; i <= tempVector_size/2; ++i)
		{
			tempVector[i]=tempVector[i]*gauss_fitted_data[i] + tempVector[tempVector_size-i]*gauss_fitted_data[tempVector_size-i];
			tempVector[tempVector_size-i]=0.0;
		}
		
		//Compute inverse fourier transformation
		gsl_fft_halfcomplex_radix2_inverse (&(*tempVector.begin()), 1, tempVector_size);
		
		tempVector.resize(tempVector_size/2);
		
		//Interpolate the current autocorrelation peak
		std::vector<DoubleReal> x(tempVector.size(),0);
		std::vector<DoubleReal> y(tempVector.size(),0);
		DoubleReal shift=expectedMzShift - maxMzDeviation;
		for (Size i=0;i<tempVector.size();++i)
		{
			x[i]=shift;
			y[i]=tempVector[i];
			shift+=stepwidth;
		}
		
		gsl_interp_accel* acc_correlation = gsl_interp_accel_alloc();
		gsl_spline* spline_correlation = gsl_spline_alloc(gsl_interp_cspline, x.size());
		gsl_spline_init(spline_correlation, &(*x.begin()), &(*y.begin()), x.size());
		
		//Search for first maximum in + and - direction and check preconditions
		for (DoubleReal act_position=0.0;act_position<=maxMzDeviation; act_position+=stepwidth)
		{
			DoubleReal last_intensity=gsl_spline_eval (spline_correlation, expectedMzShift+act_position-stepwidth, acc_correlation);
			DoubleReal act_intensity=gsl_spline_eval (spline_correlation, expectedMzShift+act_position, acc_correlation);
			DoubleReal next_intensity=gsl_spline_eval (spline_correlation, expectedMzShift+act_position+stepwidth, acc_correlation);
			// search for an actual m/z shift larger than the expected one
			if (act_intensity > last_intensity && act_intensity > next_intensity && act_intensity > 1000 && gsl_spline_eval (SILACFiltering::spline_lin, mz+expectedMzShift+act_position, SILACFiltering::current_lin) > SILACFiltering::intensity_cutoff) // Why fixed intensity cutoffs?
			{
				gsl_spline_free(spline_correlation);
				gsl_interp_accel_free(acc_correlation);
				return expectedMzShift+act_position;
			}
			last_intensity=gsl_spline_eval(spline_correlation, expectedMzShift-act_position-stepwidth, acc_correlation);
			act_intensity=gsl_spline_eval(spline_correlation, expectedMzShift-act_position, acc_correlation);
			next_intensity=gsl_spline_eval(spline_correlation, expectedMzShift-act_position+stepwidth, acc_correlation);
			// search for an actual m/z shift smaller than the expected one
			if (act_intensity > last_intensity && act_intensity > next_intensity && act_intensity > 1000 && gsl_spline_eval (SILACFiltering::spline_lin, mz+expectedMzShift-act_position, SILACFiltering::current_lin) > SILACFiltering::intensity_cutoff)
			{
				gsl_spline_free(spline_correlation);
				gsl_interp_accel_free(acc_correlation);
				return expectedMzShift-act_position;
			}
		}
		gsl_spline_free(spline_correlation);
		gsl_interp_accel_free(acc_correlation);
		return -1;		
	}
	
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
	std::vector<DoubleReal> exact_positions;    // exact m/z positions of isotopic peaks as a flat vector 
	for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
	{
		for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++)
		{
			exact_positions.push_back( current_mz + exact_shifts[peptide][isotope] );
		}
	}
	return exact_positions;
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
	
Int SILACFilter::getIsotopesPerPeptide()
{
	return isotopes_per_peptide;
}
	
}
