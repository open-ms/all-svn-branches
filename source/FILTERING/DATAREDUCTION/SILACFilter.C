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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
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
	if (isotopes_per_peptide > 1)
	{
		for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
		{

			IsotopeDistribution isoDistribution;    // isotope distribution of an averagine peptide
			isoDistribution.estimateFromPeptideWeight((mz + exact_shifts[peptide][0])*charge);    // mass of averagine peptide
			DoubleReal averagineIntensity_mono = isoDistribution.getContainer()[0].second;    // intensity of monoisotopic peak of the averagine model
			DoubleReal intensity_mono = exact_intensities[peptide][0];    // intensity around the (potential) monoisotopic peak in the real data
			for (Int isotope = 1; isotope < isotopes_per_peptide; isotope++)
			{
				DoubleReal averagineIntensity = isoDistribution.getContainer()[isotope].second;
				DoubleReal intensity = exact_intensities[peptide][isotope];
				std::cout << "theory = " << averagineIntensity /averagineIntensity_mono << "   real = " << intensity /intensity_mono << "   real/theory = " << (intensity /intensity_mono)/(averagineIntensity /averagineIntensity_mono) << "              ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
	
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

	
DoubleReal SILACFilter::computeActualMzShift(DoubleReal mz, DoubleReal expectedMzShift, DoubleReal maxMzDeviation)
{
	if (expectedMzShift <= 0.0)
	{
		return 0;
	}
	else
	{
		//--------------------------------------------------
		// compute autocorrelation
		//--------------------------------------------------
		
		std::vector<DoubleReal> akimaMz;
		DoubleReal stepwidth = maxMzDeviation / 18;
		
		// n must be a power of two for gsl fast fourier transformation; take next higher size for n, which is a power of two and fill the rest with zeros
		Size akimaMz_size = pow(2,(ceil(log2((3*maxMzDeviation)/stepwidth))));
			
		akimaMz.clear();
		akimaMz.resize(akimaMz_size);
			
		// check to not leave experiment
		DoubleReal starting_offset = std::min(mz - SILACFiltering::mz_min, maxMzDeviation);
			
		// calculate akima interpolation for region around mz (+- maxMzDeviation) and store in vector akimaMz
		// starting position: mz - maxMzDeviation
		// ending position: mz + maxMzDeviation
		// 19 steps, stepwidth = maxMzDeviation / 18
		// akima interpolation for position x
		Size i=0;
		for (DoubleReal x = mz - starting_offset; x <= mz + maxMzDeviation; x += stepwidth)
		{
			akimaMz[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x, SILACFiltering::current_lin);
			++i;
		}
		
		// calculate akima interpolation for region around mz + expectedMzShift (- 2 * maxMzDeviation, + maxMzDeviation) and store in vector akimaMzShift
		// starting position: mz - 2 * maxMzDeviation
		// ending position: mz + maxMzDeviation
		// 19 steps, maxMzDeviation / 18 each
		// akima interpolation for position (x + expectedMzShift)
		std::vector<DoubleReal> akimaMzShift(akimaMz_size, 0);
		i=0;
		for (DoubleReal x = mz - starting_offset - maxMzDeviation; x <= mz + maxMzDeviation; x += stepwidth)
		{
			akimaMzShift[i] = gsl_spline_eval_deriv2 (SILACFiltering::spline_lin, x + expectedMzShift, SILACFiltering::current_lin);
			++i;
		}
			
		// fourier transformation of the values (see Master Thesis Steffen page 28 formula 1)
		gsl_fft_real_radix2_transform (&(*akimaMz.begin()), 1, akimaMz_size);
		gsl_fft_real_radix2_transform (&(*akimaMzShift.begin()), 1, akimaMz_size);
			
		// create vector "autoCorrelations" to store autocorrelations
		std::vector<DoubleReal> autoCorrelations;
		autoCorrelations.resize(akimaMz_size);
			
		// multiply the fourier transformed complex values with the complex conjugate (see Master Thesis Steffen page 28 formula 2)
		// have a look at the GNU Scientific Library reference manual for a description of the data structure.
		autoCorrelations[0] = akimaMz[0] * akimaMzShift[0];     // special case for i = 0
			
		autoCorrelations[akimaMz_size / 2] = akimaMz[akimaMz_size / 2] * akimaMzShift[akimaMz_size / 2];      // special case for i = akimaMz_size / 2
			
		for (i = 1; i <= akimaMz_size / 2; ++i)
		{
			autoCorrelations[i] = akimaMz[i] * akimaMzShift[i] + akimaMz[akimaMz_size - i] * akimaMzShift[akimaMz_size - i];      // cases for (0 < i < akimaMz_size/ 2)
			autoCorrelations[akimaMz_size - i]=0.0;     // for (akimaMz_size / 2 < i akimaMz_size) fill the vector "autoCorrelations" with zeros
		}
			
		// compute inverse fourier transformation (see Master Thesis Steffen page 28 formula 3)
		gsl_fft_halfcomplex_radix2_inverse (&(*autoCorrelations.begin()), 1, akimaMz_size);
			
		autoCorrelations.resize(akimaMz_size / 2);     // cut the vector in half to erase zero entries
			
			
		//--------------------------------------------------
		// compute exact position
		//--------------------------------------------------
			
		// create two vectors for interpolation
		std::vector<DoubleReal> mzShifts(autoCorrelations.size(), 0);      // vector "mzShifts" contains potetial mz shifts (starting by expectedMzShift - maxMzDeviation, moving by maxMzDeviation / 18, ending by the end of vector "autocorrelations")
		std::vector<DoubleReal> correspondingAutoCorrelations(autoCorrelations.size(), 0);      // vector "correspondingAutoCorrelations" contains autocorrelations that corresponds to potetial mz shifht in vector "mzShifts"
			
		DoubleReal shift = expectedMzShift - maxMzDeviation;      // caculate first potential mz shift
			
		// fill vectors "mzShifts" and "correspondingAutoCorrelations"
		for (Size i = 0; i < autoCorrelations.size(); ++i)
		{
			mzShifts[i] = shift;
			correspondingAutoCorrelations[i] = autoCorrelations[i];
			shift += stepwidth;
		}
			
		// interpolate the current autocorrelation peak
		gsl_interp_accel* acc_correlation = gsl_interp_accel_alloc();
		gsl_spline* spline_correlation = gsl_spline_alloc(gsl_interp_cspline, mzShifts.size());
		gsl_spline_init(spline_correlation, &(*mzShifts.begin()), &(*correspondingAutoCorrelations.begin()), mzShifts.size());
			
			
		for (DoubleReal current_position = 0.0; current_position <= maxMzDeviation; current_position += stepwidth)
		{
			// search for first maximum in + direction and check preconditions
			DoubleReal last_intensity = gsl_spline_eval (spline_correlation, expectedMzShift + current_position - stepwidth, acc_correlation);
			DoubleReal current_intensity = gsl_spline_eval (spline_correlation, expectedMzShift + current_position, acc_correlation);
			DoubleReal next_intensity = gsl_spline_eval (spline_correlation, expectedMzShift + current_position + stepwidth, acc_correlation);
			
			// search for a current m/z shift larger than the expected one
			// conditions are: current intensity > 1000 (current intensity calculated with cubic interpolation based on autocorrelation)
			// intensity at position (mz + expectedMzShift + current_position) > intesity_cutoff (intensity calculated with akima interpolation based on "intensities_vec" from SILACFiltering)
			if (current_intensity > last_intensity && current_intensity > next_intensity && gsl_spline_eval (SILACFiltering::spline_lin, mz + expectedMzShift + current_position, SILACFiltering::current_lin) > SILACFiltering::intensity_cutoff) // Why fixed intensity cutoffs?
			{
				gsl_spline_free(spline_correlation);      // free interpolation object
				gsl_interp_accel_free(acc_correlation);     // free accelerator object
				return expectedMzShift + current_position;      // return exact position
			}
			
			// search for first maximum in - direction and check preconditions
			last_intensity = gsl_spline_eval(spline_correlation, expectedMzShift - current_position - stepwidth, acc_correlation);
			current_intensity = gsl_spline_eval(spline_correlation, expectedMzShift - current_position, acc_correlation);
			next_intensity = gsl_spline_eval(spline_correlation, expectedMzShift - current_position + stepwidth, acc_correlation);
			
			// search for an current m/z shift smaller than the expected one
			// conditions are: current intensity > 1000 (current intensity calculated with cubic interpolation based on autocorrelation)
			// intensity at position (mz + expectedMzShift - current_position) > intesity_cutoff (intensity calculated with akima interpolation based on "intensities_vec" from SILACFiltering)
			if (current_intensity > last_intensity && current_intensity > next_intensity && gsl_spline_eval (SILACFiltering::spline_lin, mz + expectedMzShift - current_position, SILACFiltering::current_lin) > SILACFiltering::intensity_cutoff)
			{
				gsl_spline_free(spline_correlation);      // free interpolation object
				gsl_interp_accel_free(acc_correlation);     // free accelerator object
				return expectedMzShift - current_position;      // return exact position
			}
		}
		gsl_spline_free(spline_correlation);      // free interpolation object
		gsl_interp_accel_free(acc_correlation);     // free accelerator object
		return -1;      // return -1 if no autocorrelation exists for expectedMzShift
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
