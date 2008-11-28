// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_PEAKALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_PEAKALIGNMENT_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/DSpectrum.h>

#include <cfloat>
#include <cmath>
#include <vector>
#include <map>

namespace OpenMS
{

	/**
		@brief make a PeakAlignment of two PeakSpectra

		The alignment is done according to the Needleman-Wunsch Algorithm (local alignment considering gaps).

		@ref PeakAlignment_Parameters are explained on a separate page.

		@ingroup SpectraComparison
	*/

	class OPENMS_DLLAPI PeakAlignment : public PeakSpectrumCompareFunctor
	{
	public:

		/// default constructor
		PeakAlignment();

		/// copy constructor
		PeakAlignment(const PeakAlignment& source);

		/// destructor
		virtual ~PeakAlignment();

		/// assignment operator
		PeakAlignment& operator = (const PeakAlignment& source);

		/** function call operator, calculates the similarity of the given arguments

			@param spec1 First spectrum given in a binned representation
			@param spec2 Second spectrum ginve in a binned representation
		*/
		double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;

		/// function call operator, calculates self similarity
		double operator () (const PeakSpectrum& spec) const;

	///
	static PeakSpectrumCompareFunctor* create() { return new PeakAlignment(); }

	/// get the identifier for this FactoryProduct
	static const String getProductName()
	{
		return "PeakAlignment";
	}

	/// make alignment and get the traceback
	std::vector< std::pair<UInt,UInt> > getAlignmentTraceback (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;

	private:

		/// calculates the score for aligning two peaks
	double peakPairScore_(double& pos1, double& intens1, double& pos2, double& intens2, const double& sigma) const;


	};

}
#endif //OPENMS_COMPARISON_SPECTRA_PEAKALIGNMENT_H
