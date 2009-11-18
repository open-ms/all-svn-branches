// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>

using namespace std;

namespace OpenMS
{
	BinnedSharedPeakCount::BinnedSharedPeakCount()
	  : BinnedSpectrumCompareFunctor()
	{
			setName(BinnedSharedPeakCount::getProductName());
			defaults_.setValue("normalized", 1, "is set 1 if the similarity-measurement is normalized to the range [0,1]");
			defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks for which they are supposed to be from different peptides");
			defaultsToParam_();
	}

	BinnedSharedPeakCount::BinnedSharedPeakCount(const BinnedSharedPeakCount& source)
	  : BinnedSpectrumCompareFunctor(source)
	{
	}

	BinnedSharedPeakCount::~BinnedSharedPeakCount()
	{
	}

	BinnedSharedPeakCount& BinnedSharedPeakCount::operator = (const BinnedSharedPeakCount& source)
	{
		if (this != &source)
		{
	  		BinnedSpectrumCompareFunctor::operator = (source);
		}
	  	return *this;
	}

	DoubleReal BinnedSharedPeakCount::operator () (const BinnedSpectrum<>& spec) const
	{
		return operator () (spec, spec);
	}

	DoubleReal BinnedSharedPeakCount::operator () (const BinnedSpectrum<>& spec1, const BinnedSpectrum<>& spec2, const bool lookahead) const
	{
		if(!spec1.checkCompliance(spec2))
		{
			cout << "incompatible" << endl;
			throw BinnedSpectrumCompareFunctor::IncompatibleBinning(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
		}

		// shortcut similarity calculation by comparing PrecursorPeaks (PrecursorPeaks more than delta away from each other are supposed to be from another peptide)
		DoubleReal pre_mz1 = 0.0;
		if (!spec1.getPrecursors().empty())
		{
			int c_1= spec1.getPrecursors().front().getCharge();
			/// @attention if the precursor charge is unknown, i.e. 0 best guess is its doubly charged
			(c_1==0)?c_1=2:c_1=c_1;
			pre_mz1 = (spec1.getPrecursors()[0].getMZ()*c_1 + (c_1-1)*Constants::PROTON_MASS_U);
		}
		DoubleReal pre_mz2 = 0.0;
		if (!spec1.getPrecursors().empty())
		{
			int c_2 = spec2.getPrecursors().front().getCharge();
			/// @attention if the precursor charge is unknown, i.e. 0 best guess is its doubly charged
			(c_2==0)?c_2=2:c_2=c_2;
			pre_mz2 = (spec2.getPrecursors()[0].getMZ()*c_2 + (c_2-1)*Constants::PROTON_MASS_U);
		}
		/// @attention singly charged mass difference each!
		DoubleReal pm_diff = pre_mz2-pre_mz1;

		if(fabs(pm_diff)>(DoubleReal)param_.getValue("precursor_mass_tolerance") and !lookahead)
		{
			return 0;
		}

		DoubleReal score(0), sum(0);
		Size offset(0);
		if(lookahead)
		{
			offset = (floor(fabs(pm_diff)/(DoubleReal)spec1.getBinSize()));
			//~ offset = (ceil(fabs(pm_diff)/(DoubleReal)spec1.getBinSize()));
			//~ offset = fabs(spec2.getBinNumber()-spec1.getBinNumber())+1;
			/* debug std::cout << offset << std::endl; */
		}

		Size denominator(max(spec1.getFilledBinNumber(),spec2.getFilledBinNumber()));

		// all bins at equal position that have both intensity > 0 contribute positively to score
		if(pm_diff<0)
		{
			for (Size i = 0; i+offset < spec1.getBinNumber() and i< spec2.getBinNumber(); ++i)
			{
				if(spec1.getBins()[i+offset]>0 && spec2.getBins()[i]>0)
				{
					sum++;
				}
			}
		}
		else
		{ // pm_diff >= 0
			for (Size i = 0; i+offset < spec2.getBinNumber() and i< spec1.getBinNumber(); ++i)
			{
				/* debug std::cout << i << " - "<< spec1.getBins()[i] << "," << spec2.getBins()[i+offset] << "|"; */
				if(spec1.getBins()[i]>0 && spec2.getBins()[i+offset]>0)
				{
					sum++;
				}
			}
		}
		/* debug std::cout << std::endl << sum << std::endl; */

		// resulting score normalized to interval [0,1]
		score = sum / (DoubleReal) denominator;

		return score;
	}

}
