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


#ifndef OPENMS_MATH_MISC_AUTOCORRELATIONCALCULATOR_H
#define OPENMS_MATH_MISC_AUTOCORRELATIONCALCULATOR_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
class OPENMS_DLLAPI AutocorrelationCalculator : public ProgressLogger{
private:
	DoubleReal mz_max,mz_min,max_intensity;
	DoubleReal stepwidth_;
	std::vector<DoubleReal> distances_;
public:
	AutocorrelationCalculator(DoubleReal stepwidth);
	virtual ~AutocorrelationCalculator();
	void calculate(const MSExperiment<Peak1D>& input, MSExperiment<Peak1D>& output);
	void analyzeSpectrum(const MSSpectrum<Peak1D>& input, MSSpectrum<Peak1D>& output);
};
}
#endif /* AUTOCORRELATION_H_ */
