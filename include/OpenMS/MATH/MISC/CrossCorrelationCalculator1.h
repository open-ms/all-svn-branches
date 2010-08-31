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

#ifndef OPENMS_MATH_MISC_CrossCorrelationCalculator1_H
#define OPENMS_MATH_MISC_CrossCorrelationCalculator1_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
class OPENMS_DLLAPI CrossCorrelationCalculator1 : public ProgressLogger{
private:
	DoubleReal stepwidth_;
	DoubleReal gauss_mean_;
	DoubleReal gauss_sigma_;
	DoubleReal mz_max,mz_min;
	std::vector<DoubleReal> analyzeSpectrum(const MSSpectrum<Peak1D>& input, MSSpectrum<Peak1D>& output, bool gauss_fitting=false);
public:
	CrossCorrelationCalculator1(DoubleReal stepwidth,DoubleReal gauss_mean, DoubleReal gauss_sigma);
	std::vector<DoubleReal> calculate(const MSExperiment<Peak1D>& input, MSExperiment<Peak1D>& output,Size spectrum_selection_id);
	std::vector<DoubleReal> getExactPositions(const std::vector<DoubleReal>& data,std::vector<DoubleReal> positions, DoubleReal tolerance);
	virtual ~CrossCorrelationCalculator1();
};
}


#endif /* CROSSCORRELATIONCALCULATION_H_ */
