// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/RTProbability.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RTProbability, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RTProbability* ptr = 0;
RTProbability* null_ptr = 0;
START_SECTION(RTProbability())
{
	ptr = new RTProbability();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~RTProbability())
{
	delete ptr;
}
END_SECTION

ptr = new RTProbability();
START_SECTION((RTProbability(const RTProbability &rhs)))
{
  ptr->setGaussianParameters(0.5,10.);
  RTProbability rtprob2(*ptr);
  TEST_REAL_SIMILAR(rtprob2.getGaussMu(),0.5)
  TEST_REAL_SIMILAR(rtprob2.getGaussSigma(),10.)   
}
END_SECTION


START_SECTION((RTProbability& operator=(const RTProbability &rhs)))
{
  ptr->setGaussianParameters(0.5,12.);
  RTProbability rtprob2 = *ptr;
  TEST_REAL_SIMILAR(rtprob2.getGaussMu(),0.5)
  TEST_REAL_SIMILAR(rtprob2.getGaussSigma(),12.)   
}
END_SECTION

START_SECTION((void learnGaussian(FeatureMap<> &features, String rt_model_path, DoubleReal min_score)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DoubleReal getRTProbability(DoubleReal min_obs_rt, DoubleReal max_obs_rt, DoubleReal pred_rt)))
{
  ptr->setGaussianParameters(0.,40.);
  ptr->getRTProbability(2320.5,3630.,3353.);
}
END_SECTION

START_SECTION((void setGaussianParameters(DoubleReal sigma, DoubleReal mu)))
{
  ptr->setGaussianParameters(2.5,12.9);
  TEST_REAL_SIMILAR(ptr->getGaussMu(),2.5)
  TEST_REAL_SIMILAR(ptr->getGaussSigma(),12.9)   
}
END_SECTION

START_SECTION((DoubleReal getGaussSigma()))
{
  ptr->setGaussianParameters(6.5,19.9);
  TEST_REAL_SIMILAR(ptr->getGaussSigma(),19.9)
}
END_SECTION

START_SECTION((DoubleReal getGaussMu()))
{
  ptr->setGaussianParameters(7.5,19.9);
  TEST_REAL_SIMILAR(ptr->getGaussMu(),7.5)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



