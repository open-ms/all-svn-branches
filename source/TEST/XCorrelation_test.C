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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/XCorrelation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(XCorrelation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

XCorrelation<Peak1D>* ptr = 0;
START_SECTION(XCorrelation())
{
	ptr = new XCorrelation<Peak1D>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~XCorrelation())
{
	delete ptr;
}
END_SECTION

START_SECTION((XCorrelation(const XCorrelation &source)))
{
	ptr = new XCorrelation<Peak1D>();
  Param param = ptr->getParameters();
	param.setValue("peak_tolerance", 0.3);
	param.setValue("parentmass_tolerance", 0.3);
	ptr->setParameters(param);
	XCorrelation<Peak1D> copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((XCorrelation& operator=(const XCorrelation &source)))
{
	XCorrelation<Peak1D> copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((void getXCorrelation(const SpectrumType &s1, const SpectrumType &s2, DoubleReal &best_score1, DoubleReal &best_score2, DoubleReal &best_shift, std::list< std::pair< Size, Size > > &best_matches, bool pm_diff_shift = false) const ))
{
  MSSpectrum<Peak1D> s1;
  Precursor p;
  p.setMZ(252.2);
  p.setIntensity(1);
  s1.getPrecursors().push_back(p);

  DoubleReal mzs_oripeaks[6] = {1.5, 59.5, 70, 119, 130, 178.5};
  DoubleReal ints_oripeaks[6] = {1, 1, 0.75, 1, 0.75, 1};
  DoubleReal mzs_shifts[4] = {0.3, -0.3, 0.60002, -0.60001};

  for(Size i = 0; i < 6; ++i)
	{
		Peak1D t;
		t.setMZ(mzs_oripeaks[i]);
		t.setIntensity(ints_oripeaks[i]);
		s1.push_back(t);
	}

  std::vector< MSSpectrum<Peak1D> > test_specs(5,s1);

  for(Size i = 1; i < test_specs.size(); ++i)
	{
    test_specs[i].getPrecursors().front().setMZ(test_specs[i].getPrecursors().front().getMZ() + mzs_shifts[i-1]);
    for(Size j = 1; j < test_specs[i].size(); ++j)
    {
      test_specs[i][j].setMZ(test_specs[i][j].getMZ() + mzs_shifts[i-1]);
    }
    test_specs[i].sortByPosition();
  }

  //~ std::cout << "s1: " << std::endl;
  //~ for(Size i = 0; i < s1.size(); ++i)
  //~ {
    //~ std::cout << i <<") " << s1[i].getMZ() << " , " << s1[i].getIntensity() << std::endl;
  //~ }
  //~ std::cout << "test_specs[3]: " << std::endl;
  //~ for(Size i = 0; i < test_specs[3].size(); ++i)
  //~ {
    //~ std::cout << i <<") " << test_specs[3][i].getMZ() << " , " << test_specs[3][i].getIntensity() << std::endl;
  //~ }

  //parameter setting see above!

  DoubleReal best_score1, best_score2, best_shift;
  std::list<std::pair<Size,Size> > best_matches;
  ptr->getXCorrelation(test_specs[0], test_specs[0], best_score1 , best_score2, best_shift, best_matches, false);
  //~ std::cout << "best shift: " << best_shift << std::endl;
  for(std::list<std::pair<Size,Size> >::const_iterator it = best_matches.begin(); it !=  best_matches.end(); ++it)
	{
    //~ std::cout << "-> " << it->first << ", " << it->second << std::endl;
    TEST_EQUAL(it->first, it->second)
  }
  TEST_EQUAL(best_matches.size(), 4)
  TEST_EQUAL(best_shift, 0)

  ptr->getXCorrelation(test_specs[0], test_specs[1], best_score1 , best_score2, best_shift, best_matches, false);
  //~ std::cout << "best shift: " << best_shift << std::endl;
  for(std::list<std::pair<Size,Size> >::const_iterator it = best_matches.begin(); it !=  best_matches.end(); ++it)
	{
    //~ std::cout << "-> " << it->first << ", " << it->second << std::endl;
    TEST_EQUAL(it->first, it->second)
  }
  TEST_EQUAL(best_matches.size(), 4)
  TEST_EQUAL(best_shift, 0)

  ptr->getXCorrelation(test_specs[0], test_specs[2], best_score1 , best_score2, best_shift, best_matches, false);
  //~ std::cout << "best shift: " << best_shift << std::endl;
  for(std::list<std::pair<Size,Size> >::const_iterator it = best_matches.begin(); it !=  best_matches.end(); ++it)
	{
    //~ std::cout << "-> " << it->first << ", " << it->second << std::endl;
    TEST_EQUAL(it->first, it->second)
  }
  TEST_EQUAL(best_matches.size(), 4)
  TEST_EQUAL(best_shift, 0)

  ptr->getXCorrelation(test_specs[0], test_specs[3], best_score1 , best_score2, best_shift, best_matches, false);
  //~ std::cout << "best shift: " << best_shift << std::endl;
  for(std::list<std::pair<Size,Size> >::const_iterator it = best_matches.begin(); it !=  best_matches.end(); ++it)
	{
    //~ std::cout << "-> " << it->first << ", " << it->second << std::endl;
    TEST_EQUAL(it->first, it->second)
  }
  TEST_EQUAL(best_matches.size(), 1)
  TEST_EQUAL(best_shift, 0)

  ptr->getXCorrelation(test_specs[0], test_specs[4], best_score1 , best_score2, best_shift, best_matches, false);
  //~ std::cout << "best shift: " << best_shift << std::endl;
  for(std::list<std::pair<Size,Size> >::const_iterator it = best_matches.begin(); it !=  best_matches.end(); ++it)
	{
    //~ std::cout << "-> " << it->first << ", " << it->second << std::endl;
    TEST_EQUAL(it->first, it->second)
  }
  TEST_EQUAL(best_matches.size(), 1)
  TEST_EQUAL(best_shift, 0)

}
END_SECTION

START_SECTION((void maxSparseMatches(const SpectrumType &s1, const SpectrumType &s2, std::vector< std::pair< Size, Size > > &all_matches, std::list< std::pair< Size, Size > > &best_matches) const ))
{
  // TODO
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



