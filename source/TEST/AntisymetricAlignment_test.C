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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/AntisymetricAlignment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AntisymetricAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AntisymetricAlignment<Peak1D>* ptr = 0;
START_SECTION(AntisymetricAlignment())
{
	ptr = new AntisymetricAlignment<Peak1D>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~AntisymetricAlignment())
{
	delete ptr;
}
END_SECTION

START_SECTION((AntisymetricAlignment(const AntisymetricAlignment &source)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~AntisymetricAlignment()))
{
  // TODO
}
END_SECTION

START_SECTION((AntisymetricAlignment& operator=(const AntisymetricAlignment &source)))
{
  // TODO
}
END_SECTION

START_SECTION((void findMatchingJumps(Size &index, SpectrumType &s1sym, std::vector< Real > &jump_masses, std::vector< int > &leftMatches, std::vector< int > &rightMatches) const ))
{
  // TODO
}
END_SECTION

START_SECTION((MSSpectrum<PeakT> getSymetricSpectrum(const MSSpectrum< PeakT > &spectrum) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void getAntisymetricAlignment(MSSpectrum< PeakT > &res_1, MSSpectrum< PeakT > &res_2, DoubleReal &score, DoubleReal &mod_pos, const MSSpectrum< PeakT > &s1, const MSSpectrum< PeakT > &s2) const ))
{
  DoubleReal score, mod_pos;
  MSSpectrum<Peak1D> res_1, res_2, aacc, aawc;
  Precursor pr;
  pr.setMZ(184.06); // (2*|A|+2*|C|+H2O+2*prot)/2
  pr.setIntensity(1);
  pr.setCharge(2);
  aacc.getPrecursors().push_back(pr);
  aawc.getPrecursors().push_back(pr);

  //~ AACC peptide
  DoubleReal mzs_peaks[6] = {73.06, 121.02, 144.09, 224.03, 247.10, 295.07};
  DoubleReal ints_peaks[6] = {10, 5, 10, 5, 10, 5};
  //~ DoubleReal ints_peaks[6] = {5, 10, 5, 10, 5, 10};
  //~ DoubleReal ints_peaks[6] = {10, 10, 10, 10, 10, 10};

  for(Size i = 0; i < 6; ++i)
	{
		Peak1D t;
		t.setMZ(mzs_peaks[i]);
		t.setIntensity(ints_peaks[i]);
		aacc.push_back(t);
		aawc.push_back(t);
	}
  aawc.getPrecursors().front().setMZ(225.60);
  for(Size i = 4; i < 6; ++i)
	{
		aawc[i].setMZ(aawc[i].getMZ()+83.07);
	}

  AntisymetricAlignment<Peak1D> asa;
  asa.getAntisymetricAlignment(res_1, res_2, score, mod_pos, aacc, aacc);
  TEST_EQUAL(res_1.size(),res_2.size())
  TEST_EQUAL(res_1.size(),3)
  TEST_EQUAL(mod_pos,-1)
  for(Size i = 0; i < res_1.size(); ++i)
  {
    TEST_EQUAL(res_1[i].getMZ(),res_2[i].getMZ())
    TEST_EQUAL(res_1[i].getMZ(), mzs_peaks[i*2])
  }

  asa.getAntisymetricAlignment(res_1, res_2, score, mod_pos, aacc, aawc);
  TEST_EQUAL(res_1.size(),res_2.size())
  TEST_EQUAL(res_1.size(),3)
  TEST_EQUAL(mod_pos,247.1)
  for(Size i = 0; i < res_1.size(); ++i)
  {
    if(res_1[i].getMZ()<mod_pos)
      TEST_EQUAL(res_1[i].getMZ(),res_2[i].getMZ())
    else
      TEST_EQUAL(res_1[i].getMZ()+83.07,res_2[i].getMZ())
  }

  aacc[5].setMZ(aacc[5].getMZ()+103.01-71.03);
  Peak1D t;
	t.setMZ(350.11);
	t.setIntensity(ints_peaks[10]);
  aacc.push_back(t);
	t.setMZ(398.08);
	t.setIntensity(ints_peaks[5]);
  aacc.push_back(t);
  aacc.getPrecursors().front().setMZ(235.57); // aacc is now aaccc

  aawc[5].setMZ(aawc[5].getMZ()+103.01-71.03);
	t.setMZ(433.18);
 	t.setIntensity(ints_peaks[10]);
  aawc.push_back(t);
	t.setMZ(481.15);
	t.setIntensity(ints_peaks[5]);
  aawc.push_back(t);
  aawc.getPrecursors().front().setMZ(277.10);

  for(Size i = 0; i < aacc.size(); ++i)
  {
    std::cout << aacc[i].getMZ() << " vs " << aawc[i].getMZ() << " vs " << (aacc[i].getMZ()+83.07) << std::endl;
  }

  asa.getAntisymetricAlignment(res_1, res_2, score, mod_pos, aacc, aawc);
  TEST_EQUAL(res_1.size(),res_2.size())
  TEST_EQUAL(res_1.size(),4)
  TEST_EQUAL(mod_pos,247.1)
  for(Size i = 0; i < res_1.size(); ++i)
  {
    if(res_1[i].getMZ()<mod_pos)
      TEST_EQUAL(res_1[i].getMZ(),res_2[i].getMZ())
    else
      TEST_EQUAL(res_1[i].getMZ()+83.07,res_2[i].getMZ())
  }


}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



