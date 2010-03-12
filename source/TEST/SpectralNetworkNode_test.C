// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/SpectralNetworkNode.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectralNetworkNode, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectralNetworkNode<Peak1D>* ptr = 0;
START_SECTION(SpectralNetworkNode())
{
	ptr = new SpectralNetworkNode<Peak1D>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~SpectralNetworkNode())
{
	delete ptr;
}
END_SECTION

START_SECTION((SpectralNetworkNode(const SpectralNetworkNode &source)))
{
	ptr = new SpectralNetworkNode<Peak1D>();
  Param param = ptr->getParameters();
	param.setValue("peak_tolerance", 0.3);
	param.setValue("parentmass_tolerance", 0.3);
	ptr->setParameters(param);
	SpectralNetworkNode<Peak1D> copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((void reverseSpectrum(MSSpectrum< PeakType > &spec, bool is_prm=false)))
{
  MSSpectrum<Peak1D> s1,r1;
  Precursor p;
  p.setMZ(252.2);
  p.setIntensity(1);
  p.setCharge(1);
  s1.getPrecursors().push_back(p);

  DoubleReal mzs_oripeaks[6] = {1.5, 59.5, 70, 119, 130, 178.5};
  DoubleReal mzs_revpeaks[6] = {74.7073, 123.207, 134.207, 183.207, 193.707, 251.707};
  DoubleReal ints_oripeaks[6] = {1, 1, 0.75, 1, 0.75, 1};
  for(Size i = 0; i < 6; ++i)
	{
		Peak1D t;
		t.setMZ(mzs_oripeaks[i]);
		t.setIntensity(ints_oripeaks[i]);
		s1.push_back(t);
	}
  for(Size i = 0; i < 6; ++i)
	{
		Peak1D t;
		t.setMZ(mzs_revpeaks[i]);
		t.setIntensity(ints_oripeaks[i]);
		r1.push_back(t);
	}

  ptr->reverseSpectrum(s1);
  TEST_EQUAL(s1.size(), r1.size())
  for(Size i = 0; i < r1.size() && i < s1.size(); ++i)
	{
    TEST_REAL_SIMILAR(s1[i].getMZ(), r1[i].getMZ());
  }
}
END_SECTION

START_SECTION((void getPropagation(AASequence &template_sequence, DoubleReal mod_pos, DoubleReal alt_mod_pos, DoubleReal mod_shift, AASequence &modified_sequence, String &sequence_correspondence)))
{
  //~ DF([3])G masses: {116.070, 117.042, 266.103};
  AASequence modseq("DF([3])PANGIER"), aaseq("DFPANGIER"), res;
  DoubleReal p1(263.103), p2(266.103), s(3);
  String bla;
  ptr->getPropagation(aaseq,p1,p2,s,res,bla);
  TEST_EQUAL(modseq, res);
}
END_SECTION

START_SECTION((void spanNetwork(MSExperiment< PeakType > &consensuses, ConsensusMap &id_consensuses, std::vector< Size > &indices_id_consensuses, std::vector< std::pair< Size, Size > > &aligned_pairs, std::vector< DoubleReal > &mod_positions, Size hops)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void createNetwork(std::vector< std::pair< Size, Size > > &aligned_pairs, MSExperiment< PeakType > &consensuses, std::map< Size, std::set< Size > > &stars)))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



