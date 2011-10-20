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
#include <OpenMS/ANALYSIS/TARGETED/PSProteinInference.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PSProteinInference, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PSProteinInference* ptr = 0;
PSProteinInference* null_ptr = 0;
START_SECTION(PSProteinInference())
{
	ptr = new PSProteinInference();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PSProteinInference())
{
	delete ptr;
}
END_SECTION

vector<String> accs;
vector<PeptideIdentification> ids;
PeptideHit hit;
PeptideIdentification id;
id.setHigherScoreBetter(true);
vector<PeptideHit> hits;
  
hit.setSequence(AASequence("PEPNE"));
hit.setScore(0.7);
accs.push_back("prot1");
hit.setProteinAccessions(accs);
id.insertHit(hit);
ids.push_back(id);

hit.setSequence(AASequence("PEPTW"));
hit.setScore(0.75);
accs.clear();  accs.push_back("prot2");
hit.setProteinAccessions(accs);
hits.push_back(hit);
id.setHits(hits);
ids.push_back(id);

hit.setSequence(AASequence("PEPTW"));
hit.setScore(0.9);
accs.clear();  accs.push_back("prot2");
hit.setProteinAccessions(accs);
hits.clear();
hits.push_back(hit);
id.setHits(hits);
ids.push_back(id);

hit.setSequence(AASequence("PEPTHREE"));
hit.setScore(0.5);
accs.clear();   accs.push_back("prot1");   accs.push_back("prot2"); accs.push_back("prot3");  accs.push_back("prot4");
hit.setProteinAccessions(accs);
hits.clear();hits.push_back(hit);
id.setHits(hits);
ids.push_back(id);

hit.setSequence(AASequence("PEPFUR"));
hit.setScore(0.8);
accs.clear();  accs.push_back("prot3");
hit.setProteinAccessions(accs);
hits.clear();hits.push_back(hit);
id.setHits(hits);
ids.push_back(id);

PSProteinInference pi;
START_SECTION((Size findMinimalProteinList(std::vector< String > &accessions, const std::vector< PeptideIdentification > &peptide_ids)))
{
  TEST_EQUAL(pi.findMinimalProteinList(accs,ids),3)
}
END_SECTION

std::vector<DoubleReal> probabilities;
START_SECTION((void calculateProteinProbabilities(const std::vector< String > &accessions, const std::vector< PeptideIdentification > &peptide_ids,
                                                  std::vector< DoubleReal > &probabilities)))
{
  pi.calculateProteinProbabilities(accs,ids,probabilities);
  TEST_REAL_SIMILAR(probabilities[0],0.85)
  TEST_EQUAL(accs[0],"prot1")
  TEST_REAL_SIMILAR(probabilities[1],0.95)
  TEST_EQUAL(accs[1],"prot2")
  TEST_REAL_SIMILAR(probabilities[2],0.9)
  TEST_EQUAL(accs[2],"prot3")


  for(Size i = 0; i < ids.size();++i)  ids[i].setHigherScoreBetter(false);
  accs.push_back(String("prot5"));
  pi.calculateProteinProbabilities(accs,ids,probabilities);
  TEST_REAL_SIMILAR(probabilities[0],0.65)
  TEST_REAL_SIMILAR(probabilities[1],0.625)
  TEST_REAL_SIMILAR(probabilities[2],0.6)
  TEST_REAL_SIMILAR(probabilities[3],0.)
}
END_SECTION


START_SECTION((DoubleReal getProteinProbability(const String & acc, const std::vector< String > &accessions, const std::vector< DoubleReal > &probabilities)))
{
  TEST_REAL_SIMILAR(pi.getProteinProbability("prot1",accs,probabilities),0.65)  
}
END_SECTION

START_SECTION((DoubleReal getProteinProbability(const String & acc)))
{
  String acc("prot1");
  TEST_REAL_SIMILAR(pi.getProteinProbability(acc),0.65)  
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



