// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Jens Allmer$
// $Authors: Jens Allmer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Chromosome, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Chromosome* ptr = new Chromosome(AASequence("TESTER"),1,0.5);;
Chromosome* null_ptr = 0;
Chromosome* ptt = 0;
START_SECTION(Chromosome())
{
	ptt = new Chromosome();
	TEST_NOT_EQUAL(ptt, null_ptr)
}
END_SECTION

START_SECTION((Chromosome(const AASequence &seq, const double val)))
{
  AASequence aas("TESTER");
  double val = 1.0;
  Chromosome chr(aas,val);
  TEST_EQUAL(chr.getSequence().toString(),"TESTER");
}
END_SECTION

START_SECTION((Chromosome(const Chromosome &other)))
{
  Chromosome chrT(AASequence("TESTER"),1,0.5);
  Chromosome chrA(chrT);
  TEST_EQUAL(chrA.getSequence().toString(),"TESTER");
  TEST_EQUAL(chrA.getCharge(),1);
  TEST_REAL_SIMILAR(chrA.getScore(),0.5);
  TEST_EQUAL(chrA.isScored(),true);
}
END_SECTION

START_SECTION((Chromosome& operator=(const Chromosome &rhs)))
{
  Chromosome chrT(AASequence("TESTER"),1,0.5);
  Chromosome chrA = chrT;
  TEST_EQUAL(chrA.getSequence().toString(),"TESTER");
  TEST_EQUAL(chrA.getCharge(),1);
  TEST_REAL_SIMILAR(chrA.getScore(),0.5);
  TEST_EQUAL(chrA.isScored(),true);
}
END_SECTION

START_SECTION((double getScore() const ))
{
  TEST_REAL_SIMILAR(ptr->getScore(),0.5);
}
END_SECTION

START_SECTION((void setScore(double score)))
{
  ptr->setScore(1.0);
  TEST_REAL_SIMILAR(ptr->getScore(),1.0);
}
END_SECTION

START_SECTION((const AASequence& getSequence() const ))
{
  TEST_STRING_EQUAL(ptr->getSequence().toString(),"TESTER");
}
END_SECTION

START_SECTION((void setSequence(const AASequence &sequence)))
{
  ptr->setSequence(AASequence("ALLMER"));
  TEST_STRING_EQUAL(ptr->getSequence().toString(),"ALLMER");
}
END_SECTION

START_SECTION((static bool sortScoreDesc(boost::shared_ptr<Chromosome> lhs, cboost::shared_ptr<Chromosome> rhs)))
{
  boost::shared_ptr<Chromosome> chrT(new Chromosome(AASequence("TESTER"),1,0.5));
  boost::shared_ptr<Chromosome> chrA(new Chromosome(AASequence("ALLMER"),1,1.0));
  TEST_EQUAL(Chromosome::sortScoreDesc(chrT, chrA),false);
}
START_SECTION(~Chromosome())
{
	delete ptr;
}
END_SECTION
END_SECTION

START_SECTION((static bool sortScoreAsc(boost::shared_ptr<Chromosome> lhs, boost::shared_ptr<Chromosome> rhs)))
{
  boost::shared_ptr<Chromosome> chrT(new Chromosome(AASequence("TESTER"),1,0.5));
  boost::shared_ptr<Chromosome> chrA(new Chromosome(AASequence("ALLMER"),1,1.0));
  TEST_EQUAL(Chromosome::sortScoreAsc(chrT,chrA),true);
}
END_SECTION

START_SECTION((void sort()))
{
  vector<boost::shared_ptr<Chromosome> > chrs;
  boost::shared_ptr<Chromosome> chrT(new Chromosome(AASequence("TESTER"),1,0.5));
  chrs.push_back(chrT);
  boost::shared_ptr<Chromosome> chrA(new Chromosome(AASequence("ALLMER"),1,1.0));
  chrs.push_back(chrA);
  std::sort(chrs.begin(),chrs.end(),Chromosome::sortScoreDesc);
  TEST_STRING_EQUAL(chrs[0]->getSequence().toString(),"ALLMER");
  TEST_STRING_EQUAL(chrs[1]->getSequence().toString(),"TESTER");
  std::sort(chrs.begin(),chrs.end(),Chromosome::sortScoreAsc);
  TEST_STRING_EQUAL(chrs[0]->getSequence().toString(),"TESTER");
  TEST_STRING_EQUAL(chrs[1]->getSequence().toString(),"ALLMER");

}
END_SECTION
START_SECTION((Chromosome(AASequence seq, const int charge)))
{
  Chromosome chrA(AASequence("TESTER"),2);
  TEST_EQUAL(chrA.getSequence().toString(),"TESTER");
  TEST_EQUAL(chrA.getCharge(),2);
  TEST_REAL_SIMILAR(chrA.getScore(),-1);
  TEST_EQUAL(chrA.isScored(),false);
}
END_SECTION

START_SECTION((Chromosome(AASequence seq, const int charge, const double val)))
{
  Chromosome chrA(AASequence("TESTER"),2,0.7);
  TEST_EQUAL(chrA.getSequence().toString(),"TESTER");
  TEST_EQUAL(chrA.getCharge(),2);
  TEST_REAL_SIMILAR(chrA.getScore(),0.7);
  TEST_EQUAL(chrA.isScored(),true);
}
END_SECTION

START_SECTION((int getCharge() const ))
{
  TEST_EQUAL(ptr->getCharge(),1);
}
END_SECTION

START_SECTION((void setCharge(int charge)))
{
  ptr->setCharge(2);
  TEST_EQUAL(ptr->getCharge(),2);
}
END_SECTION

START_SECTION((void setScored(bool isScored)))
{
  TEST_EQUAL(ptr->isScored(),true);
  TEST_REAL_SIMILAR(ptr->getScore(),0.5);
  ptr->setScored(false);
  TEST_EQUAL(ptr->isScored(),false);
  TEST_REAL_SIMILAR(ptr->getScore(),-1);
}
END_SECTION

START_SECTION((bool isScored() const ))
{
	TEST_EQUAL(ptr->isScored(),false);
}
END_SECTION

START_SECTION(~Chromosome())
{
	delete ptr;
	delete ptt;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



