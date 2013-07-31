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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <boost/shared_ptr.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

struct TestMater :
	public Mater
{
		TestMater(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList) :
			Mater(precursorMass, precursorMassTolerance, aaList)
		{}

	    std::vector<boost::shared_ptr<Chromosome> > mate(boost::shared_ptr<Chromosome> lhs, boost::shared_ptr<Chromosome> rhs)
	    {
	    	std::vector<boost::shared_ptr<Chromosome> > ret;
	    	return(ret);
	    }
};
START_TEST(Mater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Mater* ptr = 0;
Mater* null_ptr = 0;
START_SECTION(Mater())
{
	double precursorMass=10;
	double precursorMassTolerance=2;
	std::vector<const Residue*> aaList;

	ptr = new TestMater(precursorMass, precursorMassTolerance, aaList);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Mater())
{
	delete ptr;
}
END_SECTION

START_SECTION((Mater(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList)))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<Chromosome> mate(const Chromosome &lhs, const Chromosome &rhs)=0))
{
  // TODO
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  tm.seed(10);
  TEST_EQUAL(10,tm.getSeed());
}
END_SECTION

START_SECTION((unsigned int getSeed() const ))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  tm.seed(10);
  TEST_EQUAL(10,tm.getSeed());
}
END_SECTION

START_SECTION((const std::vector<const Residue*>& getAAList() const ))
{
	std::vector<const Residue*> aaList1a;
	aaList1a.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList1a.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList1a.push_back(ResidueDB::getInstance()->getResidue("N"));
	std::vector<const Residue*> aaList1b;
	aaList1b.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList1b.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList1b.push_back(ResidueDB::getInstance()->getResidue("N"));
	std::vector<const Residue*> aaList2;
	aaList2.push_back(ResidueDB::getInstance()->getResidue("D"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("C"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("E"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("Q"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("G"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("H"));
	TestMater tm(100.1, 0.3, aaList1a);
	TEST_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList1b[1]->getModifiedOneLetterCode());
	TEST_NOT_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList2[1]->getModifiedOneLetterCode());
	tm.setAAList(aaList2);
	TEST_NOT_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList1b[1]->getModifiedOneLetterCode());
	TEST_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList2[1]->getModifiedOneLetterCode());
}
END_SECTION

START_SECTION((void setAAList(const std::vector< const Residue * > &aaList)))
{
	std::vector<const Residue*> aaList1a;
	aaList1a.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList1a.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList1a.push_back(ResidueDB::getInstance()->getResidue("N"));
	std::vector<const Residue*> aaList1b;
	aaList1b.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList1b.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList1b.push_back(ResidueDB::getInstance()->getResidue("N"));
	std::vector<const Residue*> aaList2;
	aaList2.push_back(ResidueDB::getInstance()->getResidue("D"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("C"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("E"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("Q"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("G"));
	aaList2.push_back(ResidueDB::getInstance()->getResidue("H"));
	TestMater tm(100.1, 0.3, aaList1a);
	TEST_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList1b[1]->getModifiedOneLetterCode());
	TEST_NOT_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList2[1]->getModifiedOneLetterCode());
	tm.setAAList(aaList2);
	TEST_NOT_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList1b[1]->getModifiedOneLetterCode());
	TEST_EQUAL(tm.getAAList()[1]->getModifiedOneLetterCode(),aaList2[1]->getModifiedOneLetterCode());
}
END_SECTION

START_SECTION((double getPrecursorMass() const ))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  TEST_EQUAL(100.1,tm.getPrecursorMass());
  tm.setPrecursorMass(200.2);
  TEST_EQUAL(200.2,tm.getPrecursorMass());
}
END_SECTION

START_SECTION((void setPrecursorMass(double precursorMass)))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  tm.setPrecursorMass(200.2);
  TEST_EQUAL(200.2,tm.getPrecursorMass());
}
END_SECTION

START_SECTION((double getPrecursorMassTolerance() const ))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  TEST_EQUAL(0.3,tm.getPrecursorMassTolerance());
  tm.setPrecursorMassTolerance(0.4);
  TEST_EQUAL(0.4,tm.getPrecursorMassTolerance());
}
END_SECTION

START_SECTION((void setPrecursorMassTolerance(double precursorMassTolerance)))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  tm.setPrecursorMassTolerance(0.4);
  TEST_EQUAL(0.4,tm.getPrecursorMassTolerance());
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



