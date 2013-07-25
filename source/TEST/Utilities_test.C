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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Utilities, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Utilities* ptr = 0;
Utilities* null_ptr = 0;
START_SECTION(Utilities())
{
	ptr = new Utilities();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Utilities())
{
	delete ptr;
}
END_SECTION

START_SECTION((static const AASequence getRandomSequence(const int seed, const int len, const double weight, const double tolerance, const std::vector< const Residue * > aaList)))
{
    AASequence seq("ALLMER");
	std::vector<const Residue*> aaList;
	aaList.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("N"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("D"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("C"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("E"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("Q"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("G"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("H"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("I"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("L"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("K"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("M"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("F"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("P"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("S"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("T"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("W"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("Y"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("V"));

	int test[]{0,3000,10000,12000,18000,20000};
	String 	expSeq[]{"DGIISID","VVDVVVC","VVVCVDV","HHHGMGG","IGIDIDS",""};
	int		expLen[]{7,7,7,7,7,0};
	bool	expRes[]{true,true,true,true,true,false};
	for(int i=0; i<6; i++) {
		AASequence rand = Utilities::getRandomSequence(test[i],7,seq.getMonoWeight(Residue::Full),1.5,aaList);
		TEST_EQUAL(expRes[i],(std::abs(rand.getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5));
		TEST_EQUAL(expLen[i],rand.size());
		TEST_NOT_EQUAL(seq.toString(),rand.toString());
		TEST_EQUAL(expSeq[i],rand.toString());
	}
}
END_SECTION

START_SECTION((static const String getRandomAA(const int seed, const std::vector< const Residue * > aaList)))
{
    AASequence seq("ALLMER");
	std::vector<const Residue*> aaList;
	aaList.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("N"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("D"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("C"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("E"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("Q"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("G"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("H"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("I"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("L"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("K"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("M"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("F"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("P"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("S"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("T"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("W"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("Y"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("V"));

	int test[]{0,6000,10000,14000,15000};
	String res[]  {"D","A","V","C","D"};
	for(int i=0; i<5; i++)
	{
		TEST_EQUAL(Utilities::getRandomAA(test[i],aaList),res[i]);
	}
}
END_SECTION

START_SECTION((static bool adjustToFitMass(const int seed, AASequence &sequence, const double weight, const double tolerance, const std::vector< const Residue * > aaList)))
{
    AASequence seq("ALLMER");
	std::vector<const Residue*> aaList;
	aaList.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("N"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("D"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("C"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("E"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("Q"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("G"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("H"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("I"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("L"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("K"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("M"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("F"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("P"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("S"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("T"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("W"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("Y"));
	aaList.push_back(ResidueDB::getInstance()->getResidue("V"));

	int test[]{0,6000,10000,14000,15000};
	String 	res[]  {"DPPPFGC","GEPQVFG","VGASNFH","VQQAASE","DDNIALA"};
	String 	inSeq[]{"DPPAFGC","GEPIVFG","VGASNFK","VQQAASR","DDNIALR"};
	String 	exSeq[]{"DQPAFGC","GEPQVFG","VGASNRK","VQMAAGR","DDNIALA"};
	bool 	expRes[] {false,true,true,true,true,true};
	for(int i=0; i<5; i++) {
		AASequence rand(inSeq[i]);
		Utilities::adjustToFitMass(test[i],rand,seq.getMonoWeight(Residue::Full),1.5,aaList);
		TEST_EQUAL(expRes[i],(std::abs(rand.getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5));
		TEST_EQUAL(7,rand.size());
		TEST_NOT_EQUAL(inSeq[i],rand.toString());
		TEST_EQUAL(exSeq[i],rand.toString());
	}

	AASequence m("AAAAAAAAAA");
	AASequence c("AAGGGLLLLL");
	Utilities::adjustToFitMass(50,c,m.getMonoWeight(Residue::Full),1.5,aaList);
	TEST_EQUAL((std::abs(c.getMonoWeight(Residue::Full)-m.getMonoWeight(Residue::Full))<1.5),false);
}
END_SECTION

START_SECTION((static int editDistance(const AASequence &lhs, const AASequence &rhs)))
{
  AASequence seq1("ALLMER");
  AASequence seq2("ALLMER");
  AASequence seq3("TELLER");
  AASequence seq4("TELLLLER");
  TEST_EQUAL(0,Utilities::editDistance(seq1,seq1));
  TEST_EQUAL(0,Utilities::editDistance(seq1,seq2));
  TEST_EQUAL(3,Utilities::editDistance(seq1,seq3));
  TEST_EQUAL(2,Utilities::editDistance(seq4,seq3));
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



