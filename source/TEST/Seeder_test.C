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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

///////////////////////////

using namespace OpenMS;
using namespace std;
struct TestSeeder :
	public Seeder
{
		TestSeeder(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList) :
			Seeder(precursorMass, precursorMassTolerance, aaList)
		{}

	    ~TestSeeder()
	    {}

	    Chromosome createIndividual() const
	    {
	    	return(Chromosome(AASequence("ALLMER"), 1.0));
	    }
};
START_TEST(Seeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Seeder* ptr = 0;
Seeder* null_ptr = 0;
START_SECTION(Seeder())
{
	double precursorMass=10;
	double precursorMassTolerance=2;
	std::vector<const Residue*> aaList;

	ptr = new TestSeeder(precursorMass,precursorMassTolerance,aaList);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Seeder())
{
	delete ptr;
}
END_SECTION

START_SECTION((Seeder(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~Seeder()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual Chromosome createIndividual()=0))
{
  // TODO
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  // TODO
}
END_SECTION

START_SECTION((const AASequence getRandomSequence(const int len, const double weight, const double tolerance)))
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

	int test[]{0,3000,6000,10000,14000,15000,16000};
	String res[]{"DPPPFGC","VVTRGSI","GEPQVFG","VGASNFH","VQQAASE","DDNIALA","GYAGAYM"};
	for(int i=0; i<7; i++) {
		TestSeeder ts(seq.getMonoWeight(Residue::Full),1.5,aaList);
		ts.seed(test[i]);
		AASequence rand = ts.getRandomSequence(7,seq.getMonoWeight(Residue::Full),1.5);
		TEST_EQUAL((std::abs(rand.getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5),true);
		TEST_EQUAL(7,rand.size());
		TEST_NOT_EQUAL(seq.toString(),rand.toString());
		TEST_EQUAL(res[i],rand.toString());
	}
}
END_SECTION

START_SECTION((const String getRandomAA() const ))
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
	for(int i=0; i<5; i++) {
		TestSeeder ts(seq.getMonoWeight(Residue::Full),1.5,aaList);
		ts.seed(test[i]);
		TEST_EQUAL(ts.getRandomAA(),res[i]);
	}
}
END_SECTION

START_SECTION((bool adjustToFitMass(AASequence &sequence, const double weight, const double tolerance) const ))
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
	String res[]  {"DPPPFGC","GEPQVFG","VGASNFH","VQQAASE","DDNIALA"};
	String inSeq[]{"DPPAFGC","GEPIVFG","VGASNFK","VQQAASR","DDNIALR"};
	String exSeq[]{"DAIAQIC","GEPQVFG","VGASNRK","VQQSGRG","DDNIALA"};
	for(int i=0; i<5; i++) {
		TestSeeder ts(seq.getMonoWeight(Residue::Full),1.5,aaList);
		ts.seed(test[i]);
		AASequence rand(inSeq[i]);
		ts.adjustToFitMass(rand,seq.getMonoWeight(Residue::Full),1.5);
		TEST_EQUAL((std::abs(rand.getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5),true);
		TEST_EQUAL(7,rand.size());
		TEST_NOT_EQUAL(inSeq[i],rand.toString());
		TEST_EQUAL(exSeq[i],rand.toString());
	}
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



