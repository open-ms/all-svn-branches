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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomMater.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <stdlib.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RandomMater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RandomMater* ptr = 0;
RandomMater* null_ptr = 0;
START_SECTION((RandomMater(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList)))
{
	std::vector<const Residue*> aaList;
	ptr = new RandomMater(100.0,0.3,aaList);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~RandomMater())
{
	delete ptr;
}
END_SECTION

START_SECTION((std::vector<Chromosome> mate(const Chromosome &lhs, const Chromosome &rhs)))
{
	AASequence m("AAAAAAAAAA");
	AASequence f("GGGGGLLLLL");
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
	RandomMater rm(m.getMonoWeight(Residue::Full),1.5,aaList);
	rm.seed(150);
	vector<boost::shared_ptr<Chromosome> > res = rm.mate(boost::shared_ptr<Chromosome>(new Chromosome(m,0)), boost::shared_ptr<Chromosome>(new Chromosome(f,0)));
	TEST_STRING_EQUAL("AAAAAAAAAA",m.toString());
	TEST_STRING_EQUAL("GGGGGLLLLL",f.toString());
	TEST_EQUAL(2,res.size());
	TEST_STRING_EQUAL("AAAAGGAAGL",res[0]->getSequence().toString());
	TEST_STRING_EQUAL("GGGGGALLGA",res[1]->getSequence().toString());
	TEST_EQUAL(true,(std::abs(m.getMonoWeight(Residue::Full)-res[0]->getSequence().getMonoWeight(Residue::Full)) < 1.5));
	TEST_EQUAL(true,(std::abs(m.getMonoWeight(Residue::Full)-res[1]->getSequence().getMonoWeight(Residue::Full)) < 1.5));

	AASequence bt("AAAAAA");
	AASequence bc("WWWWW");
	res = rm.mate(boost::shared_ptr<Chromosome>(new Chromosome(bt,0)), boost::shared_ptr<Chromosome>(new Chromosome(bc,0)));
	TEST_EQUAL(0,res.size()); //a different seed may lead to 2 valid children of m and f.
}
END_SECTION

START_SECTION((const std::vector<double> getWeights() const ))
{
	std::vector<const Residue*> aaList;
    AASequence aas("WLQSEVIHAR");
    RandomMater rm(aas.getMonoWeight(),0.3,aaList);
    vector<double> res = rm.getWeights();
    TEST_REAL_SIMILAR(res[0],0.4);
    TEST_REAL_SIMILAR(res[1],0.3);
    TEST_REAL_SIMILAR(res[2],0.3);
}
END_SECTION

START_SECTION((void setWeights(const std::vector< double > &weights)))
{
	std::vector<const Residue*> aaList;
    AASequence aas("WLQSEVIHAR");
    RandomMater rm(aas.getMonoWeight(),0.3,aaList);
    vector<double> res = rm.getWeights();
    res[0] = 0.1;
    res[1] = 0.5;
    res[2] = 0.4;
    RandomMater rm1(aas.getMonoWeight(),0.3,aaList);
    rm1.setWeights(res);
    res = rm1.getWeights();
    TEST_REAL_SIMILAR(res[0],0.1);
    TEST_REAL_SIMILAR(res[1],0.5);
    TEST_REAL_SIMILAR(res[2],0.4);
	res = rm.getWeights();
    TEST_REAL_SIMILAR(res[0],0.4);
    TEST_REAL_SIMILAR(res[1],0.3);
    TEST_REAL_SIMILAR(res[2],0.3);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



