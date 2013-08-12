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
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Utilities, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

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

Utilities* ptr = 0;
Utilities* null_ptr = 0;

START_SECTION(Utilities())
{
	ptr = new Utilities();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((static const AASequence getRandomSequence(const int seed, const int len, const double weight, const double tolerance, const std::vector< const Residue * > aaList)))
{

	int test[] = {10000,12000,14000,66666};
	String 	expSeq[] = {"","EGHYIGG","CIDDPLG","HSGVSMD"};
	int		expLen[] = {0,7,7,7};
	bool	expRes[] = {false,true,true,true};
	for(int i=0; i<4; i++) {
		ptr->seed(test[i]);
		AASequence rand = ptr->getRandomSequence(7,seq.getMonoWeight(Residue::Full),1.5,aaList);
		TEST_STRING_EQUAL(rand.toString(), expSeq[i]);
		TEST_NOT_EQUAL(rand.toString(), seq.toString());
		TEST_EQUAL(expLen[i],rand.size());
		TEST_EQUAL(expRes[i],(std::abs(rand.getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5));
	}
}
END_SECTION

START_SECTION((static const String getRandomAA(const int seed, const std::vector< const Residue * > aaList)))
{

	int test[] = {0,6000,10000,14000,15000};
	String res[] = {"L","Q","L","V","D"};
	for(int i=0; i<5; i++)
	{
		ptr->seed(test[i]);
		TEST_EQUAL(ptr->getRandomAA(aaList),res[i]);
	}
}
END_SECTION

START_SECTION((static bool adjustToFitMass(const int seed, AASequence &sequence, const double weight, const double tolerance, const std::vector< const Residue * > aaList)))
{
    
	int test[] = {0,6000,10000,14000,15000};
	String 	inSeq[] = {"DPPAFGC","GEPIVFG","VGASNFK","VQQAASR","DDNIALR"};
	String 	exSeq[] = {"DPPPFGC","GQIIVFG","VGAPNFK","VQQAASE","SGNDALR"};
	bool 	expRes[] = {true,true,true,true,true,true};
	for(int i=0; i<5; i++) {
		ptr->seed(test[i]);
		AASequence rand(inSeq[i]);
		ptr->adjustToFitMass(rand,seq.getMonoWeight(Residue::Full),1.5,aaList);
		TEST_EQUAL(expRes[i],(std::abs(rand.getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5));
		TEST_EQUAL(7,rand.size());
		TEST_NOT_EQUAL(inSeq[i],rand.toString());
		TEST_EQUAL(exSeq[i],rand.toString());
	}

	AASequence m("AAAAAAAAAA");
	AASequence c("AAGGGLLLLL");
	ptr->seed(5080);
	bool success = ptr->adjustToFitMass(c,m.getMonoWeight(Residue::Full),0.05,aaList);
	TEST_EQUAL(success,false);
}
END_SECTION

START_SECTION((static int editDistance(const AASequence &lhs, const AASequence &rhs)))
{
  AASequence seq1("ALLMER");
  AASequence seq2("ALLMER");
  AASequence seq3("TELLER");
  AASequence seq4("TELLLLER");
  TEST_EQUAL(0,ptr->editDistance(seq1,seq1));
  TEST_EQUAL(0,ptr->editDistance(seq1,seq2));
  TEST_EQUAL(3,ptr->editDistance(seq1,seq3));
  TEST_EQUAL(2,ptr->editDistance(seq4,seq3));
}
END_SECTION

START_SECTION((static double getSummedIntensity(const MSSpectrum<> *ms)))
{
	TheoreticalSpectrumGenerator tsg;
	RichPeakSpectrum theoMS;
	tsg.getSpectrum(theoMS,AASequence("ALLMER"),1);
	TEST_EQUAL(theoMS.size(),9); //just b,y peaks, but all of them.
	MSSpectrum<> ms;
	for(std::vector<RichPeak1D>::const_iterator iter = theoMS.begin(); iter != theoMS.end(); iter++)
	{
		ms.push_back(Peak1D(*iter));
	}
	Param param;
	double tic = ptr->getSummedIntensity(&ms);
	TEST_REAL_SIMILAR(tic,9);
}
END_SECTION

START_SECTION(~Utilities())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



