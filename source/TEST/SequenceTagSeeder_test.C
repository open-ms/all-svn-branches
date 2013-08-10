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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SequenceTagSeeder.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

///////////////////////////

using namespace OpenMS;
using namespace std;
START_TEST(SequenceTagSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

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

double precursorMass = 1500.5;
double precursorMassTolerance = 2;
double fragmentMassTolerance = 0.5;
MSSpectrum<>* msms = new MSSpectrum<>();
//pm: 732.4073
//b: 185.1285	 298.2125	 429.2530	 558.2956
//y: 661.3702	 548.2861	 435.2020	 304.1615	 175.1190
Peak1D p1; p1.setMZ(185.1285); p1.setIntensity(3);
Peak1D p2; p2.setMZ(298.2125); p2.setIntensity(3);
Peak1D p3; p3.setMZ(429.2530); p3.setIntensity(3);
Peak1D p4; p4.setMZ(558.2956); p4.setIntensity(3);
Peak1D p5; p5.setMZ(661.3702); p5.setIntensity(5);
Peak1D p6; p6.setMZ(548.2861); p6.setIntensity(5);
Peak1D p7; p7.setMZ(435.2020); p7.setIntensity(5);
Peak1D p8; p8.setMZ(304.1615); p8.setIntensity(5);
Peak1D p9; p9.setMZ(175.1190); p9.setIntensity(5);
Peak1D p10; p10.setMZ(100.5000); p10.setIntensity(20);
msms->push_back(Peak1D(p1));
msms->push_back(Peak1D(p2));
msms->push_back(Peak1D(p3));
msms->push_back(Peak1D(p4));
msms->push_back(Peak1D(p5));
msms->push_back(Peak1D(p6));
msms->push_back(Peak1D(p7));
msms->push_back(Peak1D(p8));
msms->push_back(Peak1D(p9));
msms->push_back(Peak1D(p10));
msms->sortByPosition();

SequenceTagSeeder* ptr = 0;
SequenceTagSeeder* null_ptr = 0;

START_SECTION((SequenceTagSeeder(const MSSpectrum<> * spec, const double precursorMass, const double precursorMassTolerance, const double fagmentMassTolerance, std::vector< const Residue * > aaList)))
{
	ptr = new SequenceTagSeeder(msms,precursorMass,precursorMassTolerance,fragmentMassTolerance,aaList);
	ptr->seed(1000);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual boost::shared_ptr<Chromosome> createIndividual() const = 0))
{
	boost::shared_ptr<Chromosome> chr(ptr->createIndividual());
	TEST_STRING_EQUAL(chr->getSequence().toString(),"ALLMER");
}
END_SECTION

START_SECTION((virtual std::vector<boost::shared_ptr<Chromosome> > createIndividuals(const Size num) const))
{
	vector<boost::shared_ptr<Chromosome> > chrs(ptr->createIndividuals(10));
	TEST_EQUAL(chrs.size(),10);
	for(int i=0; i<10; i++)
		TEST_STRING_EQUAL(chrs[1]->getSequence().toString(),"ALLMER");
}
END_SECTION

START_SECTION(~SequenceTagSeeder())
{
	delete ptr;
	delete msms;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



