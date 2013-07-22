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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomSequenceSeeder.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RandomSequenceSeeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RandomSequenceSeeder* ptr = 0;
RandomSequenceSeeder* null_ptr = 0;
START_SECTION(RandomSequenceSeeder())
{
	double precursorMass=10;
	double precursorMassTolerance=2;
	std::vector<const Residue*> aaList;
	ptr = new RandomSequenceSeeder(precursorMass, precursorMassTolerance, aaList);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~RandomSequenceSeeder())
{
	delete ptr;
}
END_SECTION


START_SECTION(Chromosome RandomSequenceSeeder::createIndividual())
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
	String res[]{"IQWSIS","VVTRGSI","AGRWSR","VGASNFH","VQQAASE","DGDIALAG","GYAGAYM"};
	for(int i=0; i<7; i++) {
		RandomSequenceSeeder rss(seq.getMonoWeight(Residue::Full),1.5,aaList);
		rss.seed(test[i]);
		Chromosome rand = rss.createIndividual();
		TEST_EQUAL(abs(rand.getSequence().getMonoWeight(Residue::Full)-seq.getMonoWeight(Residue::Full))<1.5,true);
		TEST_EQUAL(rand.getSequence().toString(),res[i]);
	}

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



