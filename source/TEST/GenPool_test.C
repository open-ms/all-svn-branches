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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomSequenceSeeder.h>
#include <time.h>
#include <boost/shared_ptr.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GenPool, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GenPool* ptr = 0;
GenPool* null_ptr = 0;
START_SECTION((GenPool(const int maxPoolSize, const double precursorMassTolerance)))
{
	ptr = new GenPool(10,0.1);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~GenPool())
{
	delete ptr;
}
END_SECTION

START_SECTION((~GenPool()))
{
  // TODO
}
END_SECTION

START_SECTION((void initGenPool(const int maxPoolSize)))
{
	double pm = 1500.5;
	double pmt = 2;
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

	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	TEST_STRING_EQUAL("LNSFCGHWKPNAAG",gp.getIndividual(0)->getSequence().toString());
	TEST_STRING_EQUAL("YIQGQTSSIGGGSMD",gp.getIndividual(1)->getSequence().toString());
	TEST_STRING_EQUAL("WSVQMNKFTDMI",gp.getIndividual(2)->getSequence().toString());
	TEST_STRING_EQUAL("QQPEHRGEIIPVV",gp.getIndividual(3)->getSequence().toString());
	TEST_STRING_EQUAL("CNGTQQFETMIND",gp.getIndividual(4)->getSequence().toString());
}
END_SECTION

START_SECTION((void mutate()))
{
  // TODO
}
END_SECTION

START_SECTION((void setMutater(const Mutater &mutater)))
{
  // TODO
}
END_SECTION

START_SECTION((void crossover()))
{
  // TODO
}
END_SECTION

START_SECTION((void setMater(const Mater &mater)))
{
  // TODO
}
END_SECTION

START_SECTION((void kill()))
{
  // TODO
}
END_SECTION

START_SECTION((void setKiller(const Killer &killer)))
{
  // TODO
}
END_SECTION

START_SECTION((void setSeeder(const Seeder &seeder)))
{
	double pm = 1500.5;
	double pmt = 2;
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

	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	//gp.setSeeder(RandomSequenceSeeder(pm,pmt,aaList)); //Apparently not allowed
}
END_SECTION

START_SECTION((Size getPopulationSize()))
{
  // TODO
}
END_SECTION

START_SECTION((void sort(const int sortMethod)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPool(std::vector< Chromosome * > newPool)))
{
  // TODO
}
END_SECTION

START_SECTION((bool addIndividual(const Chromosome & individual)))
{
  GenPool gp(10,0.5);
  gp.setMSMSSpectrum(AASequence("ALLMER").getMonoWeight(Residue::Full));

  gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("ALLMER"),1.0)));
  gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTER"),1.0)));
  TEST_EQUAL(gp.getPopulationSize(),1);
}
END_SECTION

START_SECTION((void replenish(const int targetSize)))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<Chromosome *>::iterator begin()))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<Chromosome *>::iterator end()))
{
  // TODO
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  // TODO
}
END_SECTION

START_SECTION((unsigned int getSeed() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<const Residue*>& getAaList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setAaList(const std::vector< const Residue * > &aaList)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Chromosome*>& getGenPool() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setGenPool(const std::vector< Chromosome * > &genPool)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::map<String, Chromosome*>& getKnownIndividuals() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setKnownIndividuals(const std::map< String, Chromosome * > &knownIndividuals)))
{
  // TODO
}
END_SECTION

START_SECTION((int getMaxPoolSize() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setMaxPoolSize(int maxPoolSize)))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorMassTolerance() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorMassTolerance(double precursorMassTolerance)))
{
  // TODO
}
END_SECTION

START_SECTION((unsigned int getRandomSeed() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRandomSeed(unsigned int randomSeed)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



