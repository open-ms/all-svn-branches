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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SubstitutingMutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleMater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleDecreasingKiller.h>
#include <time.h>
#include <boost/shared_ptr.hpp>
#include <map>
#include <utility>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/NormShrAbuScorer.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GenPool, "$Id$")

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

START_SECTION((void initGenPool(const int maxPoolSize)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
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
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	SubstitutingMutater sm(pm,pmt,aaList);
	gp.setMutater(&sm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	gp.mutate();
	TEST_STRING_EQUAL("LNSFCGYWKANAAG",gp.getIndividual(5)->getSequence().toString()); //This one is actually mutated from LNSFCGHWKPNAAG
}
END_SECTION

START_SECTION((void setMutater(Mutater * mutater)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	SubstitutingMutater sm(pm,pmt,aaList);
	gp.setMutater(&sm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	gp.mutate();
	TEST_STRING_EQUAL("LNSFCGYWKANAAG",gp.getIndividual(5)->getSequence().toString()); //This one is actually mutated from LNSFCGHWKPNAAG
}
END_SECTION

START_SECTION((void crossover()))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	SimpleMater sm(pm,pmt,aaList);
	sm.seed(rs);
	gp.setMater(&sm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	gp.crossover();
	TEST_EQUAL(gp.getPopulationSize(),13);
	TEST_STRING_EQUAL("LNSGCGMWKGGGSMD",gp.getIndividual(5)->getSequence().toString()); //This one of the children
}
END_SECTION

START_SECTION((void setMater(Mater * mater)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	SimpleMater sm(pm,pmt,aaList);
	sm.seed(rs);
	gp.setMater(&sm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	gp.crossover();
	TEST_EQUAL(gp.getPopulationSize(),13);
	//for(int i=0; i < gp.getPopulationSize(); i++) {
	//	cerr << gp.getIndividual(i)->getSequence().toString() << endl;
	//}
	TEST_STRING_EQUAL("LNSGCGMWKGGGSMD",gp.getIndividual(5)->getSequence().toString()); //This one of the children
}
END_SECTION

START_SECTION((void kill()))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(20,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	SimpleMater sm(pm,pmt,aaList);
	sm.seed(rs);
	gp.setMater(&sm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	gp.crossover();
	TEST_EQUAL(gp.getPopulationSize(),13);
	gp.setPreviousPoolSize(13);
	SimpleDecreasingKiller sdk;
	gp.setKiller(&sdk);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),8);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),3);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),1);
}
END_SECTION

START_SECTION((void setKiller(Killer * killer)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(20,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	SimpleMater sm(pm,pmt,aaList);
	sm.seed(rs);
	gp.setMater(&sm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	gp.crossover();
	TEST_EQUAL(gp.getPopulationSize(),13);
	gp.setPreviousPoolSize(13);
	SimpleDecreasingKiller sdk;
	gp.setKiller(&sdk);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),8);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),3);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),1);
}
END_SECTION

START_SECTION((void setSeeder(Seeder * seeder)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
	//gp.setSeeder(RandomSequenceSeeder(pm,pmt,aaList)); //Apparently not allowed
}
END_SECTION

START_SECTION((Size getPopulationSize()))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(15,pmt);
	TEST_EQUAL(gp.getPopulationSize(),0);
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
}
END_SECTION

START_SECTION((void sort(const int sortMethod)))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	gp.sort();
	TEST_STRING_EQUAL(gp.getIndividual(0)->getSequence().toString(),"QQPEHRGEIIPVV");
	gp.sort(Chromosome::sortScoreAscending);
	TEST_STRING_EQUAL(gp.getIndividual(0)->getSequence().toString(),"CNGTQQFETMIND");
	gp.sort(Chromosome::sortScoreDescending);
    TEST_STRING_EQUAL(gp.getIndividual(0)->getSequence().toString(),"QQPEHRGEIIPVV");
}
END_SECTION

START_SECTION((void setGenPool(std::vector<boost::shared_ptr<Chromosome> > newPool)))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	std::vector<boost::shared_ptr<Chromosome> > ngp;
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	gp.setGenPool(ngp);
	gp.sort();
	TEST_STRING_EQUAL(gp.getIndividual(0)->getSequence().toString(),"QQPEHRGEIIPVV");
	TEST_EQUAL(gp.getPopulationSize(),5);
}
END_SECTION

START_SECTION((bool addIndividual(const Chromosome & individual)))
{
  GenPool gp(10,0.5);
  gp.setMSMSSpectrum(AASequence("ALLMER").getMonoWeight(Residue::Full));
  gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("ALLMER"),1.0)));
  gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTER"),1.0))); //Won't be added since not within precursor mass tolerance.
  TEST_EQUAL(gp.getPopulationSize(),1);
  TEST_STRING_EQUAL(gp.getIndividual(0)->getSequence().toString(),"ALLMER");
}
END_SECTION

START_SECTION((void replenish(const int targetSize)))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	gp.sort();
	SimpleDecreasingKiller sdk;
	gp.setKiller(&sdk);
	gp.kill();
	TEST_EQUAL(gp.getPopulationSize(),1);
	gp.replenish(5);
	TEST_EQUAL(gp.getPopulationSize(),5);
}
END_SECTION

START_SECTION((std::vector<boost::shared_ptr<Chromosome> >::iterator begin()))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	boost::shared_ptr<Chromosome> input[5];
	input[0] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3));
	input[1] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4));
	input[2] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2));
	input[3] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5));
	input[4] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1));
	gp.addIndividual(input[0]);
	gp.addIndividual(input[1]);
	gp.addIndividual(input[2]);
	gp.addIndividual(input[3]);
	gp.addIndividual(input[4]);
	int i = 0;
	for(std::vector<boost::shared_ptr<Chromosome> >::iterator iter=gp.begin(); iter!=gp.end(); iter++) {
		TEST_STRING_EQUAL((*iter)->getSequence().toString(),input[i++]->getSequence().toString());
	}
}
END_SECTION

START_SECTION((std::vector<boost::shared_ptr<Chromosome> >::iterator end()))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	boost::shared_ptr<Chromosome> input[5];
	input[0] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3));
	input[1] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4));
	input[2] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2));
	input[3] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5));
	input[4] = boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1));
	gp.addIndividual(input[0]);
	gp.addIndividual(input[1]);
	gp.addIndividual(input[2]);
	gp.addIndividual(input[3]);
	gp.addIndividual(input[4]);
	int i = 0;
	for(std::vector<boost::shared_ptr<Chromosome> >::iterator iter=gp.begin(); iter!=gp.end(); iter++) {
		TEST_STRING_EQUAL((*iter)->getSequence().toString(),input[i++]->getSequence().toString());
	}
	gp.sort();
	std::vector<boost::shared_ptr<Chromosome> >::iterator te = gp.end();
	--te;
	TEST_STRING_EQUAL("CNGTQQFETMIND",te->get()->getSequence().toString());
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  GenPool gp(15,1500.4);
  Size cs = gp.getSeed();
  gp.seed(100);
  TEST_NOT_EQUAL(gp.getSeed(),cs);
}
END_SECTION

START_SECTION((unsigned int getSeed() const ))
{
  GenPool gp(15,1500.4);
  Size cs = gp.getSeed();
  gp.seed(100);
  TEST_NOT_EQUAL(gp.getSeed(),cs);
  cs = gp.getSeed();
  TEST_EQUAL(gp.getSeed(),cs);
}
END_SECTION

START_SECTION((const std::vector<const Residue*>& getAaList() const ))
{
	GenPool gp(15,1500.4);
	TEST_EQUAL(gp.getAaList().size(),0);
	gp.setAaList(aaList);
	TEST_EQUAL(gp.getAaList().size(),20);
	TEST_STRING_EQUAL(gp.getAaList()[5]->getOneLetterCode(),"E");
}
END_SECTION

START_SECTION((void setAaList(const std::vector< const Residue * > &aaList)))
{
	GenPool gp(15,1500.4);
	TEST_EQUAL(gp.getAaList().size(),0);
	gp.setAaList(aaList);
	TEST_EQUAL(gp.getAaList().size(),20);
	TEST_STRING_EQUAL(gp.getAaList()[5]->getOneLetterCode(),"E");
	std::vector<const Residue*> aaList1;
	aaList1.push_back(ResidueDB::getInstance()->getResidue("A"));
	aaList1.push_back(ResidueDB::getInstance()->getResidue("R"));
	aaList1.push_back(ResidueDB::getInstance()->getResidue("N"));
	aaList1.push_back(ResidueDB::getInstance()->getResidue("D"));
	gp.setAaList(aaList1);
	TEST_EQUAL(gp.getAaList().size(),4);
	TEST_STRING_EQUAL(gp.getAaList()[2]->getOneLetterCode(),"N");
}
END_SECTION

START_SECTION((const std::vector<boost::shared_ptr<Chromosome> >& getGenPool() const ))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	std::vector<boost::shared_ptr<Chromosome> > ngp;
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	gp.setGenPool(ngp);
	std::vector<boost::shared_ptr<Chromosome> > pool = gp.getGenPool();
	TEST_STRING_EQUAL(pool[0]->getSequence().toString(),"LNSFCGHWKPNAAG");
	TEST_STRING_EQUAL(pool[1]->getSequence().toString(),"YIQGQTSSIGGGSMD");
	TEST_STRING_EQUAL(pool[2]->getSequence().toString(),"WSVQMNKFTDMI");
	TEST_STRING_EQUAL(pool[3]->getSequence().toString(),"QQPEHRGEIIPVV");
	TEST_STRING_EQUAL(pool[4]->getSequence().toString(),"CNGTQQFETMIND");
}
END_SECTION

START_SECTION((void setGenPool(const std::vector<boost::shared_ptr<Chromosome> > &genPool)))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	std::vector<boost::shared_ptr<Chromosome> > ngp;
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	gp.setGenPool(ngp);
	std::vector<boost::shared_ptr<Chromosome> > pool = gp.getGenPool();
	TEST_STRING_EQUAL(pool[0]->getSequence().toString(),"LNSFCGHWKPNAAG");
	TEST_STRING_EQUAL(pool[1]->getSequence().toString(),"YIQGQTSSIGGGSMD");
	TEST_STRING_EQUAL(pool[2]->getSequence().toString(),"WSVQMNKFTDMI");
	TEST_STRING_EQUAL(pool[3]->getSequence().toString(),"QQPEHRGEIIPVV");
	TEST_STRING_EQUAL(pool[4]->getSequence().toString(),"CNGTQQFETMIND");
}
END_SECTION

START_SECTION((const std::map<String, boost::shared_ptr<Chromosome> >& getKnownIndividuals() const ))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	std::map<String, boost::shared_ptr<Chromosome> > ki = gp.getKnownIndividuals();
	TEST_STRING_EQUAL(ki.find("LNSFCGHWKPNAAG")->second->getSequence().toString(),"LNSFCGHWKPNAAG");
	TEST_STRING_EQUAL(ki.find("QQPEHRGEIIPVV")->second->getSequence().toString(),"QQPEHRGEIIPVV");
	TEST_EQUAL(5,ki.size());
}
END_SECTION

START_SECTION((void setKnownIndividuals(const std::map< String, boost::shared_ptr<Chromosome> > &knownIndividuals)))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	gp.setMSMSSpectrum(pm);
	std::map<String, boost::shared_ptr<Chromosome> > ki;
	boost::shared_ptr<Chromosome> ind1(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3));
	ki.insert(std::pair<String,boost::shared_ptr<Chromosome> >(ind1->getSequence().toString(),ind1));
	boost::shared_ptr<Chromosome> ind2(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4));
	ki.insert(std::pair<String,boost::shared_ptr<Chromosome> >(ind2->getSequence().toString(),ind2));
	boost::shared_ptr<Chromosome> ind3(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2));
	ki.insert(std::pair<String,boost::shared_ptr<Chromosome> >(ind3->getSequence().toString(),ind3));
	boost::shared_ptr<Chromosome> ind4(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5));
	ki.insert(std::pair<String,boost::shared_ptr<Chromosome> >(ind4->getSequence().toString(),ind4));
	boost::shared_ptr<Chromosome> ind5(new Chromosome(AASequence("CNGTQQFETMIND"),0.1));
	ki.insert(std::pair<String,boost::shared_ptr<Chromosome> >(ind5->getSequence().toString(),ind5));
	gp.setKnownIndividuals(ki);
	TEST_EQUAL(gp.getPopulationSize(),0);
	gp.replenish(5);
	gp.sort();
	std::vector<boost::shared_ptr<Chromosome> > pool = gp.getGenPool();
	TEST_STRING_EQUAL(pool[0]->getSequence().toString(),"LNSFCGHWKPNAAG"); //duplicates may occur if the known individuals are so few and
	TEST_STRING_EQUAL(pool[1]->getSequence().toString(),"LNSFCGHWKPNAAG"); //a random selection is done from them.
	TEST_STRING_EQUAL(pool[2]->getSequence().toString(),"WSVQMNKFTDMI");
	TEST_STRING_EQUAL(pool[3]->getSequence().toString(),"WSVQMNKFTDMI");
	TEST_STRING_EQUAL(pool[4]->getSequence().toString(),"CNGTQQFETMIND");
}
END_SECTION

START_SECTION((int getMaxPoolSize() const ))
{
	GenPool gp(15,1500.4);
	TEST_EQUAL(gp.getMaxPoolSize(),15);
	gp.setMaxPoolSize(100);
	TEST_NOT_EQUAL(gp.getMaxPoolSize(),15);
	TEST_EQUAL(gp.getMaxPoolSize(),100);
}
END_SECTION

START_SECTION((void setMaxPoolSize(int maxPoolSize)))
{
	GenPool gp(15,1500.4);
	TEST_EQUAL(gp.getMaxPoolSize(),15);
	gp.setMaxPoolSize(100);
	TEST_NOT_EQUAL(gp.getMaxPoolSize(),15);
	TEST_EQUAL(gp.getMaxPoolSize(),100);
}
END_SECTION

START_SECTION((double getPrecursorMassTolerance() const ))
{
	  GenPool gp(15,1500.4);
	  double pmt = gp.getPrecursorMassTolerance();
	  gp.setPrecursorMassTolerance(100);
	  TEST_EQUAL(abs(gp.getPrecursorMassTolerance()-pmt)<0.001,false);
	  pmt = gp.getPrecursorMassTolerance();
	  TEST_REAL_SIMILAR(gp.getPrecursorMassTolerance(),pmt);
}
END_SECTION

START_SECTION((void setPrecursorMassTolerance(double precursorMassTolerance)))
{
	  GenPool gp(15,1500.4);
	  double pmt = gp.getPrecursorMassTolerance();
	  gp.setPrecursorMassTolerance(100);
	  TEST_EQUAL(abs(gp.getPrecursorMassTolerance()-pmt)<0.001,false);
	  pmt = gp.getPrecursorMassTolerance();
	  TEST_REAL_SIMILAR(gp.getPrecursorMassTolerance(),pmt);
}
END_SECTION

START_SECTION((void setPreviousPoolSize(unsigned int previousPoolSize)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	TEST_EQUAL(0,gp.getPreviousPoolSize());
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
    TEST_EQUAL(5,gp.getPreviousPoolSize());
    gp.setPreviousPoolSize(20);
    TEST_EQUAL(20,gp.getPreviousPoolSize());
}
END_SECTION

START_SECTION((unsigned int getPreviousPoolSize()))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	TEST_EQUAL(0,gp.getPreviousPoolSize());
	gp.setRandomSeed(rs);
	gp.setMSMSSpectrum(pm);
	RandomSequenceSeeder rss(pm,pmt,aaList);
	rss.seed(gp.getRandomSeed());
	gp.setSeeder(&rss);
	gp.initGenPool(5);
    TEST_EQUAL(5,gp.getPreviousPoolSize());
    gp.setPreviousPoolSize(20);
    TEST_EQUAL(20,gp.getPreviousPoolSize());
}
END_SECTION

START_SECTION((void GenPool::setScorer(const Scorer* scoer)))
{
	double pm = 1500.5;
	double pmt = 2;
	int rs = 84864;
	GenPool gp(5,pmt);
	NormShrAbuScorer scorer(0.5);
	gp.setScorer(&scorer);
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



