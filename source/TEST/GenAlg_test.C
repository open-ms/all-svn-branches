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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenAlg.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SubstitutingMutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomSequenceSeeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/NormShrAbuScorer.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleDecreasingKiller.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleMater.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/HomologyKiller.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <utility>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GenAlg, "$Id$")

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

MSSpectrum<> msms;
Precursor prec;
double pm = 1500.5;
double pmt = 2.0;
double fmt = 0.5;
prec.setMZ(pm);
prec.setCharge(1);
prec.setIntensity(100.0);
vector<Precursor> precursors;
precursors.push_back(prec);
msms.setPrecursors(precursors);

MSSpectrum<> msAll;
Precursor prec1;
double pm1 = AASequence("ALLMER").getMonoWeight(Residue::Full);
prec1.setMZ(pm1);
prec1.setCharge(1);
prec1.setIntensity(100.0);
vector<Precursor> precursors1;
precursors1.push_back(prec1);
msAll.setPrecursors(precursors1);

GenAlg* ptr = 0;
GenAlg* null_ptr = 0;

START_SECTION((GenAlg(const MSSpectrum<> *spec, const std::vector< const Residue * > &aaList, const int maxPoolSize, const double precursorMassTolerance, const double fragmentMassTolerance)))
{
	ptr = new GenAlg(&msAll, aaList, 20, pmt, fmt);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION
	
START_SECTION((void startEvolution(const Size numGenerations)))
{
  boost::shared_ptr<SubstitutingMutater> sm;
  ptr->setMutater(sm);	// break the default mutater set in GenAlg
  ptr->startEvolution(1);
  TEST_EQUAL(ptr->getKnownIndividuals().size(),0); // make sure it doesn't run if not all functions are set.
}
END_SECTION

START_SECTION((void setMutater(boost::shared_ptr< const Mutater > mutater)))
{
	boost::shared_ptr<SubstitutingMutater> sm(new SubstitutingMutater(pm1,pmt,aaList));
	ptr->setMutater(sm);
}
END_SECTION

START_SECTION((void setMater(boost::shared_ptr< const Mater > mater)))
{
	boost::shared_ptr<SimpleMater> mat(new SimpleMater(pm1,pmt,aaList));
	ptr->setMater(mat);
}
END_SECTION

START_SECTION((void setKiller(boost::shared_ptr< const Killer > killer)))
{
	boost::shared_ptr<SimpleDecreasingKiller> sk(new SimpleDecreasingKiller());
	ptr->setKiller(sk);
}
END_SECTION

START_SECTION((void setSeeder(boost::shared_ptr< const Seeder > seeder)))
{
	boost::shared_ptr<RandomSequenceSeeder> rss(new RandomSequenceSeeder(&msAll,pm1,pmt,fmt,aaList));
	ptr->setSeeder(rss);
}
END_SECTION

START_SECTION((void setScorer(boost::shared_ptr< const Scorer > scorer)))
{
	boost::shared_ptr<NormShrAbuScorer> ns(new NormShrAbuScorer(fmt));
	ptr->setScorer(ns);
}
END_SECTION

START_SECTION((void setHomologyKiller(boost::shared_ptr< const HomologyKiller > killer)))
{
	boost::shared_ptr<HomologyKiller> hk(new HomologyKiller());
	ptr->setHomologyKiller(hk);
}
END_SECTION

START_SECTION((void seed(const Size seed)))
{
  ptr->seed(10000);
}
END_SECTION

START_SECTION((const std::vector<const Residue*>& getAAList() const ))
{
  std::vector<const Residue*> got = ptr->getAAList();
  TEST_EQUAL(got.size(),aaList.size());
  TEST_EQUAL(got[0]->getOneLetterCode(),ptr->getAAList()[0]->getOneLetterCode());
}
END_SECTION

START_SECTION((const std::map<String, boost::shared_ptr<Chromosome> > getKnownIndividuals() const ))
{
	std::map<String, boost::shared_ptr<Chromosome> > ki = ptr->getKnownIndividuals();
	TEST_EQUAL(ki.size(),0);
}
END_SECTION

START_SECTION((void setKnownIndividuals(const std::map< String, boost::shared_ptr< Chromosome > > &knownIndividuals)))
{
	std::map<String, boost::shared_ptr<Chromosome> > ni;
	ni.insert(std::pair<String,boost::shared_ptr<Chromosome> >("ALLMER",boost::shared_ptr<Chromosome>(new Chromosome("ALLMER",1,1.0))));
	ni.insert(std::pair<String,boost::shared_ptr<Chromosome> >("ELLMER",boost::shared_ptr<Chromosome>(new Chromosome("ELLMER",1,0.7))));
	ptr->setKnownIndividuals(ni);
	TEST_EQUAL(ptr->getKnownIndividuals().size(),2);
}
END_SECTION

START_SECTION((double getPrecursorIntensity() const ))
{
	TEST_REAL_SIMILAR(ptr->getPrecursorIntensity(),100.0);
}
END_SECTION

START_SECTION((double getPrecursorMH() const ))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMH(),pm1);
}
END_SECTION

START_SECTION((double getPrecursorMass() const ))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMass(),(pm1-1));
}
END_SECTION

START_SECTION((double getPrecursorMz() const ))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMZ(),pm1);
}
END_SECTION

START_SECTION((void setRejuvenate(const Size afterNGen)))
{
  ptr->setRejuvenate(100);
}
END_SECTION
	
START_SECTION((void startEvolution(const Size numGenerations)))
{
  ptr->startEvolution(1);
  TEST_EQUAL(ptr->getKnownIndividuals().size(),43);
  TEST_EQUAL(ptr->getGenPool()->getPopulationSize(),15);
}
END_SECTION

START_SECTION((void startEvolution(const Size numGenerations, const Size endIfStableForNumGenerations)))
{
  ptr->startEvolution(1,2);
  TEST_EQUAL(ptr->getKnownIndividuals().size(),96);
  TEST_EQUAL(ptr->getGenPool()->getPopulationSize(),15);
}
END_SECTION

START_SECTION(~GenAlg())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




