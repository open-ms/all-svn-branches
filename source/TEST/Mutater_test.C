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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <boost/shared_ptr.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

struct TestMutater :
	public Mutater
{
		TestMutater(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList) :
			Mutater(precursorMass, precursorMassTolerance, aaList)
		{}

		void mutate(boost::shared_ptr<Chromosome> chromosome) const
	    {
			AASequence seq(chromosome->getSequence());
			seq.setResidue(0,ResidueDB::getInstance()->getResidue("A"));
			chromosome->setSequence(seq);
	    }
};

START_TEST(Mutater, "$Id$")

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

	GenPool gp(500,1000,1000,1);
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTA"),1,0.1)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTR"),1,0.15)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTN"),1,0.17)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTD"),1,0.2)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTC"),1,0.25)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("RANDO"),1,0.27)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("AAAAA"),1,0.3)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("ALLME"),1,0.35)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WWGGG"),1,0.37)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTO"),1,0.4)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTI"),1,0.45)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTL"),1,0.47)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTK"),1,0.5)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTM"),1,0.55)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTF"),1,0.57)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTP"),1,0.6)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTU"),1,0.65)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTS"),1,0.67)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTT"),1,0.7)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTW"),1,0.75)));

Mutater* ptr = 0;
Mutater* null_ptr = 0;

START_SECTION((Mutater(double precursorMass, double precursorMassTolerance, std::vector< Residue & > aaList)))
{
	std::vector< const Residue * > aaList;
	ptr = new TestMutater(1000, 1000, aaList);
	TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION((double getMutationRate() const ))
{
	TEST_REAL_SIMILAR(ptr->getMutationRate(),0.2);
}
END_SECTION

START_SECTION((void setMutationRate(double mutationRate=0.2)))
{
  ptr->setMutationRate(1.0);
  TEST_REAL_SIMILAR(ptr->getMutationRate(),1.0);
}
END_SECTION

START_SECTION((void mutatePool(GenPool & pool) const))
{
	ptr->mutatePool(gp);
	TEST_EQUAL(gp.getPopulationSize(),20);
	for(Size i=0; i < gp.getPopulationSize(); i++)
	{
		TEST_EQUAL(gp.getIndividual(i)->getSequence().toString()[0],'A');
	}
}
END_SECTION

START_SECTION((void mutateAndAddToPool(GenPool & pool) const))
{
	ptr->mutateAndAddToPool(gp);
	TEST_EQUAL(gp.getPopulationSize(),38);
}
END_SECTION

START_SECTION((virtual void mutate(boost::shared_ptr<Chromosome> chromosome) const = 0))
{
	boost::shared_ptr<Chromosome> chromosome(new Chromosome(AASequence("ELLMER"),1,1.0));
	ptr->mutate(chromosome);
	TEST_STRING_EQUAL(chromosome->getSequence().toString(),"ALLMER");
}
END_SECTION

START_SECTION((boost::shared_ptr<Chromosome> mutateCpy(const boost::shared_ptr<Chromosome> chromosome) const))
{
	boost::shared_ptr<Chromosome> chromosome(new Chromosome(AASequence("ELLMER"),1,1.0));
	boost::shared_ptr<Chromosome> res(ptr->mutateCpy(chromosome));
	TEST_STRING_EQUAL(res->getSequence().toString(),"ALLMER");
	TEST_STRING_EQUAL(chromosome->getSequence().toString(),"ELLMER");
}
END_SECTION
	
START_SECTION((const std::vector<const Residue*>& getAAList() const))
{
  TEST_EQUAL(ptr->getAAList().size(),0);
}
END_SECTION
	
START_SECTION((double getPrecursorMass() const))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMass(),1000);
}
END_SECTION

START_SECTION((double getPrecursorMassTolerance() const))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMassTolerance(),1000);
}
END_SECTION
	
START_SECTION((void seed(const Size seed)))
{
  Size test[] = {0,6000,10000,14000,15000};
  String res[] = {"L","Q","L","V","D"};
  for(int i=0; i<5; i++)
  {
	  ptr->seed(test[i]);
	  const Utilities * utils = ptr->getUtils();
	  TEST_EQUAL(utils->getRandomAA(aaList),res[i]);
  }
}
END_SECTION
	
START_SECTION((const Utilities * getUtils() const))
{
	const Utilities * utils = ptr->getUtils();
	TEST_EQUAL(utils->editDistance("ALLMER","ELLMER"),1);
}
END_SECTION

START_SECTION(~Mutater())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



