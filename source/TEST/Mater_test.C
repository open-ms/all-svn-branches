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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <boost/shared_ptr.hpp>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

struct TestMater :
	public Mater
{
		TestMater(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList) :
			Mater(precursorMass, precursorMassTolerance, aaList)
		{}

		std::vector<boost::shared_ptr<Chromosome> > mate(const boost::shared_ptr<Chromosome> lhs, const boost::shared_ptr<Chromosome> rhs) const
	    {
	    	std::vector<boost::shared_ptr<Chromosome> > ret;
	    	ret.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence(lhs->getSequence().toString() + rhs->getSequence().toString()),1,0.1)));
	    	ret.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence(rhs->getSequence().toString() + lhs->getSequence().toString()),1,0.1)));
	    	return(ret);
	    }
};
START_TEST(Mater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
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

Mater* ptr = 0;
Mater* null_ptr = 0;

START_SECTION((Mater(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList)))
{
	double precursorMass=1000;
	double precursorMassTolerance=20000;
	std::vector<const Residue*> aaList;
	ptr = new TestMater(precursorMass, precursorMassTolerance, aaList);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual std::vector<boost::shared_ptr<Chromosome> > mate(boost::shared_ptr< Chromosome > lhs, const boost::shared_ptr< Chromosome > rhs) const =0))
{
	std::vector<boost::shared_ptr<Chromosome> > res = ptr->mate(gp.getIndividual(0),gp.getIndividual(1));
	TEST_STRING_EQUAL(res[0]->getSequence().toString(),"TESTATESTR");
	TEST_STRING_EQUAL(res[1]->getSequence().toString(),"TESTRTESTA");
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  tm.seed(10);
  TEST_EQUAL(10,tm.getSeed());
}
END_SECTION

START_SECTION((unsigned int getSeed() const ))
{
  std::vector<const Residue*> aaList;
  TestMater tm(100.1, 0.3, aaList);
  tm.seed(10);
  TEST_EQUAL(10,tm.getSeed());
}
END_SECTION

START_SECTION((std::vector<boost::shared_ptr<Chromosome> > tournament(GenPool &genPool) const ))
{
	ptr->seed(10000);
	std::vector<boost::shared_ptr<Chromosome> > res = ptr->tournament(gp);
	TEST_EQUAL(res.size(),40);
}
END_SECTION

START_SECTION((void tournamentAndAddToPool(GenPool &genPool) const ))
{
	ptr->tournamentAndAddToPool(gp);
	TEST_EQUAL(gp.getPopulationSize(),56);
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  TEST_NOT_EQUAL(ptr->getSeed(),50);
  ptr->seed(50);
  TEST_EQUAL(ptr->getSeed(),50);
}
END_SECTION

START_SECTION((unsigned int getSeed() const ))
{
	TEST_EQUAL(ptr->getSeed(),50);
}
END_SECTION

START_SECTION(~Mater())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



