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
#include <boost/shared_ptr.hpp>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleDecreasingKiller.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GenPool, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GenPool* ptr = 0;
GenPool* null_ptr = 0;

START_SECTION((GenPool(const int maxPoolSize, const double precursorMassTolerance, const double precursorMass, const int precursorCharge)))
{
  ptr = new GenPool(10,2.0,1500.5,1);
  TEST_NOT_EQUAL(ptr,null_ptr);
}
END_SECTION

START_SECTION((Size getPopulationSize()))
{
	TEST_EQUAL(10,ptr->getMaxPoolSize());
}
END_SECTION

START_SECTION((bool addIndividual(boost::shared_ptr< Chromosome > individual)))
{
	ptr->addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("ALLMER"),1,0.3)));
	TEST_EQUAL(ptr->getPopulationSize(),0); //Above shouldn't be added since the precursor mass doesn'T match.
	ptr->addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),1,0.3)));
	TEST_EQUAL(ptr->getPopulationSize(),1);
}
END_SECTION

START_SECTION((int addIndividuals(std::vector< boost::shared_ptr< Chromosome > > individuals)))
{
	std::vector<boost::shared_ptr<Chromosome> > ngp;
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),1,0.4)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),1,0.2)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),1,0.5)));
	ngp.push_back(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),1,0.1)));
	ptr->addIndividuals(ngp);
	TEST_EQUAL(ptr->getPopulationSize(),5);
}
END_SECTION

START_SECTION((void sort(const int sortMethod=Chromosome::sortScoreDescending)))
{
  ptr->sort();
  TEST_STRING_EQUAL(ptr->getIndividual(0)->getSequence().toString(),"QQPEHRGEIIPVV");
  ptr->sort(Chromosome::sortScoreAscending);
  TEST_STRING_EQUAL(ptr->getIndividual(0)->getSequence().toString(),"CNGTQQFETMIND");
}
END_SECTION

START_SECTION((const boost::shared_ptr<Chromosome> getIndividual(const Size which)))
{
	ptr->sort();
	TEST_STRING_EQUAL(ptr->getIndividual(4)->getSequence().toString(),"CNGTQQFETMIND");
}
END_SECTION


START_SECTION((std::vector<boost::shared_ptr<Chromosome> >::iterator begin()))
{
	TEST_STRING_EQUAL((*ptr->begin())->getSequence().toString(),"QQPEHRGEIIPVV");
}
END_SECTION

START_SECTION((std::vector<boost::shared_ptr<Chromosome> >::iterator end()))
{
	std::vector<boost::shared_ptr<Chromosome> >::iterator iter = ptr->begin();
	while(iter++ != ptr->end());
	//TEST_EQUAL((*(--iter))->getSequence().toString(),"CNGTQQFETMIND");
	//TODO find out how to test end().
}
END_SECTION

START_SECTION((void seed(const Size seed)))
{
	ptr->seed(100);
}
END_SECTION

START_SECTION((const std::vector<boost::shared_ptr<Chromosome> >& getGenPool() const ))
{
	std::vector<boost::shared_ptr<Chromosome> > cp = ptr->getGenPool();
	TEST_EQUAL(cp.size(),5);
	TEST_STRING_EQUAL(cp[0]->getSequence().toString(),"QQPEHRGEIIPVV");
}
END_SECTION

START_SECTION((void setGenPool(std::vector< boost::shared_ptr< Chromosome > > genPool)))
{
	std::vector<boost::shared_ptr<Chromosome> > cp = ptr->getGenPool();
	std::vector<boost::shared_ptr<Chromosome> > ev;
	ptr->setGenPool(ev);
	TEST_EQUAL(ptr->getPopulationSize(),0);
	ptr->setGenPool(cp);
	TEST_EQUAL(ptr->getPopulationSize(),5);
}
END_SECTION

START_SECTION((const std::map<String, boost::shared_ptr<Chromosome> > getKnownIndividuals() const ))
{
	std::map<String, boost::shared_ptr<Chromosome> > ki = ptr->getKnownIndividuals();
	TEST_EQUAL(ki.size(),5);
	//TEST_NOT_EQUAL(*ki.find("QQPEHRGEIIPVV"),*ki.end());
}
END_SECTION

START_SECTION((void setKnownIndividuals(const std::map< String, boost::shared_ptr< Chromosome > > &knownIndividuals)))
{
	std::map<String, boost::shared_ptr<Chromosome> > ki = ptr->getKnownIndividuals();
	std::map<String, boost::shared_ptr<Chromosome> > ni;
	ptr->setKnownIndividuals(ni);
	TEST_EQUAL(ptr->getKnownIndividuals().size(),0);
	ptr->setKnownIndividuals(ki);
	TEST_EQUAL(ptr->getKnownIndividuals().size(),5);
}
END_SECTION

START_SECTION((unsigned int getPreviousPoolSize() const ))
{
	TEST_EQUAL(ptr->getPreviousPoolSize(),0);
}
END_SECTION

START_SECTION((void setPreviousPoolSize(unsigned int previousPoolSize)))
{
	ptr->setPreviousPoolSize(20);
	TEST_EQUAL(ptr->getPreviousPoolSize(),20);
}
END_SECTION

START_SECTION((unsigned int getMaxPoolSize() const))
{
	TEST_EQUAL(ptr->getMaxPoolSize(),10);
}
END_SECTION

START_SECTION((void replenish(const int targetSize)))
{
	std::vector<boost::shared_ptr<Chromosome> > ev;
	ptr->setGenPool(ev);
	TEST_EQUAL(ptr->getPopulationSize(),0);
	ptr->replenish(5);
	TEST_EQUAL(ptr->getPopulationSize(),5);
}
END_SECTION

START_SECTION(~GenPool())
{
	delete ptr;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



