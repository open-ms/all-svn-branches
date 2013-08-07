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
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleDecreasingKiller.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <boost/shared_ptr.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SimpleDecreasingKiller, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SimpleDecreasingKiller* ptr = 0;
SimpleDecreasingKiller* null_ptr = 0;

GenPool gp(500,1000,1000,1);
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTA"),1,0.1)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTR"),1,0.15)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTN"),1,0.17)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTD"),1,0.2)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTC"),1,0.25)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTE"),1,0.27)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTQ"),1,0.3)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTG"),1,0.35)));
gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("TESTH"),1,0.37)));
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

START_SECTION(SimpleDecreasingKiller())
{
	ptr = new SimpleDecreasingKiller();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((~SimpleDecreasingKiller()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void kill(GenPool &genPool)))
{
  ptr->kill(gp);
  TEST_EQUAL(gp.getPopulationSize(),15);
}
END_SECTION

START_SECTION((int getDecreasePerGeneration() const ))
{
  TEST_EQUAL(ptr->getDecreasePerGeneration(),5);
}
END_SECTION

START_SECTION((void setDecreasePerGeneration(int decreasePerGeneration)))
{
	ptr->setDecreasePerGeneration(10);
	TEST_EQUAL(ptr->getDecreasePerGeneration(),10);
	ptr->kill(gp);
	TEST_EQUAL(gp.getPopulationSize(),5);
	ptr->kill(gp);
	TEST_EQUAL(gp.getPopulationSize(),1);
}
END_SECTION

START_SECTION(~SimpleDecreasingKiller())
{
	delete ptr;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



