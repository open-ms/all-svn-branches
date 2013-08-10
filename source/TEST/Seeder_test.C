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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

///////////////////////////

using namespace OpenMS;
using namespace std;
struct TestSeeder :
	public Seeder
{
		TestSeeder(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList) :
			Seeder(precursorMass, precursorMassTolerance, aaList)
		{}

	    ~TestSeeder()
	    {}

	    boost::shared_ptr<Chromosome> createIndividual() const
	    {
	    	return(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("ALLMER"),1,1.0)));
	    }
};
START_TEST(Seeder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

double precursorMass=1500.5;
double precursorMassTolerance=2;
std::vector<const Residue*> aaList;

Seeder* ptr = 0;
Seeder* null_ptr = 0;

START_SECTION((Seeder(double precursorMass, double precursorMassTolerance, std::vector< const Residue * > aaList)))
{

	ptr = new TestSeeder(precursorMass,precursorMassTolerance,aaList);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual boost::shared_ptr<Chromosome> createIndividual() const = 0))
{
	boost::shared_ptr<Chromosome> chr(ptr->createIndividual());
	TEST_STRING_EQUAL(chr->getSequence().toString(),"ALLMER");
}
END_SECTION

START_SECTION((void seed(const unsigned int seed)))
{
  TEST_NOT_EQUAL(ptr->getSeed(),50);
  ptr->seed(50);
  TEST_EQUAL(ptr->getSeed(),50);
}
END_SECTION
	
START_SECTION((const std::vector<const Residue*>& getAAList() const))
{
  TEST_EQUAL(ptr->getAAList().size(),0);
}
END_SECTION
	
START_SECTION((double getPrecursorMass() const))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMass(),1500.5);
}
END_SECTION

START_SECTION((double getPrecursorMassTolerance() const))
{
  TEST_REAL_SIMILAR(ptr->getPrecursorMassTolerance(),2);
}
END_SECTION

START_SECTION((unsigned int getSeed() const))
{
  TEST_EQUAL(ptr->getSeed(),50);
}
END_SECTION

START_SECTION(~Seeder())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



