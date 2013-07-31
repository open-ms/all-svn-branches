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

		void mutate(boost::shared_ptr<Chromosome> chromosome)
	    {

	    }
};

START_TEST(Mutater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Mutater* ptr = 0;
Mutater* null_ptr = 0;
START_SECTION(Mutater())
{

}
END_SECTION

START_SECTION((Mutater(double precursorMass, double precursorMassTolerance, std::vector< Residue & > aaList)))
{
	std::vector< const Residue * > aaList;
	ptr = new TestMutater(0.1, 0.2, aaList);
	TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(~Mutater())
{
	delete ptr;
}
END_SECTION

START_SECTION((Mutater(const Mutater &other)))
{
  // TODO
}
END_SECTION

START_SECTION((Mutater& operator=(const Mutater &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((double getMutationRate() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setMutationRate(double mutationRate=0.2)))
{
  // TODO
}
END_SECTION

START_SECTION((GenPool& mutate(GenPool &genPool)))
{
  // TODO
}
END_SECTION

START_SECTION((Chromosome& mutate(Chromosome &chromosome)=0))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



