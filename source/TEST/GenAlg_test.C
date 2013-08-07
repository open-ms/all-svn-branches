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
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GenAlg, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GenAlg* ptr = 0;
GenAlg* null_ptr = 0;
START_SECTION(GenAlg())
{
	ptr = new GenAlg();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~GenAlg())
{
	delete ptr;
}
END_SECTION

START_SECTION((GenAlg(const MSSpectrum<> *spec, const std::vector< const Residue * > &aaList, const int maxPoolSize, const double precursorMassTolerance, const double fragmentMassTolerance)))
{
  // TODO
}
END_SECTION

START_SECTION((~GenAlg()))
{
  // TODO
}
END_SECTION

START_SECTION((void startEvolution(const int numGenerations)))
{
  // TODO
}
END_SECTION

START_SECTION((void startEvolution(const int numGenerations, const int endIfStableForNumGenerations)))
{
  // TODO
}
END_SECTION

START_SECTION((void setMutater(boost::shared_ptr< const Mutater > mutater)))
{
  // TODO
}
END_SECTION

START_SECTION((void setMater(boost::shared_ptr< const Mater > mater)))
{
  // TODO
}
END_SECTION

START_SECTION((void setKiller(boost::shared_ptr< const Killer > killer)))
{
  // TODO
}
END_SECTION

START_SECTION((void setSeeder(boost::shared_ptr< const Seeder > seeder)))
{
  // TODO
}
END_SECTION

START_SECTION((void setScorer(boost::shared_ptr< const Scorer > scorer)))
{
  // TODO
}
END_SECTION

START_SECTION((void setHomologyKiller(boost::shared_ptr< const HomologyKiller > killer)))
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

START_SECTION((const std::map<String, boost::shared_ptr<Chromosome> > getKnownIndividuals() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setKnownIndividuals(const std::map< String, boost::shared_ptr< Chromosome > > &knownIndividuals)))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorIntensity() const ))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorMH() const ))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorMass() const ))
{
  // TODO
}
END_SECTION

START_SECTION((double getPrecursorMz() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRejuvenate(const unsigned int afterNGen)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



