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
// $Maintainer: Jens Allmer $
// $Authors: Jens Allmer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <stdexcept>
#include <algorithm>

namespace OpenMS
{

  GenPool::GenPool() :
    maxPoolSize(200), precursorMassTolerance(0.3)
  {
	  // mutater = new DefaultMutater(precursorMass,precursorMassTolerance,aaList);
  }

  void GenPool::setMutater(const Mutater& mutater)
  {
  }

  void GenPool::sort(const int sortMethod = Chromosome::sortScoreDescending)
  {
	  switch(sortMethod) {
	    case Chromosome::sortScoreAscending :
	      std::sort(*genPool.begin(), *genPool.end(), Chromosome::sortScoreAsc);
	      break;
	    case Chromosome::sortScoreDescending :
		  std::sort(*genPool.begin(), *genPool.end(), Chromosome::sortScoreDesc);
		  break;
	  }
  }

  Size GenPool::getPopulationSize() {
    return(genPool.size());
  }

  void GenPool::setPool(std::vector<Chromosome *> newPool)
  {
    genPool = newPool;
  }

  void GenPool::initGenPool(const int maxPoolSize)
  {
	Size maxTries = 10 * maxPoolSize;
    while((maxTries > 0) && (genPool.size() < maxPoolSize))
    {
      Chromosome ni = seeder->createIndividual();
      addIndividual(ni);
      maxTries--;
    }
  }

	void GenPool::setMater(const Mater& mater)
	{
	}

	void GenPool::setKiller(const Killer& killer)
	{
	}

	void GenPool::setSeeder(const Seeder& seeder)
	{
	  this->seeder = &seeder;
	}

  bool GenPool::addIndividual(Chromosome individual) {
	  try
	  {
		Chromosome *known = knownIndividuals.at(individual.getSequence().toString());
	  }
	  catch(const std::out_of_range& oor)
	  {
		genPool.push_back(&individual);
		knownIndividuals.insert(std::pair<String,Chromosome *>(individual.getSequence().toString(),&individual));
		return true;
	  }
	  return false;
  }
} // namespace
