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
#include <time.h>
#include <stdlib.h>
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h"
#include <vector>

using std::vector;

namespace OpenMS
{

  GenPool::GenPool(const int maxPoolSize, const double pmt, const double pm, const int charge) :
    maxPoolSize_(maxPoolSize), previousPoolSize_(0),
    precursorMassTolerance_(pmt),
    precursorMass_(pm), precursorCharge_(charge)
  {
	  seed(time(0));
  }

  void GenPool::sort(const int sortMethod)
  {
	  switch(sortMethod) {
	    case Chromosome::sortScoreAscending :
	      std::sort(genPool_.begin(), genPool_.end(), Chromosome::sortScoreAsc);
	      break;
	    case Chromosome::sortScoreDescending :
		  std::sort(genPool_.begin(), genPool_.end(), Chromosome::sortScoreDesc);
		  break;
	  }
  }

  Size GenPool::getPopulationSize() {
    return(genPool_.size());
  }

  bool GenPool::addIndividual(boost::shared_ptr<Chromosome> individual)
  {
	  if (knownIndividuals_.find(individual->getSequence().toString()) == knownIndividuals_.end())
	  {
		if(std::abs(precursorMass_ - individual->getSequence().getMonoWeight(Residue::Full)) < precursorMassTolerance_)
		{
			individual->setCharge(precursorCharge_);
			genPool_.push_back(individual);
			knownIndividuals_.insert(std::pair<String,boost::shared_ptr<Chromosome> >(individual->getSequence().toString(),individual));
			return true;
		}
	  }
	  return false;
  }

  int GenPool::addIndividuals(std::vector<boost::shared_ptr<Chromosome> > individuals)
  {
	int added = 0;
    for(std::vector<boost::shared_ptr<Chromosome> >::iterator i = individuals.begin(); i != individuals.end(); i++)
    {
    	if(addIndividual(*i))
    		added++;
    }
    return added;
  }

  const boost::shared_ptr<Chromosome> GenPool::getIndividual(const Size which)
  {
    if(which < genPool_.size())
      return genPool_[which];
    boost::shared_ptr<Chromosome> ret(new Chromosome(AASequence(""),0));
    return(ret);
  }



  void GenPool::replenish(const int targetSize)
  {
	unsigned int ts = targetSize;
	if(ts > previousPoolSize_)
		ts = previousPoolSize_;
    int toAdd = ts - getPopulationSize();
    if(toAdd <= 0)
    	return;
    std::vector<boost::shared_ptr<Chromosome> > ki;
    if(knownIndividuals_.size() > 0)
    {
      for(std::map<String, boost::shared_ptr<Chromosome> >::iterator i = knownIndividuals_.begin(); i != knownIndividuals_.end(); i++) {
    	 ki.push_back(i->second);
      }
    }
    for(int i=0; i<toAdd; i++)
    {
    	Size pos = rand() % ki.size();
    	genPool_.push_back(ki[pos]);
    }
    sort(Chromosome::sortScoreDescending);
  }

  void GenPool::setGenPool(std::vector<boost::shared_ptr<Chromosome> > genPool)
  {
	  genPool_ = genPool;
  }

	void GenPool::addToGenPool(std::vector<boost::shared_ptr<Chromosome> > genPool)
	{
		addIndividuals(genPool);
	}

} // namespace
