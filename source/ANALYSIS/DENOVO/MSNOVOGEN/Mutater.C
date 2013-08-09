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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <stdlib.h>
#include <time.h>

namespace  OpenMS
{

  Mutater::Mutater(double pm, double pmt, std::vector<const Residue*> al)
    : mutationRate_(0.2), precursorMass_(pm), precursorMassTolerance_(pmt), aaList_(al)
  {
	//mutationRate = 0.2;
    //this->precursorMass = precursorMass;
    //this->precursorMassTolerance = precursorMassTolerance;
    //this->aaList = aaList;

    // make sure random is initialized
	  seed((unsigned int)time(0));
  }

  Mutater::~Mutater()
  {
  }

	void Mutater::mutatePool(GenPool& pool) const
	{
		for(std::vector<boost::shared_ptr<Chromosome> >::iterator iter = pool.begin(); iter!= pool.end(); ++iter)
		{
		  double rv = (rand() % 101)/100.0;
		  if(rv > getMutationRate())
			continue;
		  this->mutate(*iter);
		}
	}

	void Mutater::mutateAndAddToPool(GenPool& pool) const
	{
		std::vector<boost::shared_ptr<Chromosome> > add;
		for(std::vector<boost::shared_ptr<Chromosome> >::iterator iter = pool.begin(); iter!= pool.end(); ++iter)
		{
		  double rv = (rand() % 101)/100.0;
		  if(rv > getMutationRate())
			continue;
		  boost::shared_ptr<Chromosome> ni(new Chromosome((*iter)->getSequence(),(*iter)->getCharge()));
		  this->mutate(ni);
		  add.push_back(ni);
		}
		pool.addIndividuals(add);
	}

	boost::shared_ptr<Chromosome> Mutater::mutateCpy(const boost::shared_ptr<Chromosome> chromosome) const
	{
		boost::shared_ptr<Chromosome> cpy(new Chromosome(AASequence(chromosome->getSequence()),chromosome->getScore()));
		this->mutate(cpy);
		return cpy;
	}

  void Mutater::seed(unsigned int seed)
  {
	randomSeed_ = seed;
	srand(randomSeed_);
  }

} // namespace
