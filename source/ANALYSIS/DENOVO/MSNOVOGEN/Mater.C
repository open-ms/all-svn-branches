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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>
#include <time.h>
#include <stdlib.h>


namespace OpenMS
{
	Mater::Mater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList) :
	aaList_(aaList), precursorMass_(precursorMass), precursorMassTolerance_(precursorMassTolerance)
	{
	  seed((unsigned int)time(0));
	}

	Mater::~Mater()
	{}

	boost::shared_ptr<Chromosome> Mater::getPartner(GenPool& genPool, boost::shared_ptr<Chromosome> exclude) const
	{
		boost::shared_ptr<Chromosome> ret;
		boost::shared_ptr<Chromosome> test;
		Size rv = 0;
		int dist = 0;
		while(!ret)
		{
			rv = rand() % genPool.getPopulationSize();
			test = genPool.getGenPool()[rv];
			dist = Utilities::editDistance(exclude->getSequence(),test->getSequence());
			if(dist < 3)
				continue;
			ret = test;
		}
		return(ret);
	}

	void Mater::seed(const unsigned int seed)
	{
		randomSeed_ = seed;
		srand(randomSeed_);
	}

	std::vector<boost::shared_ptr<Chromosome> > Mater::tournament(GenPool & pool) const
	{
		std::vector<boost::shared_ptr<Chromosome> > ret;
		for(std::vector<boost::shared_ptr<Chromosome> >::iterator iter = pool.begin(); iter!= pool.end(); ++iter)
		{
		  boost::shared_ptr<Chromosome> p1 = this->getPartner(pool,*iter);
		  if(!p1)
			  continue;
		  boost::shared_ptr<Chromosome> p2 = *iter;
		  std::vector<boost::shared_ptr<Chromosome> > children = this->mate(p1, p2);
		  for(std::vector<boost::shared_ptr<Chromosome> >::iterator c = children.begin(); c != children.end(); ++c)
			  ret.push_back(*c);
		}
		return ret;
	}

	void Mater::tournamentAndAddToPool(GenPool & pool) const
	{
		std::vector<boost::shared_ptr<Chromosome> > children = tournament(pool);
		pool.addIndividuals(children);
	}
} // namespace
