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

#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/HomologyKiller.h>
#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>
#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{

  HomologyKiller::HomologyKiller(const int maxPopulation, const int initialPopulation) :
    Killer(maxPopulation, initialPopulation), minEditDistance(2)
  {}

  void HomologyKiller::kill(GenPool& genPool)
  {
	  genPool.sort(Chromosome::sortScoreDescending);
	  std::vector<boost::shared_ptr<Chromosome> > ngp;
	  std::vector<boost::shared_ptr<Chromosome> >::iterator i = genPool.begin();
	  boost::shared_ptr<Chromosome> b = * i;
	  ngp.push_back(b);
	  b->getSequence().toString();
	  while(i != genPool.end())
	  {
		boost::shared_ptr<Chromosome> o = *i;
		boost::shared_ptr<Chromosome> p;
		bool homolog = false;
		for(std::vector<boost::shared_ptr<Chromosome> >::iterator n = ngp.begin(); n != ngp.end(); n++) {
			p = *n;
			if(Utilities::editDistance(p->getSequence(),o->getSequence()) < minEditDistance)
			{
			  homolog = true;
			  break;
			}
		}
		if(!homolog)
		  ngp.push_back(p);
		i++;
	  }
	  genPool.setPool(ngp);
	  genPool.replenish(Killer::getPreviousPopulation());
  }

} // namespace
