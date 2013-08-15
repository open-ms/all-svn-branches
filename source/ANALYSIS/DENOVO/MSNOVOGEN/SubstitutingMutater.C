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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SubstitutingMutater.h>
#include <stdlib.h>
#include <time.h>
#include <boost/random/uniform_int_distribution.hpp>

namespace OpenMS
{

  SubstitutingMutater::SubstitutingMutater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList)
    : Mutater(precursorMass, precursorMassTolerance, aaList)
  {}

  void SubstitutingMutater::mutate(boost::shared_ptr<Chromosome> chromosome) const
  {
	AASequence as = chromosome->getSequence();
	double seqMass = as.getMonoWeight();
	boost::random::uniform_int_distribution<Size> seqPos(0,(as.size()-1));
	Size rp = seqPos(rng_);
    boost::random::uniform_int_distribution<Size> listPos(0, (getAAList().size()-1));
    int sa = listPos(rng_);
    const Residue * nr = getAAList()[sa];
    const Residue & replaced = as.setResidue(rp, nr);

    std::vector<std::pair<Size,const Residue*> > possRep;
    double mdiff = replaced.getMonoWeight(Residue::Full) - nr->getMonoWeight(Residue::Full);
    for(Size i = 0; i < as.size(); i++)
    {
    	if(i == rp)
    	  continue;
    	for(Size p = 0; p < getAAList().size(); p++)
    	{
		  if(p == rp)
		    continue;
    	  const Residue& cu = as.getResidue(i);
		  const Residue* pr = getAAList()[p];
		  if(pr->getOneLetterCode() == replaced.getOneLetterCode())
			  continue;
    	  double pdiff = cu.getMonoWeight(Residue::Full) - pr->getMonoWeight(Residue::Full);
    	  if(std::abs(getPrecursorMass() - (seqMass - mdiff - pdiff)) <= getPrecursorMassTolerance())
    	  {
    		std::pair<Size, const Residue*> p(i,pr);
    	    possRep.push_back(p);
    	  }
    	}
    }
    if(possRep.size() > 0) {
	  boost::random::uniform_int_distribution<Size> posDist(0, (possRep.size()-1));
	  Size w = posDist(rng_);
      as.setResidue(possRep[w].first, possRep[w].second);
      chromosome->setSequence(as);
    }
  }

} // namespace
