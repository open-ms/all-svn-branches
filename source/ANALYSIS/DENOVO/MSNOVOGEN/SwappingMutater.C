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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SwappingMutater.h>
#include <stdlib.h>
#include <boost/random/uniform_int_distribution.hpp>

namespace OpenMS
{

  SwappingMutater::SwappingMutater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList) :
    Mutater(precursorMass,precursorMassTolerance,aaList)
  {
	  seed((unsigned int)time(0));
  }

  void SwappingMutater::mutate(boost::shared_ptr<Chromosome> chromosome) const
  {
    AASequence as = chromosome->getSequence();
    Size mid = as.size()/2;
	boost::random::uniform_int_distribution<> ldist(0, (int)(mid-1));
	boost::random::uniform_int_distribution<> rdist(0, (int)((mid+mid)-1));
	Size lpos = ldist(rng);
    Size rpos = rdist(rng);
    const Residue & r = as.getResidue(rpos);
    const Residue & nr = as.setResidue(lpos, &r);
    as.setResidue(rpos, &nr);
    chromosome->setSequence(as);
  }

  SwappingMutater::~SwappingMutater()
  {}

} // namespace
