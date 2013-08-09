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

#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/RandomMater.h>
#include <stdlib.h>
#include <time.h>

namespace OpenMS
{

  RandomMater::RandomMater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList) :
    Mater(precursorMass, precursorMassTolerance, aaList),
    sm(precursorMass, precursorMassTolerance, aaList),
    dm(precursorMass, precursorMassTolerance, aaList),
    zm(precursorMass, precursorMassTolerance, aaList)
  {
	weights_.push_back(0.4);
	weights_.push_back(0.7);
	weights_.push_back(1.0);
	seed((unsigned int)time(0));
  }

  std::vector<boost::shared_ptr<Chromosome> > RandomMater::mate(boost::shared_ptr<Chromosome> lhs, boost::shared_ptr<Chromosome> rhs) const
  {
	double rv = (double)(rand() % 101) / (double)100; //Random number between 0 and 1 (inclusive).
    for(unsigned int i=0; i<weights_.size(); i++)
    {
      if(weights_[i] <= rv)
      {
        switch(i) {
          case RandomMater::simpleMater :
			   sm.mate(lhs, rhs);
			   break;
          case RandomMater::zipMater :
			   zm.mate(lhs, rhs);
			   break;
          default:
			   dm.mate(lhs, rhs);
			   break;
        }
      }
    }
    return dm.mate(lhs, rhs);
  }

  void RandomMater::seed(unsigned int seed)
  {
	  sm.seed(seed);
	  zm.seed(seed);
	  dm.seed(seed);
  }
} // namespace
