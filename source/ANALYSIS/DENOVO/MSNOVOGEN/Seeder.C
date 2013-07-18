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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <time.h>
#include <stdlib.h>

namespace OpenMS
{

  Seeder::Seeder(double pm, double pmt, std::vector<const Residue*> al) :
    precursorMass(pm), precursorMassTolerance(pmt), aaList(al)
  {
	seed(time(0));
  }

  Seeder::~Seeder()
  {
  }

  void Seeder::seed(const unsigned int seed)
  {
	randomSeed = seed;
	srand(randomSeed);
  }

  const AASequence Seeder::getRandomSequence(const int len, const double weight, const double tolerance)
  {
    String str;
    for(int i=0; i<len; i++)
      str += getRandomAA();
    AASequence seq(str);
    adjustToFitMass(seq, weight, tolerance);
    return(seq);
  }

  bool Seeder::adjustToFitMass(AASequence& sequence, const double weight, const double tolerance) const
  {
	double curWeight = sequence.getMonoWeight();
	double diff = std::abs(curWeight-weight);
	double nDiff;
	double minDiff;
	Size pos;
	if(diff <= tolerance)
	  return true;
	int mi = 10;	//max iterations for while loop
	while(diff > tolerance)
    {
	  pos = rand() % sequence.size();
	  minDiff = diff;
	  const Residue * replace;
	  for(Size lp=0; lp<aaList.size(); lp++)
	  {
	    nDiff = std::abs((curWeight - sequence[pos].getMonoWeight(Residue::Full) + aaList[lp]->getMonoWeight(Residue::Full)) - weight);
	    if(nDiff < minDiff)
	    {
		  minDiff = nDiff;
		  replace = aaList[lp];
	    }
	  }
	  sequence.setResidue(pos,replace);
	  curWeight = sequence.getMonoWeight(Residue::Full);
	  diff = std::abs(curWeight - weight);
	  if(--mi == 0)
		  return(false);
    }
	return(true);
  }

  const String Seeder::getRandomAA() const
  {
    Size i = rand() % aaList.size();
    return(aaList[i]->getModifiedOneLetterCode());
  }

} // namespace
