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

#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/ZipMater.h>
#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>

namespace OpenMS
{

  ZipMater::ZipMater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList) :
    Mater(precursorMass,precursorMassTolerance, aaList)
  {}

  std::vector<boost::shared_ptr<Chromosome> > ZipMater::mate(boost::shared_ptr<Chromosome> lhs, boost::shared_ptr<Chromosome> rhs) const
  {
    String uc("");
    String lc("");
    Size i = 0;
    bool swap = false;
    while(i < lhs->getSequence().size() && i < rhs->getSequence().size())
    {
    	if(swap)
    	{
			uc += rhs->getSequence().getResidue(i).getModifiedOneLetterCode();
			lc += lhs->getSequence().getResidue(i).getModifiedOneLetterCode();
    	}
    	else
    	{
			uc += lhs->getSequence().getResidue(i).getModifiedOneLetterCode();
			lc += rhs->getSequence().getResidue(i).getModifiedOneLetterCode();
    	}
    	swap = !swap;
    	++i;
    }
    // since either lhs or rhs ran out of aas the remainder of the longer sequence can be appended
    // without test which on is actually longer.
    while(i < lhs->getSequence().size())
    {
    	lc += lhs->getSequence().getResidue(i).getModifiedOneLetterCode();
    	++i;
    }
    while(i < lhs->getSequence().size())
	{
		uc += rhs->getSequence().getResidue(i).getModifiedOneLetterCode();
    	++i;
	}
    std::vector<boost::shared_ptr<Chromosome> > ret;
    AASequence ucaa(uc);
    if(Utilities::adjustToFitMass(getSeed(),ucaa,getPrecursorMass(),getPrecursorMassTolerance(),getAAList()))
    	ret.push_back(boost::shared_ptr<Chromosome>(new Chromosome(ucaa,lhs->getCharge())));
    AASequence lcaa(lc);
    if(Utilities::adjustToFitMass(getSeed(),lcaa,getPrecursorMass(), getPrecursorMassTolerance(), getAAList()))
		ret.push_back(boost::shared_ptr<Chromosome>(new Chromosome(lcaa,lhs->getCharge())));
    return ret;
  }
} // namespace
