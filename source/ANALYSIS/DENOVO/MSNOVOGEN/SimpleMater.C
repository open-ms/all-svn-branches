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

#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/SimpleMater.h>
#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>
#include <stdlib.h>


namespace OpenMS
{

  SimpleMater::SimpleMater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList) :
    Mater(precursorMass, precursorMassTolerance, aaList)
  {}

  std::vector<boost::shared_ptr<Chromosome> > SimpleMater::mate(const boost::shared_ptr<Chromosome> lhs, const boost::shared_ptr<Chromosome> rhs) const
  {
	Size min = std::min(lhs->getSequence().size(), rhs->getSequence().size());
	std::vector<boost::shared_ptr<Chromosome> > ret;
	boost::random::uniform_int_distribution<int> int_distribution(0, ((int)(min-1)));
	int cop = int_distribution(rng);
	String uc(lhs->getSequence().getSubsequence(0,cop).toString() + rhs->getSequence().getSubsequence(cop,rhs->getSequence().size()-cop).toString());
	AASequence ucaa(uc);
	if(getUtils()->adjustToFitMass(ucaa,getPrecursorMass(),getPrecursorMassTolerance(),getAAList()))
	{
		ret.push_back(boost::shared_ptr<Chromosome>(new Chromosome(ucaa,lhs->getCharge())));
	}
	String lc(rhs->getSequence().getSubsequence(0,cop).toString() + lhs->getSequence().getSubsequence(cop,lhs->getSequence().size()-cop).toString());
	AASequence lcaa(lc);
	if(getUtils()->adjustToFitMass(lcaa,getPrecursorMass(),getPrecursorMassTolerance(),getAAList()))
		ret.push_back(boost::shared_ptr<Chromosome>(new Chromosome(lcaa,lhs->getCharge())));
	return ret;
  }

} // namespace
