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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_HyperScorer_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_HyperScorer_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Scorer.h>
#include <vector>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
	/**
	* @brief The HyperScorer applys the X!Tandem hyperscore calculation.
	*/
  class OPENMS_DLLAPI HyperScorer : public Scorer
  {
private:
	/// Calculates faculty of i recursively.
	Size fac(const Size i) const;

	/// Assignment operator shouldn't be used
	HyperScorer & operator=(const HyperScorer& rhs);

	/// Copy c'tor shouldn't be used
	HyperScorer(const HyperScorer& other);

public:
	/// Default c'tor providing all necessary information.
	HyperScorer(const double fragmentMassTolerance);

	/// Implementation of virtual method Scorer::score.
    void score(const MSSpectrum<> * msms, boost::shared_ptr<Chromosome> & chromosome) const;

  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_HyperScorer_H
