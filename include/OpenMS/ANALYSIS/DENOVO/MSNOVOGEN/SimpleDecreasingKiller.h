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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_SIMPLEDECREASINGKILLER_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_SIMPLEDECREASINGKILLER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Killer.h>

namespace OpenMS
{
  /**
  * @brief The SimpleDecreasingKiller aims to decrease the population size by
  * a set number of individuals in each generation (if possible).
  * If less individuals exist than should be killed, less are being killed to account for that.
  * The last individual always survives.
  */
  class OPENMS_DLLAPI SimpleDecreasingKiller : public Killer
  {
private:
	  /// By how many absolute members a population (genPool) should shrink per generation (default: 5).
	  unsigned int decreasePerGeneration;

public:
    /// Default c'tor
    SimpleDecreasingKiller();

    /// Implementation of virtual ~Killer()
    ~SimpleDecreasingKiller();

    /// Implementing virtual function from Killer to actually delete Chromosomes from the current genPool.
    void kill(GenPool& genPool) const;

    /// Returns how many individuals should be killed per generation.
	unsigned int getDecreasePerGeneration() const
	{
		return decreasePerGeneration;
	}

	/// Allows setting of how many individuals should be killed per generation.
	void setDecreasePerGeneration(unsigned int decreasePerGeneration)
	{
		this->decreasePerGeneration = decreasePerGeneration;
	}
};
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_SIMPLEDECREASINGKILLER_H
