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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SCORER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SCORER_H

#include <OpenMS/config.h>
#include <boost/shared_ptr.hpp>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>

namespace OpenMS
{
	/**
	* @brief The Scorer base class provides a framework for
	* derived classes to implement the scoring differently.
	* Scorer may not be instantiated and kill must be implemented.
	*/
  class GenPool;

  class OPENMS_DLLAPI Scorer
  {
private:
	/// the fragment mass tolerance to use when searching for agreement among spectra.
	double fragmentMassTolerance;

	/// Copy c'tor shouldn't be used.
	Scorer(const Scorer& other);

	/// Assignment operator shouldn't be used.
	Scorer & operator=(const Scorer& rhs);

public:
    /// Default c'tor accepting all necessary parameters.
    Scorer(const double fragmentMassTolerance = 0.5);

	/// d'tor.
    virtual ~Scorer();

	/// Implemented in the base class.
	/// Calls score for each individual in the pool that hasn't been scored.
    void scorePool(const MSSpectrum<> * msms, GenPool & pool) const;

	/// scoring method must be implemented in deriving classes.
    virtual void score(const MSSpectrum<> * msms, boost::shared_ptr<Chromosome> & chromosome) const = 0;

	/// Allows to get the currently set fragment mass tolerance (needed in derived classes).
	double getFragmentMassTolerance() const {
		return fragmentMassTolerance;
	}
};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SCORER_H
