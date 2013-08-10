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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SEEDER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SEEDER_H

#include <OpenMS/config.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
	/**
	* @brief The Seeder base class provides a framework for
	* derived classes to implement the seeding differently.
	* Seeder may not be instantiated and createIndividual must be implemented.
	*/
  class OPENMS_DLLAPI Seeder
  {
private:
	/// Seed to initialize random processes.
	unsigned int randomSeed_;
	/// A list of amino acids that can form a sequence is needed for some mutation processes.
	std::vector<const Residue*> aaList_;
	/// it is necessary to know the precursor mass to propose suitable sequences
    double precursorMass_;
    /// it is necessary to know the precursor mass to propose suitable sequences
    double precursorMassTolerance_;
	/// it is necessary to know the fragment mass tolerance for some derived classes.
	double fragmentMassTolerance_;
	/// Mass spectrum which may be used in some deriving classes.
	const MSSpectrum<> * msms_;

private:
	/// To forbid copy construction
	Seeder(const Seeder& other);
	/// To forbid assignment
	Seeder & operator=(const Seeder& rhs);

public:
    /// Default c'tor accepting all necessary input
    Seeder(const MSSpectrum<> * spec, const double precursorMass, const double precursorMassTolerance, const double fagmentMassTolerance, const std::vector<const Residue*> aaList);

	/// d'tor
    virtual ~Seeder();

    /// Creates a new individual and must be implemented in derived classes.
    virtual boost::shared_ptr<Chromosome> createIndividual() const = 0;
	
    /// Creates new individuals and can be overridden in derived classes.
	/// tries to create as many as num elements but gives up after 20 * num tries and returns what it was able to create.
    virtual std::vector<boost::shared_ptr<Chromosome> > createIndividuals(const Size num) const;

	/// to change the seed or to fix it for unit tests.
	void seed(const unsigned int seed);

	/// Allows the retrieval of the currently set seed (needed in some derived classes).
	unsigned int getSeed() const
	{
		return randomSeed_;
	}
	
	/// Allows to retrieve the amino acid list (needed in derived classes).
	const std::vector<const Residue*>& getAAList() const {
		return aaList_;
	}
	
	/// Allows to retrieve the set precursor mass (needed in derived classes).
	double getPrecursorMass() const {
		return precursorMass_;
	}
	
	/// Allows to retrieve the set precursor mass tolerance (needed in derived classes)
	double getPrecursorMassTolerance() const {
		return precursorMassTolerance_;
	}

	/// Allows to retrieve the set precursor mass tolerance (needed in derived classes)
	double getFragmentMassTolerance() const {
		return fragmentMassTolerance_;
	}

	/// Allows to retrieve the spectrum (needed in derived classes)
	const MSSpectrum<> * getMassSpectrum() const {
		return msms_;
	}
};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SEEDER_H
