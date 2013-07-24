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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_MUTATER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_MUTATER_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>

namespace OpenMS
{
  class GenPool;

  class OPENMS_DLLAPI Mutater
  {
private:
	/// Seed to initialize rand
	unsigned int randomSeed_;
	/// Mutation rate (0-1) defines how many percent of the sequences will be mutated in each attempt (def: 0.2).
	double mutationRate_;
	/// For some mutation processes it is necessary to know the precursor mass.
    double precursorMass_;
    /// For some mutation processes it is necessary to know the mass tolerance.
    double precursorMassTolerance_;
    /// A list of amino acids that can form a sequence is needed for some mutation processes.
    std::vector<const Residue*> aaList_;

private:
    /// To forbid copy construction
    Mutater(const Mutater & rhs);
    /// To forbid assignment
    Mutater & operator=(const Mutater & rhs);

public:
    /// Default c'tor
    Mutater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

    virtual ~Mutater();

    /// Given a GenPool iterates over the members (Chromosomes) and given the mutationRate randomly selects Chromosomes to be mutated using mutate(Chromosome).
    void mutate(GenPool& genPool);

    /// The Chromosome that is submitted to this function will be mutated and changed.
    virtual void mutate(Chromosome& chromosome) = 0;

    /// The Chromosome submitted to this function will be copied
    /// and changes will be performed on the copy which will then be returned.
    Chromosome mutate(const Chromosome & chromosome) const;

	double getMutationRate() const {
		return mutationRate_;
	}

	/// Adjusts the mutation rate (def: 0.2)
	void setMutationRate(double mutationRate = 0.2)
	{
		OPENMS_PRECONDITION(mutationRate_ <=0 || mutationRate_ > 1, "Value must be between 0 (exclusive) and 1 (inclusive)")
		this->mutationRate_ = mutationRate;
	}

	/// to change the seed or to fix it for unit tests.
	void seed(const unsigned int seed);

	const std::vector<const Residue*>& getAaList() const {
		return aaList_;
	}

	void setAaList(const std::vector<const Residue*>& aaList) {
		aaList_ = aaList;
	}

	double getPrecursorMass() const {
		return precursorMass_;
	}

	void setPrecursorMass(double precursorMass) {
		precursorMass_ = precursorMass;
	}

	double getPrecursorMassTolerance() const {
		return precursorMassTolerance_;
	}

	void setPrecursorMassTolerance(double precursorMassTolerance) {
		precursorMassTolerance_ = precursorMassTolerance;
	}

	unsigned int getRandomSeed() const {
		return randomSeed_;
	}
};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_MUTATER_H
