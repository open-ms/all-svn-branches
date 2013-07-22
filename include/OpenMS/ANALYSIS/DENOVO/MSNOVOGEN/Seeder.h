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

namespace OpenMS
{
  class OPENMS_DLLAPI Seeder
  {
private:
	/// Seed to initialize rand
	unsigned int randomSeed;

protected:
	/// A list of amino acids that can form a sequence is needed for some mutation processes.
	std::vector<const Residue*> aaList;
	/// it is necessary to know the precursor mass to propose suitable sequences
    double precursorMass;
    /// it is necessary to know the precursor mass to propose suitable sequences
    double precursorMassTolerance;

private:
	/// To forbid copy construction
	Seeder(const Seeder& other);
	/// To forbid assignment
	Seeder & operator=(const Seeder& rhs);

public:
    /// Default c'tor
    Seeder(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

    virtual ~Seeder();

    /// Creates a new individual.
    virtual Chromosome createIndividual() const = 0;

	/// to change the seed or to fix it for unit tests.
	void seed(const unsigned int seed);

	/// Generates a random sequence of the desired length and then adjusts the
	/// amino acids until the weight matches the desired criteria.
	const AASequence getRandomSequence(const int len, const double weight, const double tolerance);

    /// from the list of available amino acids selects a random one.
    const String getRandomAA() const;

    /// The random string may not fit to the expected precursor mass and is then adjusted to fit (if possible with a maximum of 3 changes).
    /// The passed in sequence is directly modified.
    bool adjustToFitMass(AASequence & sequence, const double weight, const double tolerance) const;
  };
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SEEDER_H
