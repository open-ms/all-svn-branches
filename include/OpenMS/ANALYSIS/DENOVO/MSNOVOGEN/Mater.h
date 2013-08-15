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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_MATER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_MATER_H

#include <OpenMS/config.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>

namespace OpenMS
{
	/**
	* @brief The Mater base class provides a framework for
	* derived classes to implement the mating differently.
	* MAter may not be instantiated and mate must be implemented in derived classes.
	*/
  class GenPool;

  class OPENMS_DLLAPI Mater
  {
protected:
	/// random number generator must be mutable since Mater must be const * in GenAlg.
	mutable boost::mt19937 rng_;
private:
	/// A list of amino acids that can form a sequence is needed for some mutation processes.
	std::vector<const Residue*> aaList_;
	/// it is necessary to know the precursor mass to propose suitable sequences
	double precursorMass_;
	/// it is necessary to know the precursor mass to propose suitable sequences
	double precursorMassTolerance_;
	/// Utilities class necessary for some calculations.
	mutable Utilities utils_;

private:
	/// Copy c'tor shouldn't be used.
	Mater(const Mater& other);

	/// Assignment operator shouldn't be used.
	Mater & operator=(const Mater& rhs);

	/// Randomly selects a mating partner, not using exclude, with an edit distance of at least 2.
	/// if no mating partner can be found returns a null pointer (can be tested by: if(!getPartner(..))).
	boost::shared_ptr<Chromosome> getPartner(GenPool & genPool, boost::shared_ptr<Chromosome> exclude) const;

public:
    /// Default c'tor providing all necessary parameters.
    Mater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

    /// d'tor
    virtual ~Mater();

    /// virtual mate method needs to be implemented in all deriving classes.
    /// This differentiates between differnet possibilities to perform crossover.
	/// May return an empty vector if mating was not possible or if no suitable child could be created.
    virtual std::vector<boost::shared_ptr<Chromosome> > mate(const boost::shared_ptr<Chromosome> lhs, const boost::shared_ptr<Chromosome> rhs) const = 0;

    /// This method performs mating for all individuals in a GenPool.
    /// It calls the deriving classes mate method and is implemented in this base class.
    std::vector<boost::shared_ptr<Chromosome> > tournament(GenPool & genPool) const;

    /// Same as tournament but instead of returning the new individuals directly adds them to the GenPool.
    void tournamentAndAddToPool(GenPool & genPool) const;

	/// to change the seed or to fix it for unit tests.
	void seed(const Size seed) const;

	/// Allows derived classes to access the amino acid list.
	const std::vector<const Residue*>& getAAList() const {
		return aaList_;
	}

	/// Allows derived classes to access the set precursor mass.
	double getPrecursorMass() const {
		return precursorMass_;
	}

	/// Allows derived classes to access the set precursor mass tolerance.
	double getPrecursorMassTolerance() const {
		return precursorMassTolerance_;
	}

	/// Allows derived classes to use the utilities in the base class.
	const Utilities * getUtils() const {
		return(&utils_);
	}

};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_MATER_H
