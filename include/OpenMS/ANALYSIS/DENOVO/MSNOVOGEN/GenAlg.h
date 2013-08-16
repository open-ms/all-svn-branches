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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_GENALG_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_GENALG_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Killer.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Scorer.h>
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/HomologyKiller.h"
#include <boost/shared_ptr.hpp>
#include <vector>

namespace OpenMS
{
  class GenPool;

  class OPENMS_DLLAPI GenAlg
  {

private:
	/// The spectrum on which this genetic algorithm will be applied.
	const MSSpectrum<> * msms_;
	/// A list of expected amino acids.
	const std::vector<const Residue *> aaList_;
	/// After this amount of generations the HomologyKiller will be invoked to make the pool more dissimilar.
	Size rejuvenateAfter_;
	/// The mass tolerance for precursor ions or sequences to be acceptable for the pool
	double precursorMassTolerance_;
	/// The fragment ion tolerance for which peaks may be used during scoring.
	double fragmentMassTolerance_;
	/// The uncharged precursor mass.
	double precursorMass_;
	/// Number of charges of the precursor.
	Size precursorCharge_;
	/// The measured intensity of the precursor.
	double precursorIntensity_;
	/// The unprocessed precursor m/z.
	double precursorMZ_;
	/// random number generator must be mutable since Mater must be const * in GenAlg.
	mutable boost::mt19937 rng_;
	/// The pool holding the Chromosomes.
	GenPool* genPool_;
	/// The Mutater to be used (default: DefaultMutater).
	boost::shared_ptr<const Mutater> mutater_;
	/// The Mater to be used (default: DefaultMater).
	boost::shared_ptr<const Mater> mater_;
	/// The Killer to be used (default: DefaultKiller).
	boost::shared_ptr<const Killer> killer_;
	/// The Seeder to be used (default: DefaultSeeder).
	boost::shared_ptr<const Seeder> seeder_;
	/// The Scorer to be used (default: NormShrAbuScorer).
	boost::shared_ptr<const Scorer> scorer_;
	/// HomologyKiller kills all individuals that are too similar sequence wise (default: edit distance 1).
	boost::shared_ptr<const HomologyKiller> homologyKiller_;

	/// Copy c'tor
	GenAlg(const GenAlg& other);

	/// Assignment operator
	GenAlg & operator=(const GenAlg& rhs);

	/// Checks whether the functions such as Mutater, Mater, Killer, Seeder, Scorer and HomologyKiller are available.
	/// Returns true if everything is fine and false otherwise.
	bool checkEnvironment();

	/// Initializes the GenPool with individuals using the set Seeder up to the number of max pool size.
	/// Uses a maximum of 20 * maxPoolSize tries to do so.
	void initGenPool();

	/// Calls mutate of the set Mutater on all elements of the GenPool.
    void mutate();

    /// Calls tournament of the set Mater for this GenPool.
    void crossover();

    /// Calls kill of the set Killer for this GenPool.
    void kill();

    /// Calls score for the given Scorer for all elements of the given GenPool if they have not been scored.
    void score();

    /// sets the given MSSpectrum and derives the precursor information for the most abundant precursor ion.
    void setMSMSSpectrum(const MSSpectrum<> * msms);

	/// Takes the most abundant of availbable precursors and assigns m/z, intensity and charge to this.
	void setPrecursorInfo(const MSSpectrum<> * ms);

public:
    /// Default c'tor
    GenAlg(const MSSpectrum<>* spec, const std::vector<const Residue *>& aaList, const Size maxPoolSize, const double precursorMassTolerance, const double fragmentMassTolerance);

    /// Default destructor only deletes genPool.
    ~GenAlg();

    /// Performs this genetic algorithm for numGenerations generations.
    /// Calls startEvolution(numGenerations, numGenerations) to do so.
    PeptideIdentification startEvolution(const Size numGenerations, const Size bestNHits = 10);

    /// Performs this genetic algorithm for numGenerations generations.
    /// or until the same sequence is the best sequence for the set number of generations (endIfStableForNumGenerations).
    /// calls:
    ///         crossover();
    ///			mutate();
    ///			score();
    ///			if((i % rejuvenateAfter_) == 0)
    ///			{
    ///				homologyKiller_->kill(*genPool_);
    ///				genPool_->replenish(genPool_->getPreviousPoolSize());
    ///			}
    ///			kill();
    /// for each generation.
    PeptideIdentification startEvolution(const Size numGenerations, const Size endIfStableForNumGenerations, const Size bestNHits = 10);

    /// Allows the setting of the custom Mutater derived from Mutater.
    void setMutater(boost::shared_ptr<const Mutater> mutater);

    /// Allows the setting of the custom Mater derived from Mater.
    void setMater(boost::shared_ptr<const Mater> mater);

    /// Allows the setting of the custom Killer derived from Killer.
    void setKiller(boost::shared_ptr<const Killer>  killer);

    /// Allows the setting of the custom Seeder derived from Seeder.
    void setSeeder(boost::shared_ptr<const Seeder>  seeder);

    /// Allows the setting of the custom Scorer derived from Scorer.
    void setScorer(boost::shared_ptr<const Scorer>  scorer);

    /// Allows to set a customized HomologyKiller i.e.: with larger edit distance.
    void setHomologyKiller(boost::shared_ptr<const HomologyKiller>  killer);

    /// Allows to set the random seed for testing purposes or for adding another level of randomness during the algorithm.
	void seed(const Size seed) const;

	/// Returns the amino acids that are known to this object.
	const std::vector<const Residue*>& getAAList() const {
		return aaList_;
	}

	/// Returns all different sequences that were generated during the evolution.
	const std::map<String, boost::shared_ptr<Chromosome> > getKnownIndividuals() const;

	/// Allows the possibility to seed the pool with some known possible Individuals.
	/// All individuals will first be scored using the set Scorer.
	void setKnownIndividuals(const std::map<String, boost::shared_ptr<Chromosome> >& knownIndividuals);

	/// returns the precursor intensity used for this object.
	double getPrecursorIntensity() const {
		return precursorIntensity_;
	}

	/// returns the calculated precursor MH+ for this object.
	double getPrecursorMH() const {
		return(precursorMass_ + 1);
	}

	/// returns the uncharged mass of the precursor used for this object.
	double getPrecursorMass() const {
		return precursorMass_;
	}

	/// returns the measured m/z of the precursor used for this object.
	double getPrecursorMZ() const {
		return precursorMZ_;
	}

	/// Sets after how many generations the HomologyKiller should be called such that the population doesn't become to similar.
	void setRejuvenate(const Size afterNGen)
	{
		rejuvenateAfter_ = afterNGen;
	}

	/// Retrieve the individuals in the curreny genPool
	const GenPool * getGenPool() const {
		return const_cast<const GenPool *>(genPool_);
	}
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_GENALG_H
