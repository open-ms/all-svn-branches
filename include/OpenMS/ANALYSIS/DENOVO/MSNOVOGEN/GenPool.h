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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_GENPOOL_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_GENPOOL_H

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <boost/random/mersenne_twister.hpp>

namespace OpenMS
{
  /**
  * @brief The GenPool is the collection of Chromosomes (individuals).
  * It further stores all individuals that ever lived in the pool as knownIndividuals.
  *
  */
  class OPENMS_DLLAPI GenPool
  {
private:
	  /// random number generator must be mutable since Mater must be const * in GenAlg.
	  mutable boost::mt19937 rng;
	  /// The genPool during one generation (can grow beyond maxPoolSize during mating).
	  std::vector<boost::shared_ptr<Chromosome> > genPool_;
	  /// The poos of all different individuals that were seen in this pool in any generation.
	  std::map<String, boost::shared_ptr<Chromosome> > knownIndividuals_;
	  /// The maximum pool size which will be enforced during killing at the end of any generation.
	  unsigned int maxPoolSize_;
	  /// The previous pool size so that if a decreasing killer is used the new generation size can be determined.
	  unsigned int previousPoolSize_;
	  /// The precursor MH+.
	  double precursorMass_;
	  /// The charge of the precursor.
	  int precursorCharge_;
	  /// The precursor mass tolerance to be used to accept individuals into the pool.
	  double precursorMassTolerance_;

	  ///not implemented, should not be copied
	  GenPool( const GenPool& );
	  ///not implemented, should not be assigned.
	  GenPool & operator=(const GenPool &);


public:
    /// Default c'tor
    GenPool(const int maxPoolSize, const double precursorMassTolerance, const double precursorMass, const int precursorCharge);

    /// Default d'tor.
    ~GenPool() {};

    /// Returns the size of the current population.
    Size getPopulationSize() const;

    /// Sorts the current gen pool by score either decreasing (default) or increasing (Chromosome::sortScoreAscending).
    void sort(const int sortMethod = Chromosome::sortScoreDescending);

    /// Add a Chromosome which must not have been seen before (i.e.: is not in known individuals).
    /// In the current gene pool duplicates can occur during replenish().
    /// Returns true if added and false if not.
    bool addIndividual(boost::shared_ptr<Chromosome> individual);

    /// Calls addIndividual for each individual in the supplied vector.
    int addIndividuals(std::vector<boost::shared_ptr<Chromosome> > individuals);

    /// Returns the specified individual if which is within range.
    const boost::shared_ptr<Chromosome> getIndividual(const Size which);

    /// Draws random individuals from the knownIndividuals and adds them to the current pool
    /// to reach the targetSize which must be <= previousPoolSize (otherwise it will be adjusted automatically).
    void replenish(const int targetSize);

    /// returns the begin of the current gen pool.
    std::vector<boost::shared_ptr<Chromosome> >::iterator begin()
    {
      return(genPool_.begin());
    };

    /// returns the end of the current gen pool.
    std::vector<boost::shared_ptr<Chromosome> >::iterator end() {
    	return(genPool_.end());
    };

    /// Allows to set the random seed for testing purposes or for adding another level of randomness during the algorithm.
	void seed(const Size seed)
	{
		rng.seed(seed);
	}

	/// Returns the currently set random number generator to be used in child classes.
	boost::mt19937 getRandomGenerator() const {
		return(rng);
	}

	/// returns the vector containing the current gen pool.
	const std::vector<boost::shared_ptr<Chromosome> >& getGenPool() const {
		return genPool_;
	}

	/// Allows to set the current gen pool.
	void setGenPool(std::vector<boost::shared_ptr<Chromosome> > genPool);

	/// Allows to add multiple individuals at a time.
	void addToGenPool(std::vector<boost::shared_ptr<Chromosome> > genPool);

	/// returns the map with all known individuals.
	const std::map<String, boost::shared_ptr<Chromosome> > getKnownIndividuals() const {
		return knownIndividuals_;
	}

	/// Allows the possibility to seed the pool with some known possible Individuals.
	void setKnownIndividuals(
			const std::map<String, boost::shared_ptr<Chromosome> >& knownIndividuals) {
		this->knownIndividuals_ = knownIndividuals;
	}

	/// Returns the previous pool size since during mating it can increase arbitrarily.
	unsigned int getPreviousPoolSize() const {
		return previousPoolSize_;
	}

	/// allows to set the previous pool size for decreasing Killers for instance.
	void setPreviousPoolSize(unsigned int previousPoolSize) {
		this->previousPoolSize_ = previousPoolSize;
	}

	/// returns the maximum pool size that was set during construction.
	unsigned int getMaxPoolSize() const {
		return maxPoolSize_;
	}

	/// Clear the current gene pool.
	void clear() {
		genPool_.clear();
	}
};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_GENPOOL_H
