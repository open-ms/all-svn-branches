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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Killer.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Scorer.h>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /**
  * @brief The GenPool is the collection of Chromosomes (individuals) and it
  * also allows to perform operations on the pool like mutations and crossovers.
  *
  */
  class OPENMS_DLLAPI GenPool
  {
public:
	  enum {BYSCOREDEC = 0};

private:
	  unsigned int randomSeed_;
	  std::vector<boost::shared_ptr<Chromosome> > genPool;
	  std::map<String, boost::shared_ptr<Chromosome> > knownIndividuals;
	  const Mutater* mutater;
	  const Mater* mater;
	  const Killer* killer;
	  const Seeder* seeder;
	  const Scorer* scorer;
	  unsigned int maxPoolSize;
	  unsigned int previousPoolSize;
	  double precursorMass;
	  int precursorCharge;
	  double precursorIntensity;
	  double precursorMZ;
	  double precursorMassTolerance;
	  std::vector<const Residue *> aaList;
	  const MSSpectrum<> * msms;

	  ///not implemented, should not be copied
	  GenPool( const GenPool& );
	  ///not implemented, should not be assigned.
	  GenPool & operator=(const GenPool &);

	  /// Takes the most abundant of availbable precursors and assigns m/z, intensity and charge to this.
	  void setPrecursorInfo(const MSSpectrum<> * ms);

public:
    /// Default c'tor
    GenPool(const int maxPoolSize, const double precursorMassTolerance);

    ~GenPool()
    {};

    void mutate();
    void setMutater(const Mutater* mutater);

    void crossover();
    void setMater(const Mater * mater);

    void kill();
    void setKiller(const Killer * killer);

    /// Tries to initialize the gen pool with maxPoolSize individials.
    /// If not successful within 10*maxPoolSize tries, gives up.
    void initGenPool(const int maxPoolSize);
    void setSeeder(const Seeder * seeder);

    void setScorer(const Scorer * scoer);
    void score();

    Size getPopulationSize();
    void sort(const int sortMethod = Chromosome::sortScoreDescending);
    void setPool(std::vector<boost::shared_ptr<Chromosome> > newPool);

    /// Add a Chromosome which must not have been seen before.
    /// Returns true if added and false if not.
    bool addIndividual(boost::shared_ptr<Chromosome> individual);

    int addIndividuals(std::vector<boost::shared_ptr<Chromosome> > individuals);

    const boost::shared_ptr<Chromosome> getIndividual(const Size which);

    void replenish(const int targetSize);


    std::vector<boost::shared_ptr<Chromosome> >::iterator begin()
    {
      return(genPool.begin());
    };

    std::vector<boost::shared_ptr<Chromosome> >::iterator end() {
    	return(genPool.end());
    };
	void seed(const unsigned int seed)
	{
		randomSeed_ = seed;
	}

	unsigned int getSeed() const
	{
		return randomSeed_;
	}

	const std::vector<const Residue*>& getAaList() const {
		return aaList;
	}

	void setAaList(const std::vector<const Residue*>& aaList) {
		this->aaList = aaList;
	}

	const std::vector<boost::shared_ptr<Chromosome> >& getGenPool() const {
		return genPool;
	}

	void setGenPool(std::vector<boost::shared_ptr<Chromosome> > genPool);

	const std::map<String, boost::shared_ptr<Chromosome> > getKnownIndividuals() const {
		return knownIndividuals;
	}

	void setKnownIndividuals(
			const std::map<String, boost::shared_ptr<Chromosome> >& knownIndividuals) {
		this->knownIndividuals = knownIndividuals;
	}

	int getMaxPoolSize() const {
		return maxPoolSize;
	}

	void setMaxPoolSize(int maxPoolSize) {
		this->maxPoolSize = maxPoolSize;
	}

	double getPrecursorMassTolerance() const {
		return precursorMassTolerance;
	}

	void setPrecursorMassTolerance(double precursorMassTolerance) {
		this->precursorMassTolerance = precursorMassTolerance;
	}

	unsigned int getRandomSeed() const {
		return randomSeed_;
	}

	void setRandomSeed(unsigned int randomSeed) {
		randomSeed_ = randomSeed;
	}

	void setMSMSSpectrum(const MSSpectrum<> * msms);

	unsigned int getPreviousPoolSize() const {
		return previousPoolSize;
	}

	void setPreviousPoolSize(unsigned int previousPoolSize) {
		this->previousPoolSize = previousPoolSize;
	}

	int getPrecursorCharge() const {
		return precursorCharge;
	}

	double getPrecursorIntensity() const {
		return precursorIntensity;
	}

	double getPrecursorMass() const {
		return precursorMass;
	}

	double getPrecursorMz() const {
		return precursorMZ;
	}
};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_GENPOOL_H
