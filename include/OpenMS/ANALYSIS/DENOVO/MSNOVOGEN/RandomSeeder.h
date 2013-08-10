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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_RANDOMSEEDER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_RANDOMSEEDER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomSequenceSeeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SequenceTagSeeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultSeeder.h>
#include <vector>

namespace OpenMS
{
	/**
	* @brief The RandomSeeder calls one of the classes derived from Seeder
	* and uses it to perform seeding of news individuals.
	* Which implementation is randomly chosen but can be influenced by
	* setting the weights. The number of weights must equal the number
	* of available maters (currently 3). 
	*/
  class OPENMS_DLLAPI RandomSeeder : public Seeder
  {
private:
	/// The vector holds the weights for the random decision of which Mutater to use.
	/// The weights are increasing with the size of the vector and the last double must be 1.0.
	std::vector<double> weights_;
	/// The RandomSequenceSeeder is a contained object so it doesn't need to be instantiated more than once.
	RandomSequenceSeeder rss;
	/// The SequenceTagSeeder is a contained object so it doesn't need to be instantiated more than once.
	SequenceTagSeeder sts;
	/// The DefaultSeeder is a contained object so it doesn't need to be instantiated more than once.
	DefaultSeeder ds;
	
private:
	/// To forbid copy construction
	  RandomSeeder(const RandomSeeder& other);

	/// To forbid assignment
	  RandomSeeder & operator=(const RandomSeeder& rhs);

public:
	/// identifier for SubstitutingMutater
	const static int randomSequenceSeeder = 0;
	/// identifier for SwappingMutater
	const static int sequenceTagSeeder = 1;
	/// identifier for DefaultSeeder
	const static int defaultSeeder = 2;

    /// Default c'tor
    RandomSeeder(const MSSpectrum<> * spec, const double precursorMass, const double precursorMassTolerance, const double fagmentMassTolerance, const std::vector<const Residue*> aaList);

	/// Implementation of the virtual method Seeder::createIndividual.
    boost::shared_ptr<Chromosome> createIndividual() const;

	/// Ovrridden to forward the seed to the contained objects.
	void seed(const unsigned int seed);

    /// Returns the weights currently set for the Mutaters.
	const std::vector<double> getWeights() const
	{
		if(weights_.size() < 1)
			throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		std::vector<double> ret;
		double val = weights_[0];
		ret.push_back(val);
		for(Size i = 1; i < weights_.size(); i++)
		{
	      val = weights_[i] - weights_[i-1];
		  ret.push_back(val);
		}
		return ret;
	}

	/// Sets the input weights for the decision which Mutater to use
	/// Only accepts as many weights as exist Mutater implementations and forces the last element to be 1.
	/// Weights must be given such that they sum up to 1 e.g.: {0.3,0.4,0.3}.
	void setWeights(const std::vector<double>& weights) {
		this->weights_.clear();
		this->weights_.resize(weights.size());
		this->weights_[0] = weights[0];
		for(unsigned int i=1; i<weights.size(); i++)
		{
		  this->weights_[i] = weights[i]+this->weights_[i-1];
		}
		if(this->weights_[weights.size()-1] < 1 || this->weights_[weights.size()-1] > 1)
			this->weights_[weights.size()-1] = 1.0;
	}
  };
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_RANDOMSEEDER_H
