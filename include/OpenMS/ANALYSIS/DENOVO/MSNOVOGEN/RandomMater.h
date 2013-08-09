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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_RANDOMMATER_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_RANDOMMATER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleMater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultMater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/ZipMater.h>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  class OPENMS_DLLAPI RandomMater : public Mater
  {
public:
  /// identifier for SubstitutingMutater
  const static int simpleMater = 0;
  /// identifier for SwappingMutater
  const static int defaultMater = 1;
  /// identifier for InvertingMutater
  const static int zipMater = 2;

private:
	/// The vector holds the weights for the random decision of which Mutater to use.
	/// The weights are increasing with the size of the vector and the last double must be 1.0.
    std::vector<double> weights_;
	SimpleMater sm;
	DefaultMater dm;
	ZipMater zm;

private:
	/// Copy c'tor
	RandomMater(const RandomMater& other);

	/// Assignment operator
	RandomMater & operator=(const RandomMater& rhs);
public:
    /// Default c'tor
    RandomMater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);
    std::vector<boost::shared_ptr<Chromosome> > mate(const boost::shared_ptr<Chromosome> lhs, const boost::shared_ptr<Chromosome> rhs) const;
    void seed(unsigned int seed);
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
		this->weights_[0] = weights[0];
		for(unsigned int i=1; i<this->weights_.size(); i++)
		{
		  this->weights_[i] = weights[i]+this->weights_[i-1];
		}
		if(this->weights_[this->weights_.size()-1] < 1 || this->weights_[this->weights_.size()-1] > 1)
			this->weights_[this->weights_.size()-1] = 1.0;
	}
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_RANDOMMATER_H
