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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_RANDOMMUTATER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_RANDOMMUTATER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mutater.h>
#include <map>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI RandomMutater : Mutater
  {
private:
	/// The vector holds the weights for the random decision of which Mutater to use.
	/// The weights are increasing with the size of the vector and the last double must be 1.0.
    std::vector<double> weights;

public:
    /// identifier for SubstitutingMutater
    const static int substitutingMutater = 0;
    /// identifier for SwappingMutater
    const static int swappingMutater = 1;
    /// identifier for InvertingMutater
    const static int invertingMutater = 2;

    RandomMutater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

    /// Generates a random number between 0 and 1 and then chooses one of the defined Mutaters to
    /// mutate the Chromosome. Returns the mutated Chromosome.
    virtual void mutate(Chromosome& chromosome);

    /// Returns the weights currently set for the Mutaters.
	const std::vector<double>& getWeights() const
	{
		std::vector<double> ret;
		for(int i = weights.size()-2; i >= 0; i--)
		{
	      double val = weights[i+1] - weights[i];
	      std::vector<double>::iterator it = ret.begin();
		  ret.insert(it, val);
		}
		return ret;
	}

	/// Sets the input weights for the decision which Mutater to use
	/// Only accepts as many weights as exist Mutater implementations and forces the last element to be 1.
	/// Weights must be given such that they sum up to 1 e.g.: {0.3,0.4,0.3}.
	void setWeights(const std::vector<double>& weights) {
		this->weights[0] = weights[0];
		for(int i=1; i<this->weights.size(); i++)
		{
		  this->weights[i] = weights[i]+this->weights[i-1];
		}
		if(this->weights[this->weights.size()-1] < 1 || this->weights[this->weights.size()-1] > 1)
			this->weights[this->weights.size()-1] = 1.0;
	}
};
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_RANDOMMUTATER_H
