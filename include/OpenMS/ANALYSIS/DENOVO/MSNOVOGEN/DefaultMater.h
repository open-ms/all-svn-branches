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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_DEFAULTMATER_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_DEFAULTMATER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleMater.h>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
	/**
	* @brief The DefaultMater calls one of the classes derived from Mater
	* and uses it to perform crossover of the two provided individuals.
	* Currently the DefaultMater is a SimpleMater.
	*/
  class OPENMS_DLLAPI DefaultMater : public Mater
  {
private:
	/// The implementation of Mater that is the current default (SimpleMater).
	SimpleMater dm_;

	/// Copy c'tor shouldn't be used.
	DefaultMater(const DefaultMater& other);

	/// Assignment operator shouldn't be used.
	DefaultMater & operator=(const DefaultMater& rhs);

public:
    /// Default c'tor providing all necessary input.
    DefaultMater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

    /// The implementation of mate from Mater.
    std::vector<boost::shared_ptr<Chromosome> > mate(boost::shared_ptr<Chromosome> lhs, const boost::shared_ptr<Chromosome> rhs) const;

    /// Overriding seed since it must be forwarded to the contained mating implementation (dm_).
    void seed(const unsigned int seed);

    /// Overriding getSeed since it must be requested from the mating implementation (dm_).
    unsigned int getSeed();
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_DEFAULTMATER_H
