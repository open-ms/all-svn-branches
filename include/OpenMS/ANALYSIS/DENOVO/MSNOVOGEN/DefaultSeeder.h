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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_DEFAULTSEEDER_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_DEFAULTSEEDER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Seeder.h>
#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/RandomSequenceSeeder.h>

namespace OpenMS
{
	/**
	* @brief The DefaultSeeder calls one of the classes derived from Seeder
	* and uses it to perform seeding of an individual.
	* Currently the DefaultSeeder is a RandomSequenceSeeder.
	*/
  class OPENMS_DLLAPI DefaultSeeder : public Seeder
  {
private:
	/// The currently set DefaultSeeder (RandomSequenceSeeder).
	RandomSequenceSeeder ds_;

	/// To forbid copy construction
	DefaultSeeder(const DefaultSeeder& other);

	/// To forbid assignment
	DefaultSeeder & operator=(const DefaultSeeder& rhs);

public:
	/// Default c'tor providing all necessary input.
	DefaultSeeder(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

	/// Implementation of the virtual method Seeder::createIndividual.
	boost::shared_ptr<Chromosome> createIndividual() const;

	/// Overridden to forward the seed to the default seeder.
	void seed(const unsigned int seed);

	/// Overridden to request the seed from the default seeder.
	unsigned int getSeed() const
	{
		return ds_.getSeed();
	}
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_DEFAULTSEEDER_H
