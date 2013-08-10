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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_ZIPMATER_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_ZIPMATER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
	/**
	* @brief The ZipMater takes two Chromosomes and swaps every second amino acid.
	* Leaves overhangs unmodified.
	* The mass of the resulting Chromosomes may not confirm to the
	* required precursor MH+ +/- tolerance so Utilities::adjustToFitMass is called
	* on both children in order to try to change the sequence such that it suits the mass constraint.
	* May return an empty vector.
	*/
  class OPENMS_DLLAPI ZipMater : public Mater
  {
private:
	/// Copy c'tor shouln't be used.
	ZipMater(const ZipMater& other);

	/// Assignment operator shouldn't be sued.
	ZipMater & operator=(const ZipMater& rhs);

public:
    /// Default c'tor providing all necessary information.
    ZipMater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);

	/// Implementation of the virtual method Mater::mate.
    std::vector<boost::shared_ptr<Chromosome> > mate(const boost::shared_ptr<Chromosome> lhs, const boost::shared_ptr<Chromosome> rhs) const;

  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_ZIPMATER_H
