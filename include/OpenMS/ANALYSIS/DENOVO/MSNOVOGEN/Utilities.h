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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_UTILITIES_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_UTILITIES_H

#include <OpenMS/config.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI Utilities
  {
private:
	/// Copy c'tor
	Utilities(const Utilities& other);

	/// Assignment operator
	Utilities & operator=(const Utilities& rhs);

  public:
    /// Default c'tor
    Utilities();

	/// Generates a random sequence of the desired length and then adjusts the
	/// amino acids until the weight matches the desired criteria.
	static const AASequence getRandomSequence(const int seed, const int len, const double weight, const double tolerance, const std::vector<const Residue *> aaList);

    /// from the list of available amino acids selects a random one.
    static const String getRandomAA(const int seed, const std::vector<const Residue *> aaList);

    /// The random string may not fit to the expected precursor mass and is then adjusted to fit (if possible with a maximum of 3 changes).
    /// The passed in sequence is directly modified.
    static bool adjustToFitMass(const int seed, AASequence & sequence, const double weight, const double tolerance, const std::vector<const Residue *> aaList);

    static int editDistance(const AASequence & lhs, const AASequence & rhs);
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_UTILITIES_H
