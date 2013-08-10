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

#ifndef OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SWAPPINGMUTATER_H
#define OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SWAPPINGMUTATER_H

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mutater.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
	/**
	* @brief The SwappingMutater selects an amino acid
	* from the first part of the sequence and another
	* one from the second part of the sequence.
	* It then swaps the two amino acids.
	*/
  class OPENMS_DLLAPI SwappingMutater : public Mutater
  {
	  
private:
    /// To forbid copy construction
    SwappingMutater(const SwappingMutater & rhs);
    /// To forbid assignment
    SwappingMutater & operator=(const SwappingMutater & rhs);

public:
	/// constructor providing all neccessary information.
	SwappingMutater(double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList);
	
	///d'tor
	virtual ~SwappingMutater();

	/// Implementation of virtual method Mutater::mutate.
	void mutate(boost::shared_ptr<Chromosome> chromosome) const;
  };
} // namespace

#endif // OPENMS_ANALYSIS_DENOVO_MSNOVOGEN_SWAPPINGMUTATER_H
