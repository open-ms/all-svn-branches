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

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_RESIDUESERVER_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_RESIDUESERVER_H

#include <OpenMS/config.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <vector>
#include <map>

namespace OpenMS
{
	/**
	* @brief The ResidueServer provides Residues fast.
	*/
  class OPENMS_DLLAPI ResidueServer
  {

private: 
	/**
	*/
	std::multimap<const Size, const Residue*> massResidueMap_;
	std::multimap<const String, const Residue*> nameResidueMap_;
	std::set<const Residue*> residues_;

	/// Copy c'tor shouldn't be used.
	ResidueServer(const ResidueServer& other);

	/// Assignment operator shouldn't be used.
	ResidueServer & operator=(const ResidueServer& rhs);

public:
    /// Default c'tor providing all necessary input.
    ResidueServer();
	~ResidueServer(){};
	
	std::vector<const Residue *> getResidues(const double mass, const double tolerance) const;

	void add(const Residue * residue)
	{
		residues_.insert(residue);
		massResidueMap_.insert(std::pair<const Size, const Residue*>((Size)residue->getMonoWeight(),residue));
		nameResidueMap_.insert(std::pair<const String, const Residue*>(residue->getOneLetterCode(),residue));
	};

	void add(const std::vector<const Residue *> & residues)
	{
		for(Size i=0; i<residues.size(); i++)
			add(residues[i]);
	};
    
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_RESIDUESERVER_H
