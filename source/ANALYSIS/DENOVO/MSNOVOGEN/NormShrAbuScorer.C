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

#include <OpenMS//ANALYSIS/DENOVO/MSNOVOGEN/NormShrAbuScorer.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>

namespace OpenMS
{

  NormShrAbuScorer::NormShrAbuScorer(const double fmt) :
    Scorer(fmt)
  {}

  std::vector<Peak1D> NormShrAbuScorer::getPeaksInRange(
		const MSSpectrum<> & msms,
		const double rangeStart,
		const double rangeEnd,
		std::vector<Peak1D>::const_iterator start,
		std::vector<Peak1D>::const_iterator end        ) const
  {
	std::vector<Peak1D> ret;

    return ret;
  }

  void NormShrAbuScorer::score(const MSSpectrum<> & msms, boost::shared_ptr<Chromosome> chromosome) const
  {
	  double score = 0;
	  double ems = 0; //abundance of peaks in experimental spectrum not shared with theoretical spectrum.
	  int epue = 0;   //number of peaks in experimental spectrum not in theoretical spectrum.
	  int tpms = 0;   //number of peaks in theoretical spectrum not in experimental specetrum.
	  double sms = 0; //abundance of shared peaks (taken from theoretical spectrum).
	  double rangeStart;
	  double rangeEnd;
	  AASequence peptide = chromosome->getSequence();
	  TheoreticalSpectrumGenerator tsg;
	  RichPeakSpectrum rms;
	  tsg.getSpectrum(rms,peptide,2);
	  std::vector<Peak1D>::const_iterator curExp = msms.begin();
	  std::vector<Peak1D>::const_iterator end = msms.end();
	  for(std::vector<RichPeak1D>::const_iterator iter = rms.begin(); iter != rms.end(); iter++)
	  {
		  rangeStart = iter->getMZ() - getFragmentMassTolerance();
		  rangeEnd = iter->getMZ() + getFragmentMassTolerance();
		  std::vector<Peak1D> sp = getPeaksInRange(msms,rangeStart,rangeEnd,curExp, end);
	  }
	  chromosome->setScore(score);
  }
} // namespace
