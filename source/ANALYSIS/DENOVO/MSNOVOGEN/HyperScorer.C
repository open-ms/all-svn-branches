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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/HyperScorer.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>

namespace OpenMS
{

  HyperScorer::HyperScorer(const double fmt) :
    Scorer(fmt)
  {}

  void HyperScorer::score(const MSSpectrum<> * expMS, boost::shared_ptr<Chromosome> & chromosome) const
  {
	
	TheoreticalSpectrumGenerator generator;
	Param p(generator.getParameters());
	p.setValue("add_metainfo", "true");
	generator.setParameters(p);
	RichPeakSpectrum sp;
	// generate spectrum for b- and y-ions single and doubly charged
	for (int c = 1; c<3; c++) {
		generator.addPeaks(sp, chromosome->getSequence(), Residue::BIon, c);
		generator.addPeaks(sp, chromosome->getSequence(), Residue::YIon, c);
	}

	std::vector<Peak1D>::const_iterator p_it = expMS->begin();
	std::vector<RichPeak1D>::const_iterator tp_it = sp.begin();
	double score = 0.0;
	Size yN = 0;
	Size bN = 0;
	while (p_it != expMS->end() && tp_it != sp.end()) {
		// check if peaks match in m/z
		const double p_mz = p_it->getMZ();
		const double tp_mz = tp_it->getMZ();
		double dist = abs(p_mz - tp_mz);
		double dist_old;
		if (abs(p_mz - tp_mz) < tp_mz * getFragmentMassTolerance() * 1e-6) {
			// find the closest peak
			double s_tmp = 0.0;
			do {
				dist_old = dist;
				s_tmp = p_it->getIntensity();
				++p_it;
				if(p_it == expMS->end()) {
					--p_it;
					break;
				}
				dist = abs(p_it->getMZ() - tp_mz);
			} while (dist < dist_old && p_it != expMS->end());
			score += s_tmp;
			if (tp_it->getMetaValue("IonName").toString()[0]
					== 'y')
				++yN;
			else
				++bN;
			++tp_it;
			continue;
		}
		// peaks do not match, increase respective peak
		if (p_mz < tp_mz) {
			++p_it;
		} else {
			++tp_it;
		}
	}
	// found hit...
	if (score > 0) {
		score *= fac(bN) * fac(yN);
	}
	chromosome->setScore(score);
  }

	Size HyperScorer::fac(const Size i) const
	{
		if (i < 2)
			return 1;
		return i * fac(i - 1);
	}

} // namespace
