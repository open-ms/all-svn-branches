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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SequenceTagSeeder.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Utilities.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <time.h>

namespace OpenMS
{
	SequenceTagSeeder::SequenceTagSeeder(const MSSpectrum<> * ms, const double pm, const double pmt, const double fmt, const std::vector<const Residue*> al) :
	  Seeder(ms, pm, pmt,fmt, al)
	{
		seed((unsigned int)time(0));
	}

	boost::shared_ptr<Chromosome> SequenceTagSeeder::createIndividual() const
	{
		std::vector<OpenMS::SequenceTagSeeder::SeqTag> tags = createSequenceTags();
		return createIndividual(tags);
	}

	boost::shared_ptr<Chromosome> SequenceTagSeeder::createIndividual(std::vector<OpenMS::SequenceTagSeeder::SeqTag> tags) const
	{
		Size rs = rand() % tags.size();
		OpenMS::SequenceTagSeeder::SeqTag rst = tags[rs];
		int len = rst.before_ / 110;
		AASequence beg = Utilities::getRandomSequence(getSeed(), len, rst.before_+19 /* method assumes full sequence */, getPrecursorMassTolerance(), getAAList());
		len = rst.after_ / 110;
		AASequence end = Utilities::getRandomSequence(getSeed(), len, rst.after_+19 /* method assumes full sequence */, getPrecursorMassTolerance(), getAAList());
		String fullSeq = beg.toString() + rst.seq_ + end.toString();
		AASequence seq(fullSeq);
		if(Utilities::adjustToFitMass(getSeed(),seq,getPrecursorMass(),getPrecursorMassTolerance(),getAAList()))
			return boost::shared_ptr<Chromosome>(new Chromosome(seq,1));
		return boost::shared_ptr<Chromosome>(new Chromosome());
	}

	std::vector<boost::shared_ptr<Chromosome> > SequenceTagSeeder::createIndividuals(const Size num) const
	{
		std::vector<OpenMS::SequenceTagSeeder::SeqTag> tags = createSequenceTags();
		std::vector<boost::shared_ptr<Chromosome> > ret;
		Size maxTries = 20 * num;
		Size curIter = 0;
		while(curIter++ < maxTries) 
		{
			boost::shared_ptr<Chromosome> ni = createIndividual(tags);
			if(ni)
				ret.push_back(ni);
			if(ret.size() >= num)
				break;
		}
		return ret;
	}

	std::vector<Peak1D> SequenceTagSeeder::getMatchingPeaks(const MSSpectrum<> * msms, const Peak1D * p1, const Residue * aa1) const
	{
		std::vector<Peak1D> ret;
		double rangeStart = p1->getMZ() + aa1->getMonoWeight(Residue::Internal) - getFragmentMassTolerance();
		double rangeEnd = rangeStart + 2 * getFragmentMassTolerance();
		MSSpectrum<>::ConstIterator wndBeg = msms->MZBegin(msms->begin(), rangeStart, msms->end());
		MSSpectrum<>::ConstIterator wndEnd = msms->MZEnd(wndBeg, rangeEnd, msms->end());
		for(MSSpectrum<>::ConstIterator ri=wndBeg; ri != wndEnd; ri++)
		{
			ret.push_back(*ri);
		}
		return ret;
	}

	std::vector<OpenMS::SequenceTagSeeder::SeqTag> SequenceTagSeeder::createSequenceTags() const
	{
		std::vector<OpenMS::SequenceTagSeeder::SeqTag> ret;
		MSSpectrum<> * lMS = const_cast<MSSpectrum<> *>(getMassSpectrum());
		lMS->sortByPosition();
		const std::vector<const Residue*> laa = getAAList();
		double pm = getPrecursorMass();
		for(std::vector<Peak1D>::const_iterator p1 = lMS->begin(); p1 != lMS->end(); p1++)
		{
			for(std::vector<const Residue*>::const_iterator aa1 = laa.begin(); aa1 != laa.end(); aa1++)
			{
				std::vector<Peak1D> mps1 = getMatchingPeaks(lMS,&(*p1),*aa1);
				if(mps1.empty()) continue;
				for(std::vector<Peak1D>::iterator p2 = mps1.begin(); p2 != mps1.end(); p2++)
				{
					for(std::vector<const Residue*>::const_iterator aa2 = laa.begin(); aa2 != laa.end(); aa2++)
					{
						std::vector<Peak1D> mps2 = getMatchingPeaks(lMS,&(*p2),*aa2);
						if(mps2.empty()) continue;
						for(std::vector<Peak1D>::iterator p3 = mps2.begin(); p3 != mps2.end(); p3++)
						{
							for(std::vector<const Residue*>::const_iterator aa3 = laa.begin(); aa3 != laa.end(); aa3++)
							{
								std::vector<Peak1D> mps3 = getMatchingPeaks(lMS,&(*p3),*aa3);
								if(mps3.empty()) continue;
								for(std::vector<Peak1D>::iterator p4 = mps3.begin(); p4 != mps3.end(); p4++)
								{
									String seq = (*aa1)->getOneLetterCode() + (*aa2)->getOneLetterCode() + (*aa3)->getOneLetterCode();
									ret.push_back(SeqTag(p1->getMZ(),seq,(pm - p4->getMZ())));
								}
							}
						}
					}
				}
			}
		}
		return ret;
	}
} // namespace
