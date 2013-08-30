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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenAlg.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/NormShrAbuScorer.h>
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultKiller.h"
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultMutater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultMater.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultSeeder.h>
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/HomologyKiller.h"
#include <vector>
#include <time.h>
#include <boost/shared_ptr.hpp>

using std::vector;

namespace OpenMS
{

	GenAlg::GenAlg(const MSSpectrum<>* spec, const std::vector<const Residue *>& al, const Size mps, const double pmt, const double fmt) :
	  aaList_(al),
	  rejuvenateAfter_(mps),
	  precursorMassTolerance_(pmt),
	  fragmentMassTolerance_(fmt),
	  precursorCharge_(1),
	  rng_(time(NULL))
	{
		setMSMSSpectrum(spec);
		genPool_ = new GenPool(mps, precursorMassTolerance_, getPrecursorMH(), precursorCharge_);
		mutater_ = boost::shared_ptr<const Mutater>(new DefaultMutater(getPrecursorMH(),precursorMassTolerance_,aaList_));
		mater_= boost::shared_ptr<const Mater>(new DefaultMater(getPrecursorMH(),precursorMassTolerance_,aaList_));
		killer_= boost::shared_ptr<const Killer>(new DefaultKiller());
		seeder_= boost::shared_ptr<const Seeder>(new DefaultSeeder(msms_,getPrecursorMH(),precursorMassTolerance_,fragmentMassTolerance_,aaList_));
		scorer_= boost::shared_ptr<const Scorer>(new NormShrAbuScorer(getPrecursorMH()));
		homologyKiller_ = boost::shared_ptr<const HomologyKiller>(new HomologyKiller());
	}

	GenAlg::~GenAlg()
	{
		delete genPool_;
	}

	void GenAlg::seed(const Size seed) const
	{
		rng_.seed(seed);
		mutater_->seed(seed);
		mater_->seed(seed);
		seeder_->seed(seed);
	}

	const std::map<String, boost::shared_ptr<Chromosome> > GenAlg::getKnownIndividuals() const
    {
		return genPool_->getKnownIndividuals();
	}

	void GenAlg::setKnownIndividuals(const std::map<String, boost::shared_ptr<Chromosome> >& knownIndividuals)
	{
		std::map<String, boost::shared_ptr<Chromosome> > scored;
		for(std::map<String, boost::shared_ptr<Chromosome> >::const_iterator iter=knownIndividuals.begin(); iter !=knownIndividuals.end(); iter++)
		{
			std::pair<String,boost::shared_ptr<Chromosome> > p = *iter;
			scorer_->score(msms_,p.second);
			scored.insert(p);
		}
		genPool_->setKnownIndividuals(scored);
	}

	void GenAlg::setMutater(boost::shared_ptr<const Mutater> mutater)
	{
		this->mutater_ = mutater;
	}

	void GenAlg::mutate()
	{
	  mutater_->mutateAndAddToPool(*genPool_);
	}

	void GenAlg::setMater(boost::shared_ptr<const Mater> mater)
	{
		this->mater_ = mater;
	}

	void GenAlg::crossover()
	{
		mater_->tournamentAndAddToPool(*genPool_);
	}

	void GenAlg::setKiller(boost::shared_ptr<const Killer> killer)
	{
		this->killer_ = killer;
	}

	void GenAlg::kill()
	{
		killer_->kill(*genPool_);
	}

	void GenAlg::setSeeder(boost::shared_ptr<const Seeder> seeder)
	{
	  this->seeder_ = seeder;
	}

	void GenAlg::setScorer(boost::shared_ptr<const Scorer> scorer)
	{
	this->scorer_ = scorer;
	}

	void GenAlg::score()
	{
		scorer_->scorePool(msms_, *genPool_);
	}

	void GenAlg::setHomologyKiller(boost::shared_ptr<const HomologyKiller>  killer)
	{
		this->homologyKiller_ = killer;
	}

	void GenAlg::setPrecursorInfo(const MSSpectrum<>* ms)
	{
	  double abMax = -1;
	  const vector<Precursor>& precursors = ms->getPrecursors();
	  for (vector<Precursor>::const_iterator pc_it = precursors.begin(); pc_it != precursors.end(); ++pc_it) {
		  if(abMax <= pc_it->getIntensity())
		  {
			  precursorMZ_ = pc_it->getMZ();
			  precursorCharge_ = pc_it->getCharge();
			  precursorMass_ = precursorMZ_ * precursorCharge_ - precursorCharge_;
			  precursorIntensity_ = pc_it->getIntensity();
		  }
	  }
	}

	bool GenAlg::checkEnvironment()
	{
		if(!mutater_)
			return false;
		if(!seeder_)
			return false;
		if(!scorer_)
			return false;
		if(!killer_)
			return false;
		if(!mater_)
			return(false);
		return true;
	}

	PeptideIdentification GenAlg::startEvolution(const Size numGenerations, const Size bestNHits)
	{
		return startEvolution(numGenerations, numGenerations + 1, bestNHits);
	}

	void GenAlg::initGenPool()
	{
		if(genPool_->getPopulationSize() > 0)
			genPool_->clear();
		genPool_->addIndividuals(seeder_->createIndividuals(genPool_->getMaxPoolSize()));

//		Size maxTries = 20 * genPool_->getMaxPoolSize();
//	    while((maxTries > 0) && (genPool_->getPopulationSize() < genPool_->getMaxPoolSize()))
//	    {
//	      boost::shared_ptr<Chromosome> ni = seeder_->createIndividual();
//		  if(ni && ni->getSequence().size() > 0)
//	        genPool_->addIndividual(ni);
//	      maxTries--;
//	    }
	    genPool_->setPreviousPoolSize(genPool_->getPopulationSize());
	}

	PeptideIdentification GenAlg::startEvolution(const Size numGenerations,
			const Size endIfStableForNumGenerations, const Size bestNHits)
	{
		if(!checkEnvironment())
			return PeptideIdentification();
		PeptideIdentification pi;
		String score_name = "Normalized shared peak abundance ratio";
		pi.setScoreType(score_name);
		pi.setHigherScoreBetter(true);
		pi.setMetaValue("MZ", getPrecursorMZ());
		pi.setMetaValue("intensity", getPrecursorIntensity());
		pi.setMetaValue("charge", precursorCharge_);
		pi.setMetaValue("software","MSNovoGen");
		pi.setMetaValue("type","de novo prediction");

		initGenPool();
		Size i = 0;
		String bestSeq = "";
		Size sameBestRepeat = 0;
		while((i++ < numGenerations) && (sameBestRepeat < endIfStableForNumGenerations))
		{
			crossover();
			mutate();
			score();
			if((i % rejuvenateAfter_) == 0)
			{
				homologyKiller_->kill(*genPool_);
				genPool_->replenish(genPool_->getPreviousPoolSize());
			}
			kill();
			genPool_->sort();
			String cbseq = genPool_->getIndividual(0)->getSequence().toString();
			if(cbseq != bestSeq)
			{
			  sameBestRepeat = 0;
			  bestSeq = cbseq;
			}
			else
			  sameBestRepeat++;
		}

		std::vector<boost::shared_ptr<Chromosome> > results;
		for(std::map<String, boost::shared_ptr<Chromosome> >::const_iterator i = genPool_->getKnownIndividuals().begin(); i != genPool_->getKnownIndividuals().end(); i++) 
		{
			results.push_back(i->second);
		}
		std::sort(results.begin(), results.end(), Chromosome::sortScoreDesc);
		for(Size i = 0; (i < bestNHits && i < results.size()); i++)
		{
			PeptideHit ph(results[i]->getScore(),i,results[i]->getCharge(),results[i]->getSequence());
			pi.insertHit(ph);
		}
		return pi;
	}

	void GenAlg::setMSMSSpectrum(const MSSpectrum<> * msms)
	{
		this->msms_ = msms;
		setPrecursorInfo(this->msms_);
	}
} // namespace
