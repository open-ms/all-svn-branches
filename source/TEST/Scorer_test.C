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
// $Maintainer: Jens Allmer$
// $Authors: Jens Allmer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Scorer.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenPool.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Chromosome.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
struct TestScorer :
	public Scorer
{
		TestScorer(double precursorMassTolerance) :
			Scorer(precursorMassTolerance)
		{}

	    ~TestScorer()
	    {}

	    void score(const MSSpectrum<> * msms, boost::shared_ptr<Chromosome> & chromosome) const
	    {
	    	chromosome->setScore(-1);
	    }
};
struct NullScorer :
	public Scorer
{
		NullScorer(double precursorMassTolerance) :
			Scorer(precursorMassTolerance)
		{}

		~NullScorer()
		{}

		void score(const MSSpectrum<> * msms, boost::shared_ptr<Chromosome> & chromosome) const
		{
			chromosome->setScore(0);
		}
};
START_TEST(Scorer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Scorer* ptr = 0;
Scorer* null_ptr = 0;
START_SECTION(Scorer())
{
	ptr = new TestScorer(0.5);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Scorer())
{
	delete ptr;
	ptr = new TestScorer(0.5);
}
END_SECTION

START_SECTION((Scorer(const double fragmentMassTolerance=0.5)))
{
	ptr = new TestScorer(0.5);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~Scorer()))
{
	delete ptr;
	ptr = new TestScorer(0.5);
}
END_SECTION

START_SECTION((void scorePool(const MSSpectrum<> *msms, GenPool &pool) const ))
{
	double pm = 1500.5;
	double pmt = 2;
	GenPool gp(15,pmt);
	MSSpectrum<> msms;
	Precursor prec;
	prec.setMZ(pm);
	prec.setCharge(1);
	prec.setIntensity(100.0);
	vector<Precursor> precursors;
	precursors.push_back(prec);
	msms.setPrecursors(precursors);
	NullScorer s0(0.5);
	gp.setScorer(&s0);
	gp.setMSMSSpectrum(&msms);
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("LNSFCGHWKPNAAG"),0.3)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("YIQGQTSSIGGGSMD"),0.4)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("WSVQMNKFTDMI"),0.2)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("QQPEHRGEIIPVV"),0.5)));
	gp.addIndividual(boost::shared_ptr<Chromosome>(new Chromosome(AASequence("CNGTQQFETMIND"),0.1)));
	gp.sort();
	TEST_REAL_SIMILAR(gp.getIndividual(0)->getScore(),0);
	TEST_REAL_SIMILAR(gp.getIndividual(1)->getScore(),0);
	TEST_REAL_SIMILAR(gp.getIndividual(2)->getScore(),0);
	TEST_REAL_SIMILAR(gp.getIndividual(3)->getScore(),0);
	TEST_REAL_SIMILAR(gp.getIndividual(4)->getScore(),0);
	ptr->scorePool(&msms,gp);
	TEST_REAL_SIMILAR(gp.getIndividual(0)->getScore(),-1);
	TEST_REAL_SIMILAR(gp.getIndividual(1)->getScore(),-1);
	TEST_REAL_SIMILAR(gp.getIndividual(2)->getScore(),-1);
	TEST_REAL_SIMILAR(gp.getIndividual(3)->getScore(),-1);
	TEST_REAL_SIMILAR(gp.getIndividual(4)->getScore(),-1);

}
END_SECTION

START_SECTION((virtual void score(const MSSpectrum<> *msms, boost::shared_ptr< Chromosome > &chromosome) const =0))
{
	boost::shared_ptr< Chromosome > chr(new Chromosome(AASequence("TESTER"),0,1));
	MSSpectrum<> ts;
	ptr->score(&ts,chr);
	TEST_REAL_SIMILAR(chr->getScore(),-1);
}
END_SECTION

START_SECTION((double getFragmentMassTolerance() const ))
{
  TEST_REAL_SIMILAR(ptr->getFragmentMassTolerance(),0.5);
}
END_SECTION

START_SECTION((void setFragmentMassTolerance(double fragmentMassTolerance)))
{
	TEST_REAL_SIMILAR(ptr->getFragmentMassTolerance(),0.5);
	ptr->setFragmentMassTolerance(0.4);
	TEST_REAL_SIMILAR(ptr->getFragmentMassTolerance(),0.4);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



