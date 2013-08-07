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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultScorer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DefaultScorer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DefaultScorer* ptr = 0;
DefaultScorer* null_ptr = 0;

START_SECTION((DefaultScorer(const double fragmentMassTolerance)))
{
	ptr = new DefaultScorer(0.5);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((void score(const MSSpectrum<> *msms, boost::shared_ptr< Chromosome > &chromosome) const ))
{
	MSSpectrum<> msms;
	//pm: 732.4073
	//b: 185.1285	 298.2125	 429.2530	 558.2956
	//y: 661.3702	 548.2861	 435.2020	 304.1615	 175.1190
	Peak1D p1; p1.setMZ(185.1285); p1.setIntensity(3);
	Peak1D p2; p2.setMZ(298.2125); p2.setIntensity(3);
	Peak1D p3; p3.setMZ(429.2530); p3.setIntensity(3);
	Peak1D p4; p4.setMZ(558.2956); p4.setIntensity(3);
	Peak1D p5; p5.setMZ(661.3702); p5.setIntensity(5);
	Peak1D p6; p6.setMZ(548.2861); p6.setIntensity(5);
	Peak1D p7; p7.setMZ(435.2020); p7.setIntensity(5);
	Peak1D p8; p8.setMZ(304.1615); p8.setIntensity(5);
	Peak1D p9; p9.setMZ(175.1190); p9.setIntensity(5);
	Peak1D p10; p10.setMZ(100.5000); p10.setIntensity(20);
	msms.push_back(Peak1D(p1));
	msms.push_back(Peak1D(p2));
	msms.push_back(Peak1D(p3));
	msms.push_back(Peak1D(p4));
	msms.push_back(Peak1D(p5));
	msms.push_back(Peak1D(p6));
	msms.push_back(Peak1D(p7));
	msms.push_back(Peak1D(p8));
	msms.push_back(Peak1D(p9));
	msms.push_back(Peak1D(p10));
	msms.sortByPosition();
	//4*3 + 5*5 = 37/57/6
	boost::shared_ptr< Chromosome > chr(new Chromosome(AASequence("ALLMER"),1));
	ptr->score(&msms,chr);
	TEST_REAL_SIMILAR(37.0/57.0/6.0,chr->getScore());
}
END_SECTION

START_SECTION(~DefaultScorer())
{
	delete ptr;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



