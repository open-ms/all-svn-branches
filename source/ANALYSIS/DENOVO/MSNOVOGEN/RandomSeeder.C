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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomSeeder.h>
#include <stdlib.h>
#include <time.h>
#include <boost/random.hpp>

namespace OpenMS
{

  RandomSeeder::RandomSeeder(const MSSpectrum<> * ms, const double pm, const double pmt, const double fmt, const std::vector<const Residue*> al) :
	  Seeder(ms, pm, pmt,fmt, al),
	  rss(ms, pm, pmt,fmt, al), 
	  sts(ms, pm, pmt,fmt, al), 
	  ds(ms, pm, pmt,fmt, al)
  {
	weights_.push_back(0.4);
	weights_.push_back(0.7);
	weights_.push_back(1.0);
	seed((Size)time(NULL));
  }

  boost::shared_ptr<Chromosome> RandomSeeder::createIndividual() const
  {
	boost::shared_ptr<Chromosome> chr;
	boost::random::uniform_real_distribution<double> u01;
	double rv = u01(rng);
    for(unsigned int i=0; i<weights_.size(); i++)
    {
      if(weights_[i] > rv)
      {
        switch(i) {
          case RandomSeeder::randomSequenceSeeder :
			   chr = rss.createIndividual();
			   break;
          case RandomSeeder::sequenceTagSeeder :
			   chr = sts.createIndividual();
			   break;
          default:
			   chr = ds.createIndividual();
			   break;
        }
      }
    }
    return(chr);
  }

  void RandomSeeder::seed(const Size seed)
  {
	  rss.seed(seed);
	  sts.seed(seed);
	  ds.seed(seed);
  }
} // namespace
