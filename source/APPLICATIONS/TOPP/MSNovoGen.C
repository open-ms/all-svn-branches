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

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <boost/shared_ptr.hpp>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/MutaterCreator.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SeederCreator.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/ScorerCreator.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/KillerCreator.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/MaterCreator.h>
#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenAlg.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <vector>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <time.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MSNovoGen MSNovoGen

  @brief MSNovoGen

  <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MRMMapper \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> </td>
          </tr>
      </table>
  </CENTER>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MSNovoGen.cli

  <B>The algorithm parameters for the Analyzer filter are:</B>
  @htmlinclude TOPP_MSNovoGen.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSNovoGen
  : public TOPPBase
{

public:
  TOPPMSNovoGen() :
    TOPPBase("MSNovoGen", "MSNovoGen is ...", true)
  {
  }

  void registerOptionsAndFlags_()
  {
	  registerInputFile_("in", "", "", "Input map", true, false);
	  setValidFormats_("in", StringList::create("mzML,dta"));

	  registerOutputFile_("out", "", "", "Output file", true, false);
	  setValidFormats_("out", StringList::create("idXML"));

      std::vector<String> all_mods;
	  ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
	  registerStringList_("modifications", "<mods>", StringList::create(""), "Modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
	  setValidStrings_("modifications", all_mods);

	  registerDoubleOption_("precursor_mass_tolerance", "", 1.5, "", false, false);
	  setMaxFloat_("precursor_mass_tolerance", 10.0);
	  setMinFloat_("precursor_mass_tolerance", 0.00001);

	  registerDoubleOption_("fragment_mass_tolerance", "", 0.3, "", false, false);
	  setMaxFloat_("fragment_mass_tolerance", 5.0);
	  setMinFloat_("fragment_mass_tolerance", 0.00001);

	  registerIntOption_("number_of_individuals", "", 500, "", false, false);
	  setMinInt_("number_of_individuals",1);
	  setMaxInt_("number_of_individuals",1000);

	  registerIntOption_("number_of_generations", "", 100, "", false, false);
	  setMinInt_("number_of_generations",1);
	  setMaxInt_("number_of_generations",1000);

	  registerIntOption_("rejuvenate_after_generations", "", 20, "", false, false);
	  setMinInt_("rejuvenate_after_generations",5);
	  setMaxInt_("rejuvenate_after_generations",1000);

	  registerIntOption_("return_n_best_results", "", 10, "", false, false);
	  setMinInt_("return_n_best_results",1);

	  registerIntOption_("retain_n_peaks", "", 100, "", false, false);
	  setMinInt_("retain_n_peaks",50);
	  setMaxInt_("retain_n_peaks",200);
	  
	  StringList peakFilter;
	  peakFilter.push_back("yes");
	  peakFilter.push_back("no");
	  registerStringOption_("use_peak_filter", "", "no", "Filter peaks (recommended if not done in pipeline).", false, true);
	  setValidStrings_("use_peak_filter", peakFilter);
	  
	  registerIntOption_("end_after_n_stable_generations", "", 10, "", false, false);
	  setMinInt_("number_of_generations",1);
	  setMaxInt_("number_of_generations",100);

	  StringList mutators;
	  mutators.push_back("InvertingMutator");
	  mutators.push_back("SubstitutingMutator");
	  mutators.push_back("SwappingMutator");
	  mutators.push_back("RandomMutator");
	  registerStringOption_("mutators", "", "SubstitutingMutator", "Mutators used by the genetic algorithm.", false, true);
	  setValidStrings_("mutators", mutators);
	  
	  StringList seeders;
	  seeders.push_back("DefaultSeeder");
	  seeders.push_back("RandomSequenceSeeder");
	  seeders.push_back("SequenceTagSeeder");
	  seeders.push_back("RandomSeeder");
	  registerStringOption_("seeders", "", "SequenceTagSeeder", "Seeders used by the genetic algorithm to seed the gene pool.", false, true);
	  setValidStrings_("seeders", seeders);
	  
	  StringList scorers;
	  scorers.push_back("NormalizedSharedAbundanceScorer");
	  scorers.push_back("HyperScorer");
	  scorers.push_back("DefaultScorer");
	  registerStringOption_("scorers", "", "NormalizedSharedAbundanceScorer", "Scorers used by the genetic algorithm as fitness function.", false, true);
	  setValidStrings_("scorers", scorers);
	  
	  StringList killers;
	  killers.push_back("SimpleDecreasingKiller");
	  killers.push_back("HomologyKiller");
	  killers.push_back("RandomKiller");
	  killers.push_back("DefaultKiller");
	  registerStringOption_("killers", "", "SimpleDecreasingKiller", "Killers used by the genetic algorithm to remove unfit individuals.", false, true);
	  setValidStrings_("killers", killers);

	  StringList maters;
	  maters.push_back("SimpleMater");
	  maters.push_back("RandomMater");
	  maters.push_back("ZipMater");
	  maters.push_back("DefaultMater");
	  registerStringOption_("maters", "", "SimpleMater", "Maters used by the genetic algorithm to perform cross over.", false, true);
	  setValidStrings_("maters", maters);
  }

  ExitCodes main_(int, const char **)
  {
	  String input_file = getStringOption_("in");
	  String output_file = getStringOption_("out");

	  std::vector<const Residue *> aaList;
	  aaList.push_back(ResidueDB::getInstance()->getResidue("A"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("R"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("N"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("D"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("C"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("E"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("Q"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("G"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("H"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("O"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("I"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("L"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("K"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("M"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("F"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("P"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("U"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("S"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("T"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("W"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("Y"));
	  aaList.push_back(ResidueDB::getInstance()->getResidue("V"));

	  StringList allowed_modifications = getStringList_("modifications");
	  for(Size i=0; i<allowed_modifications.size(); i++)
		  aaList.push_back(ResidueDB::getInstance()->getModifiedResidue(allowed_modifications[i]));

	  MSExperiment<> exp;
	  FileHandler().loadExperiment(input_file, exp);

	  Size poolSize = getIntOption_("number_of_individuals");
	  double pmt = getDoubleOption_("precursor_mass_tolerance");
	  double fmt = getDoubleOption_("fragment_mass_tolerance");
	  Size endStable = getIntOption_("end_after_n_stable_generations");
	  Size numGenerations = getIntOption_("number_of_generations");
	  Size bestHits = getIntOption_("return_n_best_results");
	  Size numPeaks = getIntOption_("retain_n_peaks");
	  Size rejAft = getIntOption_("rejuvenate_after_generations");

	  // identifications
	  std::vector<PeptideIdentification> pep_idents;
	  bool filterPeaks = false;
	  if(getStringOption_("use_peak_filter") == "yes")
		  filterPeaks = true;
	  
	  time_t s, e;
	  s = time(NULL);
	  //do work
	  for(MSExperiment<>::const_iterator spec_it = exp.begin(); spec_it != exp.end(); ++spec_it)
	  {
		  //const MSSpectrum<> * msms = &*spec_it;
		  MSSpectrum<> * msms = const_cast<MSSpectrum<> *>(&*spec_it);
		  //MSSpectrum<> * msms(const_cast<MSSpectrum<> *>(org));
		  if(filterPeaks)
		  {
		    GoodDiffFilter::create()->apply(msms);
			if(msms->size() > numPeaks)
			{
			  NLargest lf(numPeaks);
			  lf.filterPeakSpectrum(*msms);
			}
			msms->sortByPosition();
		  }
		  GenAlg ga(msms,aaList,poolSize,pmt,fmt);
	///TODO remove below
	std::map<String, boost::shared_ptr<Chromosome> > ni;
	ni.insert(std::pair<String,boost::shared_ptr<Chromosome> >("MIFAGIKK",boost::shared_ptr<Chromosome>(new Chromosome("MIFAGIKK",1))));
	ga.setKnownIndividuals(ni);
	///TODO remove above
		  ga.setRejuvenate(rejAft);
		  ga.setMutater(MutaterCreator::getInstance(getStringOption_("mutators"),ga.getPrecursorMH(),pmt,aaList));
		  ga.setSeeder(SeederCreator::getInstance(getStringOption_("seeders"),msms,ga.getPrecursorMH(),pmt,fmt,aaList));
		  ga.setScorer(ScorerCreator::getInstance(getStringOption_("scorers"),fmt));
		  ga.setKiller(KillerCreator::getInstance(getStringOption_("killers")));
		  ga.setMater(MaterCreator::getInstance(getStringOption_("maters"),ga.getPrecursorMH(),pmt,aaList));
		  pep_idents.push_back(ga.startEvolution(numGenerations,endStable,bestHits));
	  }
	  e = time(NULL);
	  double diff = difftime(e,s);

	  // search parameters
	  ProteinIdentification::SearchParameters search_parameters;
	  search_parameters.db = "none";
	  search_parameters.db_version = "";
	  search_parameters.taxonomy = "";
	  //search_parameters.charges = getStringOption_("charges");
	  search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
	  search_parameters.fixed_modifications = allowed_modifications;
	  search_parameters.enzyme = ProteinIdentification::TRYPSIN;
	  search_parameters.missed_cleavages = 1;
	  search_parameters.peak_mass_tolerance = fmt;
	  search_parameters.precursor_tolerance = pmt;

      DateTime now;	  
      String protein_identifier = "MSNOVOGEN_" + now.get();

	  ProteinIdentification protein_identification;
	  protein_identification.setDateTime(now);
	  protein_identification.setSearchEngine("MSNOVOGEN");
	  protein_identification.setSearchEngineVersion("alpha");
	  protein_identification.setSearchParameters(search_parameters);
	  protein_identification.setIdentifier(protein_identifier);
	  protein_identification.setMetaValue("runtime",String(diff)+"s");

	  std::vector<ProteinIdentification> prot_idents;
	  prot_idents.push_back(protein_identification);

	  for(std::vector<PeptideIdentification>::iterator pepit = pep_idents.begin(); pepit != pep_idents.end(); ++pepit)
	  {
		pepit->setIdentifier(protein_identifier);
	  }

	  IdXMLFile().store(output_file,prot_idents, pep_idents);

	  return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{

  TOPPMSNovoGen tool;
  return tool.main(argc, argv);
}
