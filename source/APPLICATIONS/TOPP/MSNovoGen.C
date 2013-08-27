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

#include <OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/GenAlg.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <vector>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

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
	  setValidFormats_("in", StringList::create("mzML"));

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

	  registerIntOption_("number_of_individuals", "", 100, "", false, false);
	  setMinInt_("number_of_individuals",1);
	  setMaxInt_("number_of_individuals",10000);

	  registerIntOption_("number_of_generations", "", 100, "", false, false);
	  setMinInt_("number_of_generations",1);
	  setMaxInt_("number_of_generations",10000);

	  registerIntOption_("return_n_best_results", "", 10, "", false, false);
	  setMinInt_("return_n_best_results",1);

	  registerIntOption_("end_after_n_stable_generations", "", 10, "", false, false);
	  setMinInt_("number_of_generations",1);
	  setMaxInt_("number_of_generations",100);

	  StringList mutators;
	  mutators.push_back("InvertingMutator");
	  mutators.push_back("SwappingMutator");
	  registerStringOption_("mutators", "", "InvertingMutator", "Mutators used by the genetic algorithm.", false, true);
	  setValidStrings_("mutators", mutators);
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

	  // identifications
	  std::vector<PeptideIdentification> pep_idents;

	  // do work
	  for(MSExperiment<>::const_iterator spec_it = exp.begin(); spec_it != exp.end(); ++spec_it)
	  {
		  const MSSpectrum<> * msms = &*spec_it;
		  GenAlg ga(msms,aaList,poolSize,pmt,fmt);
		  pep_idents.push_back(ga.startEvolution(numGenerations,endStable,bestHits));
	  }

	  std::vector<ProteinIdentification> prot_idents;
	  IdXMLFile().store(output_file,prot_idents, pep_idents);

	  return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{

  TOPPMSNovoGen tool;
  return tool.main(argc, argv);
}
