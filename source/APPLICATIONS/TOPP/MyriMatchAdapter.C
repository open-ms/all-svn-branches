// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Dimitri Schachmann $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <stdio.h>

using namespace OpenMS;
using namespace std;

class MyriMatchAdapter
  : public TOPPBase
{
public:
  MyriMatchAdapter()
    : TOPPBase("MyriMatchAdapter","Annotates MS/MS spectra using MyriMatch.",false)
  {
  }

protected :

  struct MyriMatchVersion
  {
    MyriMatchVersion ()
      : myrimatch_major(0), myrimatch_minor(0), myrimatch_patch(0)
    {}

    MyriMatchVersion (Int maj, Int min, Int pat)
      : myrimatch_major(maj), myrimatch_minor(min), myrimatch_patch(pat)
    {}

    Int myrimatch_major;
    Int myrimatch_minor;
    Int myrimatch_patch;

    bool operator < (const MyriMatchVersion& v) const
    {
      if (myrimatch_major > v.myrimatch_major) return false;
      else if (myrimatch_major < v.myrimatch_major) return true;
      else
        {
          if (myrimatch_minor > v.myrimatch_minor) return false;
          else if (myrimatch_minor < v.myrimatch_minor) return true;
          else
            {
              return (myrimatch_patch < v.myrimatch_patch);
            }
        }

    }

  };

  bool getVersion_(const String& version, MyriMatchVersion& myrimatch_version_i) const
  {
    // we expect three components
    IntList nums = IntList::create(StringList::create(version,'.'));
    if (nums.size()!=3) return false;

    myrimatch_version_i.myrimatch_major =nums[0];
    myrimatch_version_i.myrimatch_minor =nums[1];
    myrimatch_version_i.myrimatch_patch =nums[2];
    return true;
  }


  void registerOptionsAndFlags_()
  {
      addEmptyLine_();
      addText_("Common Identification engine options");

      registerInputFile_("in", "<file>", "", "Input file ");
      setValidFormats_("in",StringList::create("mzML,mzXML")); //TODO: forbid mzXML
      registerOutputFile_("out", "<file>", "", "Output file ");
      setValidFormats_("out",StringList::create("idXML"));
      registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "Precursor mono mass tolerance.", false);

      registerStringOption_("precursor_mass_tolerance_unit", "<unit>", "Da", "Unit to be used for precursor mass tolerance.",false);
      setValidStrings_("precursor_mass_tolerance_unit", StringList::create("Da,ppm"));

      registerFlag_("precursor_mass_tolerance_avg", "If this flag is set, the average mass is used in the precursor mass tolerance.");
      registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "Fragment mass error in Dalton", false);

      registerStringOption_("fragment_mass_tolerance_unit", "<unit>", "Da", "Unit to be used for fragment mass tolerance.",false);
      setValidStrings_("fragment_mass_tolerance_unit", StringList::create("Da,ppm"));

      registerInputFile_("database", "<fasta-file>", "",
                         "NCBI formatted FASTA files. Only the .fasta filename should be given.",
                         true, false);
      registerIntOption_("min_precursor_charge", "<charge>", 1, "Minimum precursor ion charge", false);
      registerIntOption_("max_precursor_charge", "<charge>", 3, "Maximum precursor ion charge", false);
      vector<String> all_mods;
      ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
      registerStringList_("fixed_modifications", "<mods>", StringList::create(""),
                          "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
      setValidStrings_("fixed_modifications", all_mods);
      registerStringList_("variable_modifications", "<mods>", StringList::create(""),
                          "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
      setValidStrings_("variable_modifications", all_mods);

      addEmptyLine_();
      addText_("MyriMatch specific input options");

      registerInputFile_("myrimatch_executable", "<executable>", "myrimatch",
                         "The 'myrimatch' executable of the MyriMatch installation", true, false, StringList::create("skipexists"));
      registerIntOption_("NumChargeStates", "<num>", 3, "The number of charge states that MyriMatch will handle during all stages of the program.", false);
      registerDoubleOption_("TicCutoffPercentage", "<percentage>", 0.98, "Noise peaks are filtered out by sorting the original peaks in descending order of intensity, and then picking peaks from that list until the cumulative ion current of the picked peaks divided by the total ion current (TIC) is greater than or equal to this parameter.", false);
      registerIntList_("MonoisotopeAdjustmentSet","<set>",IntList::create(""),"This parameter defines a set of isotopes (0 being the instrument-called monoisotope) to try as the monoisotopic precursor m/z. To disable this technique, set the value to '0'.",false);
      registerIntOption_("MaxDynamicMods", "<num>", 2, "This parameter sets the maximum number of modified residues that may be in any candidate sequence.", false);
      registerIntOption_("MaxResultRank", "<rank>", 5, "This parameter sets the maximum rank of peptide-spectrum-matches to report for each spectrum.", false);
      registerStringOption_("CleavageRules", "<rule>", "", "This parameter allows the user to control the way peptides are generated from the protein database.", false);
      vector<String> all_rules;
      all_rules.push_back("Trypsin");
      all_rules.push_back("Trypsin/P");
      all_rules.push_back("Chymotrypsin");
      all_rules.push_back("TrypChymo");
      all_rules.push_back("Lys-C");
      all_rules.push_back("Lys-C/P");
      all_rules.push_back("Asp-N");
      all_rules.push_back("PepsinA");
      all_rules.push_back("CNBr");
      all_rules.push_back("Formic_acid");
      all_rules.push_back("NoEnzyme");
      setValidStrings_("CleavageRules", all_rules);

      registerIntOption_("MinTerminiCleavages", "<num>", 2, "By default, when generating peptides from the protein database, a peptide must start and end at a valid cleavage site. Setting this parameter to 0 or 1 will reduce that requirement, so that neither terminus or only one terminus of the peptide must match one of the cleavage rules specified in the CleavageRules parameter. This parameter is useful to turn a tryptic digest into a semi-tryptic digest.", false); // TODO: Description copied from MM doc
      registerIntOption_("MaxMissedCleavages", "<num>", -1, "By default, when generating peptides from the protein database, a peptide may contain any number of missed cleavages. A missed cleavage is a site within the peptide that matches one of the cleavage rules (refer to CleavageRules). Settings this parameter to some other number will stop generating peptides from a sequence if it contains more than the specified number of missed cleavages.", false); // TODO: Description copied from MM doc

      // advanced options
      registerDoubleOption_("MinPeptideMass", "<mass>", 0.0, "When preprocessing the experimental spectra, any spectrum with a precursor mass that is less than the specified mass will be disqualified.", false, true); // TODO: Description copied from MM doc
      registerDoubleOption_("MaxPeptideMass", "<mass>", 10000.0, "When preprocessing the experimental spectra, any spectrum with a precursor mass that exceeds the specified mass will be disqualified.", false, true); // TODO: Description copied from MM doc
      registerIntOption_("MinPeptideLength", "<length>", 5, "When digesting proteins, any peptide which does not meet or exceed the specified length will be disqualified.", false, true); // TODO: Description copied from MM doc
      registerIntOption_("MaxPeptideLength", "<length>", 75, "When digesting proteins, any peptide which exceeds this specified length will be disqualified.", false, true); // TODO: Description copied from MM doc
      registerFlag_("UseSmartPlusThreeModel", "When this parameter is set, then for each peptide bond, an internal calculation is done to estimate the basicity of the b and y fragment sequence. The precursors protons are distributed to those ions based on that calculation, with the more basic sequence generally getting more of the protons..",true); // TODO: Description copied from MM doc

      registerIntOption_("ProteinSampleSize", "<size>", 100, "Before beginning sequence candidate generation and scoring, MyriMatch will do a random sampling of the protein database to get an estimate of the number of comparisons that will be done by the job.", false, true); // TODO: Description copied from MM doc
      registerIntOption_("NumIntensityClasses", "<num>", 3, "Before scoring any candidates, experimental spectra have their peaks stratified into the number of intensity classes specified by this parameter.", false, true); // TODO: Description copied from MM doc
      registerDoubleOption_("ClassSizeMultiplier", "<factor>", 2.0, "When stratifying peaks into a specified, fixed number of intensity classes, this parameter controls the size of each class relative to the class above it (where the peaks are more intense). ", false, true); // TODO: Description copied from MM doc

  }


  ExitCodes main_(int , const char**)
  {
    StringList parameters;
    String logfile(getStringOption_("log"));
    String myrimatch_executable(getStringOption_("myrimatch_executable"));

    //-------------------------------------------------------------
    // get version of MyriMatch
    //-------------------------------------------------------------

    QProcess qp;
    qp.start(myrimatch_executable.toQString(), QStringList(), QIODevice::ReadOnly); // does automatic escaping etc...
    qp.waitForFinished();
    String output (QString(qp.readAllStandardOutput ()));
    String myrimatch_version;
    MyriMatchVersion myrimatch_version_i;
    String tmp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/").toQString()); // body for the tmp files

    vector<String> lines;
    vector<String> version_split;
    output.split('\n',lines);
    lines[1].split(' ',version_split);

    if (version_split.size() == 3 && getVersion_(version_split[1], myrimatch_version_i))
      {
        myrimatch_version = version_split[1].removeWhitespaces();
        writeDebug_("Setting MyriMatch version to " + myrimatch_version, 1);
      }
    else
      {
        writeLog_("Warning: MyriMatch version output (" + output + ") not formatted as expected!");
      }
    if(myrimatch_version_i.myrimatch_major != 2 && myrimatch_version_i.myrimatch_minor != 1)
      {
        writeDebug_("Warning: unsupported MyriMatch version (" + myrimatch_version + "). Tested only for MyriMatch 2.1.x",0);
      }


    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String inputfile_name = getStringOption_("in");
    writeDebug_(inputfile_name,0);
    String outputfile_name = getStringOption_("out");
    String db_name = String(getStringOption_("database"));
    FileHandler fh;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> peptide_identifications;

    // building parameter String

    // Common Identification engine options

    // translating unimod notation to MyriMatch notation of PTMs.
    ModificationDefinitionsSet mod_set(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
    if ( !getStringList_("fixed_modifications").empty())
      {
        set<String> mod_names = mod_set.getFixedModificationNames();
        StringList mod_list;
        for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
          {
            ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
            String origin = String(mod.getOrigin());
            String mass_diff = String(mod.getDiffMonoMass());
            printf("origin for %s ist %s\n",(*it).c_str(), origin.c_str());
            if(origin == "N-term")
              {
                origin = "(";
              }
            else if (origin == "C-term")
              {
                origin = ")";
              }
            else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "N-term")
              {
                origin = "(" + origin;
              }
            else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "C-term")
              {
                origin = ")" + origin;
              }
            mod_list.push_back(origin + " " + mod.getDiffMonoMass());
          }
        if (mod_list.size() > 0)
          {
            parameters << "-StaticMods" << mod_list.concatenate(" ");
          }
      }

    if ( !getStringList_("variable_modifications").empty())
      {
        char mod_chars [] =  {'@','^','$'};
        int mod_num = sizeof(mod_chars) / sizeof(mod_chars[0]);
        int mod_count = 0;
        set<String> mod_names = mod_set.getVariableModificationNames();
        StringList mod_list;
        for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
          {
            if(mod_count < mod_num)
              {
                ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
                String origin = String(mod.getOrigin());
                String mass_diff = String(mod.getDiffMonoMass());
                if(origin == "N-term")
                  {
                    origin = "(";
                  }
                else if (origin == "C-term")
                  {
                    origin = ")";
                  }
                else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "N-term")
                  {
                    origin = "(" + origin;
                  }
                else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "C-term")
                  {
                    origin = ")" + origin;
                  }
                mod_list.push_back(origin + " " + String(mod_chars[mod_count++]) + " " + mass_diff);
              }
            else
              {
                writeDebug_("Too many variable modifications provided. Please provied no more than " + String(mod_num),0);
                break;
              }
          }
        if (mod_list.size() > 0)
          {
            parameters << "-DynamicMods" << mod_list.concatenate(" ");
            writeDebug_("DynamicMods:",0); // TODO: delete line
            writeDebug_(mod_list.concatenate(" "),0); // TODO: delete line
          }
      }

    parameters << "-ProteinDatabase"  << db_name;

    String precursor_mass_tolerance_unit = getStringOption_("precursor_mass_tolerance_unit") == "Da" ? " m/z" : " ppm";

    if(getFlag_("precursor_mass_tolerance_avg"))
      {
        parameters << "-AvgPrecursorMzTolerance" << String(getDoubleOption_("precursor_mass_tolerance")) + precursor_mass_tolerance_unit;
      }
    else
      {
        parameters << "-MonoPrecursorMzTolerance" << (String(getDoubleOption_("precursor_mass_tolerance")) + precursor_mass_tolerance_unit);
      }

    String fragment_mass_tolerance_unit = getStringOption_("fragment_mass_tolerance_unit");
    if(fragment_mass_tolerance_unit == "Da")
    {
      fragment_mass_tolerance_unit = "m/z";
    }

    parameters << "-FragmentMzTolerance" << String(getDoubleOption_("fragment_mass_tolerance")) + " " + fragment_mass_tolerance_unit;
    int min_charge = getIntOption_("min_precursor_charge");
    int max_charge = getIntOption_("max_precursor_charge");
    parameters << "-SpectrumListFilters" << "chargeStatePredictor false " +  String(max_charge) + " " +  String(min_charge) + " 0.9"; // Mal gucken
    parameters << "-ThreadCountMultiplier" << String(getIntOption_("threads")); // MyriMatch does not recognise this, even though it's in the manual.


    // MyriMatch specific parameters
    parameters << "-NumChargeStates" << getIntOption_("NumChargeStates");
    parameters << "-TicCutoffPercentage" << String(getDoubleOption_("TicCutoffPercentage"));
    parameters << "-MaxDynamicMods" << getIntOption_("MaxDynamicMods");
    parameters << "-MaxResultRank" << getIntOption_("MaxResultRank");
    parameters << "-MinTerminiCleavages" << getIntOption_("MinTerminiCleavages");
    parameters << "-MaxMissedCleavages" << getIntOption_("MaxMissedCleavages");
    String cleavage_rule = getStringOption_("CleavageRules");
    if(cleavage_rule.empty())
      {
        cleavage_rule = "Trypsin/P";
      }
    parameters << "-CleavageRules" << cleavage_rule;

    // advanced parameters
    parameters << "-MinPeptideMass"   << getDoubleOption_("MinPeptideMass");
    parameters << "-MaxPeptideMass"   << getDoubleOption_("MaxPeptideMass");
    parameters << "-MinPeptideLength" << getIntOption_("MinPeptideLength");
    parameters << "-MaxPeptideLength" << getIntOption_("MaxPeptideLength");
    parameters << "-ProteinSampleSize" << getIntOption_("ProteinSampleSize");
    parameters << "-NumIntensityClasses" << getIntOption_("NumIntensityClasses");
    parameters << "-ClassSizeMultiplier" << getDoubleOption_("ClassSizeMultiplier");

    // getIntList gibt nur das erste Element der Liste her
//    IntList    adjustment_int_list = getIntList_("MonoisotopeAdjustmentSet");
//    //adjustment_int_list = IntList::create("1,2,3,4");
//    set<Int>   adjustment_int_set(adjustment_int_list.begin(),adjustment_int_list.end());
//    StringList adjustment_str_list = StringList();
//    for (set<Int>::const_iterator it = adjustment_int_set.begin(); it != adjustment_int_set.end(); ++it)
//      {
//        adjustment_str_list << String(*it);
//        writeDebug_("wtf?",0);
//      }
//    String adjustment_set_string = "[";
//    adjustment_set_string += adjustment_str_list.concatenate(",");
//    adjustment_set_string += "]";
//    parameters << "-MonoisotopeAdjustmentSet" << adjustment_set_string;

    // Constant parameters
    parameters << "-DecoyPrefix" << "rev_";
    //    parameters << "-UseMultipleProcessors" << String(true);

    // path to inputfile must be the last parameter
    //    writeDebug_(inputfile_name,0)

    parameters << inputfile_name;

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    QStringList qparam;
    for (Size i=0; i<parameters.size();++i)
      {
        qparam << parameters[i].toQString();
        writeDebug_(parameters[i].toQString(), 0);
      }

    QProcess process;
    process.setWorkingDirectory(tmp_dir.toQString());

    process.start(myrimatch_executable.toQString(), qparam, QIODevice::ReadOnly);
    String myri_msg (QString(process.readAllStandardOutput ()));
    bool success = process.waitForFinished(-1);
    if (!success)
      {
        writeLog_("Error: MyriMatch problem! (Details can be seen in the logfile: \"" + logfile + "\")");
        writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory");
        return EXTERNAL_PROGRAM_ERROR;
      }

    //-------------------------------------------------------------
    // reading MyriMatch output
    //-------------------------------------------------------------

    writeDebug_("Reading output of MyriMatch", 5);
    // String exp_name = inputfile_name;
    String exp_name = File::basename(inputfile_name);
    String pep_file =  tmp_dir + File::removeExtension(exp_name)+".pepXML";
    //String pep_file =  File::removeExtension(exp_name)+".pepXML";
    bool use_precursor_data = true;
    MSExperiment<> exp;

    fh.loadExperiment(inputfile_name, exp);

    PepXMLFile().load(pep_file, protein_identifications, peptide_identifications,
                      exp_name, exp, use_precursor_data);

    //QFile(pep_file.toQString()).remove();
    //-------------------------------------------------------------
    // writng results
    //-------------------------------------------------------------

    IdXMLFile().store(outputfile_name, protein_identifications, peptide_identifications);
    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{
  MyriMatchAdapter tool;
  return tool.main(argc,argv);
}

