// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SpectraSTAdapter SpectraSTAdapter

  @brief Identifies peptides in MS/MS spectra via SpectraST (Open Mass Spectrometry Search Algorithm).

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ SpectraSTAdapter \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
    </tr>
  </table>
</CENTER>

  @em SpectraST must be installed on the system to be able to use the @em SpectraSTAdapter. See the SpectraST site
  for further information on how to download and install @em SpectraST on your system. You might find that the latest SpectraST version
  does not run on your system (to test this, run @em SpectraST from command line). If you encounter
  an error message, try another SpectraST version

  This adapter supports relative database filenames, which (when not found in the current working directory) is looked up in
  the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

  This wrapper has been tested successfully with SpectraST, version 4.x.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SpectraSTAdapter.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPSpectraSTAdapter
  : public TOPPBase
{
  public:
    TOPPSpectraSTAdapter()
      : TOPPBase("SpectraSTAdapter","Annotates MS/MS spectra using SpectraST.")
    {
    }

  protected:

    struct SpectraSTVersion
    {
      SpectraSTVersion ()
        : SpectraST_major(0), SpectraST_minor(0), SpectraST_patch(0)
      {}

      SpectraSTVersion (Int maj, Int min, Int pat)
        : SpectraST_major(maj), SpectraST_minor(min), SpectraST_patch(pat)
      {}

      Int SpectraST_major;
      Int SpectraST_minor;
      Int SpectraST_patch;

      bool operator < (const SpectraSTVersion& v) const
      {
        if (SpectraST_major > v.SpectraST_major) return false;
        else if (SpectraST_major < v.SpectraST_major) return true;
        else // ==
        {
          if (SpectraST_minor > v.SpectraST_minor) return false;
          else if (SpectraST_minor < v.SpectraST_minor) return true;
          else
          {
            return (SpectraST_patch < v.SpectraST_patch);
          }
        }

      }
    };

    bool getVersion_(const String& version, SpectraSTVersion& SpectraST_version_i) const
    {
      // we expect three components
      IntList nums = IntList::create(StringList::create(version,'.'));
      if (nums.size()!=3) return false;

      SpectraST_version_i.SpectraST_major =nums[0];
      SpectraST_version_i.SpectraST_minor =nums[1];
      SpectraST_version_i.SpectraST_patch =nums[2];
      return true;
    }

    void registerOptionsAndFlags_()
    {

      addEmptyLine_();
      addText_("(I) Create Mode");
      addText_("General options");

      registerInputFile_("in", "<file>", "", "Input file ");
      setValidFormats_("in",StringList::create("mzML,mzXML,mzData,mgf,dta,msp"));
      registerOutputFile_("out", "<file>", "", "Output file ");
      setValidFormats_("out",StringList::create("pepXML"));

      registerInputFile_("cF", "<file>", "", "If <file> is not given, “spectrast_create.params“ is assumed.",false);
      registerStringOption_("cm", "<remark>", "", "Remark. Add a Remark=<remark> comment to all library entries created.", false);
      registerInputFile_("cT", "", "", "Use probability table in <file>. Only those peptide ions included in the table will be imported. A probability table is a text file with one peptide ion in the format AC[160]DEFGHIK/2 per line. If a probability is supplied following the peptide ion separated by a tab, it will be used to replace the original probability of that library entr", false);
      registerInputFile_("cO", "", "", "Use protein list in <file>. Only those peptide ions associated with proteins in the list will be imported.\n A protein list is a text file with one protein identifier per line. If a number X is supplied following the protein separated by a tab, then at most X peptide ions associated with that protein will be imported.", false);
      

      addEmptyLine_();
      addText_("PEPXML IMPORT OPTIONS (Applicable with .pepXML files)");
      registerDoubleOption_("cP","<float>",0,"Include all spectra identified with probability no less than <prob> in the library.",false);
      registerStringOption_("cn", "<name>", "", "Specify a dataset identifier for the file to be imported.", false);
      registerFlag_("co", "Add the originating mzXML file name to the dataset identifier. Good for keeping track of in which \nMS run the peptide is observed. (Turn off with -co!)");
      registerFlag_("cg", "Set all asparagines (N) in the motif NX(S/T) as deamidated (N[115]). Use for glycocaptured peptides. (Turn off with -cg!).");
      registerFlag_("cI", "Set the instrument and acquisition settings of the spectra.");


      addEmptyLine_();
      addText_("LIBRARY MANIPULATION OPTIONS (Applicable with .splib files)");
      registerStringOption_("cf", "<pred>", "", "Filter library. Keep only those entries satisfying the predicate <pred>. \n <pred> should be a C-style predicate in quotes.", false);
      registerFlag_("cI", "Set the instrument and acquisition settings of the spectra.");
      registerFlag_("cJU", "Union. Include all the peptide ions in all the files.");
      registerFlag_("cJI", "Intersection. Only include peptide ions that are present in all the files.");
      registerFlag_("cJS", "Subtraction. Only include peptide ions in the first file that are not present in any of the other files.");
      registerFlag_("cJH", "Subtraction of homologs. Only include peptide ions in the first file \n that do not have any homologs with same charge and similar m/z in any of the other files.");
      registerFlag_("cAB", "Best replicate. Pick the best replicate of each peptide ion.");
      registerFlag_("cAC", "Consensus. Create the consensus spectrum of all replicate spectra of each peptide ion.");
      registerFlag_("cAQ", "Quality filter. Apply quality filters to library.\nIMPORTANT: Quality filter can only be applied on a SINGLE .splib file with no peptide ion represented by more than one spectrum.");
      registerFlag_("cAN", "Sort library entries by descending number of replicates used (tie-breaking by probability).");
      registerInputFile_("cD", "<file>", "", "Refresh protein mappings of each library entry against the protein database <file> (Must be in .fasta format).",false);
      registerFlag_("cu", "Delete entries whose peptide sequences do not map to any protein during refreshing with -cD option.\nWhen off, unmapped entries will be marked with Protein=0/UNMAPPED but retained in library. (Turn off with -cu!).");
      registerFlag_("cd", "Delete entries whose peptide sequences map to multiple proteins during refreshing with -cD option. (Turn off with -cd!).");


      addEmptyLine_();
      addText_("CONSENSUS/BEST-REPLICATE OPTIONS (Applicable with -cAC and -cAB options)");
      //stdwert löschen
      registerIntOption_("cr","<num>",1,"Minimum number of replicates required for each library entry.\nPeptide ions failing to have originated from enough replicates\nwill be excluded from library when creating consensus/best-replicate library.",false);


      addEmptyLine_();
      addText_("QUALITY FILTER OPTIONS (Applicable with -cAQ option)");
      //stdwert löschen
      registerIntOption_("cr","<num>",1,"Replicate quorum. Its value affects behavior of quality filter (see below).",false);
      registerIntOption_("cL","<level>",0,"Specify the stringency of the quality filter.",false);
      setMinInt_("cL", 0);
      setMaxInt_("cL", 5);
      registerIntOption_("cl","<level>",0,"-cL specifies the level for removal, -cl specifies the level for marking.\n<level> = 0: No filter.\n<level> = 1: Remove/mark impure spectra.\n<level> = 2: Also remove/mark spectra with a spectrally similar counterpart in the library that is better.\n<level> = 3: Also remove/mark inquorate entries (defined with -cr) that share no peptide sub-sequences with any other entries in the library.\n<level> = 4: Also remove/mark all singleton entries.\n<level> = 5: Also remove/mark all inquorate entries (defined with -cr).",false);
      setMinInt_("cl", 0);
      setMaxInt_("cl", 5);


      addText_("(II) Search Mode\nUsage: spectrast [ options ] <SearchFileName1> [ <SearchFileName2> ... <SearchFileNameN> ]\nwhere: SearchFileNameX = Name(s) of file containing unknown spectra to be searched.\n\t\t(Extension specifies format of file. Supports .mzXML, .mzData, .dta and .msp)");
      registerInputFile_("sF", "<file>", "", "Read search options from file.\nIf <file> is not given, ”spectrast.params” is assumed.\nNOTE: All options set in the file will be overridden by command-line options, if specified.",false);
      registerInputFile_("sL", "<file>", "", "Specify library file.\n<file> must have .splib extension. The existence of the corresponding .spidx file of the same name\nin the same directory is assumed.",false);
      registerInputFile_("sD", "<file>", "", "Specify a sequence database file.\n<file> must be in .fasta format. This will not affect the search in any way,\nbut this information will be included in the output for any downstream data processing.",false);
      registerStringOption_("sT", "<type>", "AA", "Specify the type of the sequence database file.\n<type> must be either ”AA” or ”DNA”.", false, false);
            vector<String> valid_strings;
            valid_strings.push_back("AA");
            valid_strings.push_back("DNA");
            setValidStrings_("sT", valid_strings);

      /*registerStringOption_("databaseType", "<unit>", "sTAA", "	-sTAA (default) = protein database -sTDNA = genomic database.", false, false);
      vector<String> valid_strings;
      valid_strings.push_back("sTAA");
      valid_strings.push_back("sTDNA");
      setValidStrings_("databaseType", valid_strings);*/
      registerFlag_("sR", "Cache all entries in RAM. (Turn off with -sR!)\nRequires a lot of memory (the library will usually be loaded almost in its entirety), but speeds up search for unsorted queries.");
      registerInputFile_("sS", "<file>", "", "Only search a subset of the query spectra in the search file.\nOnly query spectra with names matching a line of <file> will be searched.",false);

      addEmptyLine_();
      addText_("CANDIDATE SELECTION AND SCORING OPTIONS");
      registerDoubleOption_("sM", "<tol>", 1.5, "Specify m/z tolerance. \n<tol> is the m/z tolerance within which candidate entries\nare retrieved from the library for spectral comparision.", false);
      registerFlag_("sA", "Use isotopically averaged mass instead of monoisotopic mass. (Turn of with -sA!)");
      registerStringOption_("sC", "<type>", "ICAT_cl", "Specify the expected kind of cysteine modification for the query spectra.\n<type> must be ”ICAT_cl” for cleavable ICAT, ”ICAT_uc”for uncleavable ICAT, or ”CAM” for CarbAmidoMethyl.\nThose library spectra with a different kind of cysteine modification will be ignored.\nThe ICAT type, if any, will also be included in the pepXML output for validation by PeptideProphet.", false, false);
      vector<String> valid_strings1;
        valid_strings1.push_back("ICAT_cl");
        valid_strings1.push_back("ICAT_uc");
        valid_strings1.push_back("CAM");
        setValidStrings_("sC", valid_strings1);

      addEmptyLine_();
      addText_("Miscellaneous Options:");
      registerFlag_("V", "Verbose mode.");
      registerFlag_("Q", "Quiet mode.");
      registerOutputFile_("L", "<file>", "", "Specify name of log file. Default is ”spectrast.log”.", false);



    }

    ExitCodes main_(int , const char**)
    {
      StringList parameters;
      // path to the log file
      String logfile(getStringOption_("log"));
      String SpectraST_executable(getStringOption_("SpectraST_executable"));
      String unique_name = QDir::toNativeSeparators(String(File::getTempDirectory() + "/" + File::getUniqueName()).toQString()); // body for the tmp files
      String unique_input_name = unique_name + "_SpectraST.mgf";
      String unique_output_name = unique_name + "_SpectraST.xml";
      String unique_usermod_name = unique_name + "_SpectraST_user_mod_file.xml";

      //-------------------------------------------------------------
      // parsing parameters
      //-------------------------------------------------------------

      // get version of SpectraST
      QProcess qp;
      qp.start(SpectraST_executable.toQString(), QStringList() << "-version", QIODevice::ReadOnly); // does automatic escaping etc...
      bool success = qp.waitForFinished();
      String output (QString(qp.readAllStandardOutput ()));
      String SpectraST_version;
      SpectraSTVersion SpectraST_version_i;
      if (!success || qp.exitStatus() != 0 || qp.exitCode()!=0)
      {
        writeLog_("Warning: unable to determine the version of SpectraST - the process returned an error. Call string was: '" + SpectraST_executable + " -version'. Make sure that the path to the SpectraST executable is correct!");
        return ILLEGAL_PARAMETERS;
      }
      else
      {
        vector<String> version_split;
        output.split(' ', version_split);
        if (version_split.size() == 2 && getVersion_(version_split[1], SpectraST_version_i))
        {
          SpectraST_version = version_split[1].removeWhitespaces();
          writeDebug_("Setting SpectraST version to " + SpectraST_version, 1);
        }
        else
        {
          writeLog_("Warning: SpectraST version output (" + output + ") not formatted as expected!");
        }
      }
      // parse arguments
      String inputfile_name = getStringOption_("in");
      String outputfile_name = getStringOption_("out");
      String db_name = String(getStringOption_("database"));
      // @todo: find DB for SpectraST (if not given) in OpenMS_bin/share/OpenMS/DB/*.fasta|.pin|...


      if (db_name.suffix('.') != "psq")
      {
        db_name += ".psq";
      }

      if (!File::readable(db_name))
      {
        String full_db_name;
        try
        {
          full_db_name = File::findDatabase(db_name);
        }
        catch (...)
        {
          printUsage_();
          return ILLEGAL_PARAMETERS;
        }
        db_name = full_db_name;
      }

      db_name = db_name.substr(0,db_name.size()-4); // SpectraST requires the filename without the .psq part

      parameters << "-d"  << String(db_name);
      parameters << "-to" << String(getDoubleOption_("fragment_mass_tolerance")); //String(getDoubleOption_("to"));
      parameters << "-hs" << String(getIntOption_("hs"));
      parameters << "-te" << String(getDoubleOption_("precursor_mass_tolerance")); //String(getDoubleOption_("te"));
      if (getFlag_("precursor_mass_tolerance_unit_ppm"))
      {
        if (SpectraST_version_i < SpectraSTVersion(2,1,8))
        {
          writeLog_("This SpectraST version (" + SpectraST_version + ") does not support the 'precursor_mass_tolerance_unit_ppm' flag."
                   +" Please disable it and set the precursor tolerance in Da."
                   +" Required version is 2.1.8 and above.\n");
          return ILLEGAL_PARAMETERS;
        }
        parameters << "-teppm"; // only from SpectraST 2.1.8 on
      }
      parameters << "-zl" << String(getIntOption_("min_precursor_charge")); //String(getIntOption_("zl"));
      parameters << "-zh" <<  String(getIntOption_("max_precursor_charge")); //String(getIntOption_("zh"));
      parameters << "-zt" <<  String(getIntOption_("zt"));
      parameters << "-zc" <<  String(getIntOption_("zc"));
      parameters << "-zcc" << String(getIntOption_("zcc"));
      parameters << "-zoh" << String(getIntOption_("zoh"));
      parameters << "-no" << String(getIntOption_("no"));
      parameters << "-nox" << String(getIntOption_("nox"));
      parameters << "-sp" << String(getIntOption_("sp"));
      parameters << "-sb1" << String(getIntOption_("sb1"));
      parameters << "-sct" << String(getIntOption_("sct"));
      parameters << "-x" << getStringOption_("x");
      parameters << "-hl" << String(getIntOption_("hl"));
      parameters << "-hm" << String(getIntOption_("hm"));
      parameters << "-ht" <<  String(getIntOption_("ht"));
      parameters << "-tex" << String(getDoubleOption_("tex"));
      parameters << "-i" << getStringOption_("i");
      parameters << "-z1" << String(getDoubleOption_("z1"));
      parameters << "-v" << String(getIntOption_("v"));
      parameters << "-e" << String(getIntOption_("e"));
      parameters << "-tez" << String(getIntOption_("tez"));


      parameters << "-tom" << String(getIntOption_("tom"));
      parameters << "-tem" << String(getIntOption_("tem"));

      parameters << "-mm" << String(getIntOption_("mm"));
      parameters << "-is" << String(getDoubleOption_("is"));
      parameters << "-ir" << String(getDoubleOption_("ir"));
      parameters << "-ii" << String(getDoubleOption_("ii"));
      parameters << "-nt" << String(getIntOption_("threads"));

      if (getFlag_("mnm"))
      {
        parameters << "-mnm";
      }

      parameters << "-fm" << unique_input_name;
      parameters << "-ox" << unique_output_name;

      if (getIntOption_("debug") == 0)
      {
        parameters << "-ni";
      }
      parameters << "-he" << String(getDoubleOption_("he"));


      // read mapping for the modifications
      String file = File::find("CHEMISTRY/SpectraST_modification_mapping");

      TextFile infile(file);
      Map<String, UInt> mods_map;
      for (TextFile::ConstIterator it = infile.begin(); it != infile.end(); ++it)
      {
        vector<String> split;
        it->split(',', split);

        if (it->size() > 0 && (*it)[0] != '#')
        {
          if (split.size() < 2)
          {
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parse mapping file line: '" + *it + "'", "");
          }
          vector<ResidueModification> mods;
          for (Size i = 2; i != split.size(); ++i)
          {
            String tmp(split[i].trim());
            if (!tmp.empty())
            {
              mods_map[tmp] = split[0].trim().toInt();
            }
          }
        }
      }

      writeDebug_("Evaluating modifications", 1);
      ModificationDefinitionsSet mod_set(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
      writeDebug_("Setting modifications", 1);
      UInt user_mod_num(119);
      vector<pair<UInt, String> > user_mods;
      // fixed modifications
      if ( !getStringList_("fixed_modifications").empty())
      {
        set<String> mod_names = mod_set.getFixedModificationNames();
        StringList mod_list;
        for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
        {
          if (mods_map.has(*it))
          {
            mod_list.push_back(String(mods_map[*it]));
          }
          else
          {
            mod_list.push_back(String(user_mod_num));
            // add this to the usermods
            user_mods.push_back(make_pair(user_mod_num++, *it));
            writeDebug_("Inserting unknown fixed modification: '" + *it + "' into SpectraST", 1);
          }
        }
        if (mod_list.size() > 0)
        {
          parameters << "-mf" << mod_list.concatenate(",");
        }
      }

      if ( !getStringList_("variable_modifications").empty() )
      {
        set<String> mod_names = mod_set.getVariableModificationNames();
        StringList mod_list;

        for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
        {
          if (mods_map.has(*it))
          {
            mod_list.push_back(String(mods_map[*it]));
          }
          else
          {
            mod_list.push_back(String(user_mod_num));
            // add this to the usermods
            user_mods.push_back(make_pair(user_mod_num++, *it));
            writeDebug_("Inserting unknown variable modification: '" + *it + "' into SpectraST", 1);
          }
        }

        if (mod_list.size() > 0)
        {
          parameters << "-mv" << mod_list.concatenate(",");
        }
      }

      String additional_user_mods_filename = getStringOption_("SpectraST_user_mods");
      // write unknown modifications to user mods file
      if ( !user_mods.empty() || additional_user_mods_filename != "")
      {
        writeDebug_("Writing usermod file to " + unique_usermod_name, 1);
        parameters << "-mux" << File::absolutePath(unique_usermod_name);
        ofstream out(unique_usermod_name.c_str());
        out << "<?xml version=\"1.0\"?>" << endl;
        out << "<MSModSpecSet xmlns=\"http://www.ncbi.nlm.nih.gov\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\" xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov SpectraST.xsd\">" << endl;

        UInt user_mod_count(1);
        for (vector<pair<UInt, String> >::const_iterator it = user_mods.begin(); it != user_mods.end(); ++it)
        {
          writeDebug_("Writing information into user mod file of modification: " + it->second, 1);
          out << "<MSModSpec>" << endl;
          out << "\t<MSModSpec_mod>" << endl;
          out << "\t\t<MSMod value=\"usermod" << user_mod_count++ << "\">" << it->first << "</MSMod>" << endl;
          out << "\t</MSModSpec_mod>" << endl;
          out << "\t<MSModSpec_type>" << endl;

          /*
              0 modaa	-  at particular amino acids
              1 modn	-  at the N terminus of a protein
              2 modnaa	-  at the N terminus of a protein at particular amino acids
              3 modc	-  at the C terminus of a protein
              4 modcaa	-  at the C terminus of a protein at particular amino acids
              5 modnp	-  at the N terminus of a peptide
              6 modnpaa	-  at the N terminus of a peptide at particular amino acids
              7 modcp	-  at the C terminus of a peptide
              8 modcpaa	-  at the C terminus of a peptide at particular amino acids
              9 modmax	-  the max number of modification types
          */

          ResidueModification::Term_Specificity ts = ModificationsDB::getInstance()->getModification(it->second).getTermSpecificity();
          String origin = ModificationsDB::getInstance()->getModification(it->second).getOrigin();
          if (ts == ResidueModification::ANYWHERE)
          {
            out << "\t\t<MSModType value=\"modaa\">0</MSModType>" << endl;
          }
          if (ts == ResidueModification::C_TERM)
          {
            if (origin == "" || origin == "X")
            {
              out << "\t\t<MSModType value=\"modcp\">7</MSModType>" << endl;
            }
            else
            {
              out << "\t\t<MSModType value=\"modcpaa\">8</MSModType>" << endl;
            }
          }
          if (ts == ResidueModification::N_TERM)
          {
            if (origin == "" || origin == "X")
            {
              out << "\t\t<MSModType value=\"modnp\">5</MSModType>" << endl;
            }
            else
            {
              out << "\t\t<MSModType value=\"modnpaa\">6</MSModType>" << endl;
            }
          }
          out << "\t</MSModSpec_type>" << endl;

          out << "\t<MSModSpec_name>" << it->second << "</MSModSpec_name>" << endl;
          out << "\t<MSModSpec_monomass>" << ModificationsDB::getInstance()->getModification(it->second).getDiffMonoMass()  << "</MSModSpec_monomass>" << endl;
          out << "\t<MSModSpec_averagemass>" << ModificationsDB::getInstance()->getModification(it->second).getDiffAverageMass() << "</MSModSpec_averagemass>" << endl;
          out << "\t<MSModSpec_n15mass>0</MSModSpec_n15mass>" << endl;

          if (origin != "")
          {
            out << "\t<MSModSpec_residues>" << endl;
            out << "\t\t<MSModSpec_residues_E>" << origin << "</MSModSpec_residues_E>" << endl;
            out << "\t</MSModSpec_residues>" << endl;

            /* TODO: Check why these are always 0
            DoubleReal neutral_loss_mono = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossMonoMass();
            DoubleReal neutral_loss_avg = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossAverageMass();
            */
            DoubleReal neutral_loss_mono = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossDiffFormula().getMonoWeight();
            DoubleReal neutral_loss_avg = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossDiffFormula().getAverageWeight();

            if (fabs(neutral_loss_mono) > 0.00001)
            {
              out << "\t<MSModSpec_neutralloss>" << endl;
              out << "\t\t<MSMassSet>" << endl;
              out << "\t\t\t<MSMassSet_monomass>" << neutral_loss_mono << "</MSMassSet_monomass>" << endl;
              out << "\t\t\t<MSMassSet_averagemass>" << neutral_loss_avg << "</MSMassSet_averagemass>" << endl;
              out << "\t\t\t<MSMassSet_n15mass>0</MSMassSet_n15mass>" << endl;
              out << "\t\t</MSMassSet>" << endl;
              out << "\t</MSModSpec_neutralloss>" << endl;
            }

            out << "</MSModSpec>" << endl;
          }
        }

        // Add additional MSModSPec subtrees to generated user mods
        ifstream additional_user_mods_file(additional_user_mods_filename.c_str());
        String line;
        if(additional_user_mods_file.is_open())
        {
          while (additional_user_mods_file.good())
          {
            getline(additional_user_mods_file, line);
            out << line << endl;
          }
          additional_user_mods_file.close();
        }
        out << "</MSModSpecSet>" << endl;
        out.close();
      }

      //-------------------------------------------------------------
      // reading input
      //-------------------------------------------------------------

      FileHandler fh;
      FileTypes::Type in_type = fh.getType(inputfile_name);
      PeakMap map;
      fh.getOptions().addMSLevel(2);
      fh.loadExperiment(inputfile_name, map, in_type, log_type_);

      ProteinIdentification protein_identification;
      vector<PeptideIdentification> peptide_ids;

      writeDebug_("Read " + String(map.size()) + " spectra from file", 5);

      vector<ProteinIdentification> protein_identifications;
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      writeDebug_("Storing input file: " + unique_input_name, 5);
      MascotInfile SpectraST_infile;
      SpectraST_infile.store(unique_input_name, map, "SpectraST search tmp file");

      // @todo find SpectraST if not given
      // executable is stored in OpenMS_bin/share/OpenMS/3rdParty/SpectraST/SpectraSTcl(.exe)
      // or PATH

      QStringList qparam;
      for (Size i=0; i<parameters.size();++i)
      {
        qparam << parameters[i].toQString();
      }
      writeDebug_("SpectraST_executable " + parameters.concatenate(" "), 5);
      Int status = QProcess::execute(SpectraST_executable.toQString(), qparam);
      if (status != 0)
      {
        writeLog_("Error: SpectraST problem! (Details can be seen in the logfile: \"" + logfile + "\")");
        writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory");
        if (getIntOption_("debug") <= 1)
        {
          QFile(unique_input_name.toQString()).remove();
          QFile(unique_output_name.toQString()).remove();
        }
        if ( !user_mods.empty() || additional_user_mods_filename!="")
        {
          QFile(unique_usermod_name.toQString()).remove();
        }
        return EXTERNAL_PROGRAM_ERROR;
      }

      // read SpectraST output
      writeDebug_("Reading output of SpectraST", 10);
//      SpectraSTXMLFile SpectraST_out_file;
//      SpectraST_out_file.setModificationDefinitionsSet(mod_set);  // TODO: add modifications from additional user mods subtree
//      SpectraST_out_file.load(unique_output_name, protein_identification, peptide_ids);

      // SpectraST does not write fixed modifications so we need to add them to the sequences
      set<String> fixed_mod_names = mod_set.getFixedModificationNames();
      vector<String> fixed_nterm_mods, fixed_cterm_mods;
      Map<String, String> fixed_residue_mods;
      writeDebug_("Splitting modification into N-Term, C-Term and anywhere specificity", 1);
      for (set<String>::const_iterator it = fixed_mod_names.begin(); it != fixed_mod_names.end(); ++it)
      {
        ResidueModification::Term_Specificity ts = ModificationsDB::getInstance()->getModification(*it).getTermSpecificity();
        if (ts == ResidueModification::ANYWHERE)
        {
          fixed_residue_mods[ModificationsDB::getInstance()->getModification(*it).getOrigin()] = *it;
        }
        if (ts == ResidueModification::C_TERM)
        {
          fixed_cterm_mods.push_back(*it);
        }
        if (ts == ResidueModification::N_TERM)
        {
          fixed_nterm_mods.push_back(*it);
        }
      }
      writeDebug_("Assigning modifications to peptides", 1);
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        vector<PeptideHit> hits = it->getHits();
        for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
        {
          AASequence seq = pit->getSequence();
          for (vector<String>::const_iterator mit = fixed_nterm_mods.begin(); mit != fixed_nterm_mods.end(); ++mit)
          {
            seq.setNTerminalModification(*mit);
          }
          for (vector<String>::const_iterator mit = fixed_cterm_mods.begin(); mit != fixed_cterm_mods.end(); ++mit)
          {
            seq.setCTerminalModification(*mit);
          }
          UInt pos = 0;
          for (AASequence::Iterator mit = seq.begin(); mit != seq.end(); ++mit, ++pos)
          {
            if (fixed_residue_mods.has(mit->getOneLetterCode()))
            {
              seq.setModification(pos, fixed_residue_mods[mit->getOneLetterCode()]);
            }
          }
          pit->setSequence(seq);
        }
        it->setHits(hits);
      }

      // delete temporary files
      if (getIntOption_("debug") <= 1)
      {
        writeDebug_("Removing temporary files", 10);
        QFile(unique_input_name.toQString()).remove();
        QFile(unique_output_name.toQString()).remove();
        if ( !user_mods.empty() )
        {
          QFile(unique_usermod_name.toQString()).remove();
        }
      }

      // handle the search parameters
      ProteinIdentification::SearchParameters search_parameters;
      search_parameters.db = getStringOption_("database");
      //search_parameters.db_version =
      search_parameters.taxonomy = getStringOption_("x");
      search_parameters.charges = "+" + String(getIntOption_("min_precursor_charge")) + "-+" + String(getIntOption_("max_precursor_charge"));
      ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;

      if (getIntOption_("tom") == 1)
      {
        mass_type = ProteinIdentification::AVERAGE;
      }
      else
      {
        if (getIntOption_("tom") != 0)
        {
          writeLog_("Warning: unrecognized mass type: " + String(getIntOption_("tom")));
        }
      }
      search_parameters.mass_type = mass_type;
      search_parameters.fixed_modifications = getStringList_("fixed_modifications");
      search_parameters.variable_modifications = getStringList_("variable_modifications");
      ProteinIdentification::DigestionEnzyme enzyme = ProteinIdentification::TRYPSIN;

      UInt e(getIntOption_("e"));
      if (e != 0)
      {
        writeLog_("Warning: cannot handle enzyme: " + getIntOption_("e"));
      }

      search_parameters.enzyme = enzyme;
      search_parameters.missed_cleavages = getIntOption_("v");
      search_parameters.peak_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
      search_parameters.precursor_tolerance = getDoubleOption_("precursor_mass_tolerance");


      protein_identification.setSearchParameters(search_parameters);
      protein_identification.setSearchEngineVersion(SpectraST_version);
      protein_identification.setSearchEngine("SpectraST");

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      protein_identifications.push_back(protein_identification);
      IdXMLFile().store(outputfile_name, protein_identifications, peptide_ids);

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPSpectraSTAdapter tool;

  return tool.main(argc,argv);
}

/// @endcond
