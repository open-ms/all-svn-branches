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
        : SpectraST_major(0), SpectraST_minor(0)
      {}

      SpectraSTVersion (Int maj, Int min)
        : SpectraST_major(maj), SpectraST_minor(min)
      {}

      Int SpectraST_major;
      Int SpectraST_minor;

      bool operator < (const SpectraSTVersion& v) const
      {
        if (SpectraST_major > v.SpectraST_major) return false;
        else if (SpectraST_major < v.SpectraST_major) return true;
        else // ==
        {
          if (SpectraST_minor > v.SpectraST_minor) return false;
          else if (SpectraST_minor < v.SpectraST_minor) return true;
        }
        return false;
      }
    };

    bool getVersion_(const String& version, SpectraSTVersion& SpectraST_version_i) const
    {
      // we expect two components separated by '.'
      IntList nums = IntList::create(StringList::create(version,'.'));
      if (nums.size() != 2)
      {
        return false;
      }
      SpectraST_version_i.SpectraST_major = nums[0];
      SpectraST_version_i.SpectraST_minor = nums[1];
      return true;
    }

    void registerOptionsAndFlags_()
    {
      registerInputFile_("spectrast_executable", "<executable>", "", "The 'spectrast' executable of the SpectraST installation", true, false, StringList::create("skipexists"));

      registerStringOption_("mode", "<type>", "search", "SpectraST mode.", true, false);
      vector<String> mode_strings;
      mode_strings.push_back("create");
      mode_strings.push_back("search");
      setValidStrings_("mode", mode_strings);


      addEmptyLine_();

      registerTOPPSubsection_("create","Create Mode");
      addText_("General options");

      registerOutputFile_("create:outputFileName", "<file>", "", "Output file base name");
//    registerInputFile_("cF", "<file>", "", "If <file> is not given, “spectrast_create.params“ is assumed.",false);

      registerInputFile_("create:useProbTable", "<file>", "", "Use probability table in <file>. Only those peptide ions included in the table will be imported. A probability table is a text file with one peptide ion in the format AC[160]DEFGHIK/2 per line. If a probability is supplied following the peptide ion separated by a tab, it will be used to replace the original probability of that library entr", false);
      registerInputFile_("create:useProteinList", "<file>", "", "Use protein list in <file>. Only those peptide ions associated with proteins in the list will be imported.\n A protein list is a text file with one protein identifier per line. If a number X is supplied following the protein separated by a tab, then at most X peptide ions associated with that protein will be imported.", false);
      //registerStringOption_("create:remark", "<remark>", "", "Remark. Add a Remark=<remark> comment to all library entries created.", false);
      registerStringOption_("create:annotatePeaks","<type>", "true", "Annotate peaks.", true, false);
      registerFlag_("create:binaryFormat", "Write library in binary format, which enables quicker search.", true);
      registerFlag_("create:annotatePeaks", "Annotate peaks.", true);
      registerFlag_("create:writeDtaFiles", "Write all library spectra as .dta files.", false);
      registerFlag_("create:writeMgfFiles", "Write all library spectra as one .mgf file.");

      vector<String> bool_strings;
      bool_strings.push_back("true");
      bool_strings.push_back("false");
      setValidStrings_("create:annotatePeaks", bool_strings);



      registerInputFile_("create:in", "<file>", "", "Input file ");
      setValidFormats_("create:in",StringList::create("mzML,mzXML,mzData,mgf,dta,msp"));

      

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
      // path to the log file
      String logfile(getStringOption_("log"));
      String spectrast_executable(getStringOption_("spectrast_executable"));
      String unique_name = QDir::toNativeSeparators(String(File::getTempDirectory() + "/" + File::getUniqueName()).toQString()); // body for the tmp files

      //-------------------------------------------------------------
      // parsing parameters
      //-------------------------------------------------------------

      // get version of SpectraST
      QProcess qp;
      qp.setProcessChannelMode(QProcess::MergedChannels);
      qp.start(spectrast_executable.toQString());
      bool success = qp.waitForFinished(-1);
      String output(QString(qp.readAllStandardOutput()));
      String SpectraST_version;
      SpectraSTVersion SpectraST_version_i;
      if (!success/* || qp.exitStatus() != 0 || qp.exitCode() != 0*/)
      {
        writeLog_("Warning: unable to determine the version of SpectraST - the process returned an error. Call string was: " + spectrast_executable + ". Make sure that the path to the SpectraST executable is correct!");
        return ILLEGAL_PARAMETERS;
      }
      else
      {
        vector<String> version_split;
        output.split("SpectraST (version ", version_split);
        String tmp = version_split[1];
        tmp.split(",", version_split);
        String version_line = version_split[0];

        if (getVersion_(version_line, SpectraST_version_i))
        {
          SpectraST_version = version_line;
          writeDebug_("Setting SpectraST version to " + SpectraST_version, 1);
        }
        else
        {
          writeLog_("Warning: SpectraST version output (" + output + ") not formatted as expected!");
        }
      }

      // TODO: create parameter file


      // create parameters
      StringList parameters;
      if (getStringOption_("mode") == "search")
      {
        parameters << "-s";
      } else
      {
        parameters << "-c";
      }

      // fill parameters
      parameters << "spectrast_create.params";

      QStringList qparam;
      for (Size i = 0 ; i < parameters.size();++i)
      {
        qparam << parameters[i].toQString();
      }

      writeDebug_("spectrast_executable " + parameters.concatenate(" "), 5);
      Int status = QProcess::execute(spectrast_executable.toQString(), qparam);
      if (status != 0)
      {
        writeLog_("Error: SpectraST problem! (Details can be seen in the logfile: \"" + logfile + "\")");
        writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory");

        if (getIntOption_("debug") <= 1)
        {
          // TODO: if no debuggin is enabled, delete temporary spectrast param file
          writeDebug_("Removing temporary files", 10);
          // QFile(bla.toQString()).remove();
        }
        return EXTERNAL_PROGRAM_ERROR;
      }

      // read SpectraST output
      writeDebug_("Reading output of SpectraST", 10);

      // TODO: delete temporary files
      if (getIntOption_("debug") <= 1)
      {
        writeDebug_("Removing temporary files", 10);
        // QFile(bla.toQString()).remove();
      }

      //-------------------------------------------------------------
      // convert pepXML output to idXML
      //-------------------------------------------------------------

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPSpectraSTAdapter tool;

  return tool.main(argc,argv);
}

/// @endcond
