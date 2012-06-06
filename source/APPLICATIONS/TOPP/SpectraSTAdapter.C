// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reirt
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
// $Maintainer: Matthias Seybold$
// $Authors: Matthias Seybold$
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
#include <stdio.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <QtGui>
#include <QApplication>

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
/*  virtual bool notify(QObject * receiver, QEvent * event) {
    try {
      return TOPPBase::QApplication::notify(receiver, event);
    } catch(std::exception& e) {
      qCritical() << "Exception thrown:" << e.what();
    }
    return false;
  }*/

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

    //int convert

    void registerOptionsAndFlags_()
    {
      registerInputFile_("spectrast_executable", "<executable>", "", "The 'spectrast' executable of the SpectraST installation", true, false, StringList::create("skipexists"));

      registerStringOption_("mode", "<type>", "", "SpectraST mode.", true, false);
      vector<String> mode_strings;
      mode_strings.push_back("search");
      mode_strings.push_back("create");
      setValidStrings_("mode", mode_strings);

      addEmptyLine_();

      registerTOPPSubsection_("create","Create Mode");
      addText_("General options");

      registerInputFileList_("create:in","<files>",StringList(),"two or more input files separated by blanks", false);
      setValidFormats_("create:in",StringList::create("mzML,mzXML,mzData,mgf,dta,msp,pepXML"));
      registerOutputFile_("create:out_splib","<file>", "", "splib Output files", false);
      registerOutputFile_("create:out_pepidx","<file>", "", "pepidx Output files", false);
      registerOutputFile_("create:out_sptxt","<file>", "", "sptxt Output files", false);
      registerOutputFile_("create:out_spidx","<file>", "", "spidx Output files", false);


      //option disabled because of problems with the output files
      //registerStringOption_("create:outputFileName", "<name>", "", "Output file base name", false);

      registerInputFile_("create:useProbTable", "<file>", "", "Use probability table in <file>. Only those peptide ions included in the table will be imported. A probability table is a text file with one peptide ion in the format AC[160]DEFGHIK/2 per line. If a probability is supplied following the peptide ion separated by a tab, it will be used to replace the original probability of that library entr", false);
      registerInputFile_("create:useProteinList", "<file>", "", "Use protein list in <file>. Only those peptide ions associated with proteins in the list will be imported.\n A protein list is a text file with one protein identifier per line. If a number X is supplied following the protein separated by a tab, then at most X peptide ions associated with that protein will be imported.", false);
      //registerStringOption_("create:remark", "<remark>", "", "Remark. Add a Remark=<remark> comment to all library entries created.", false);
      registerStringOption_("create:annotatePeaks","<type>", "true", "Annotate peaks.", false, false);
      registerStringOption_("create:binaryFormat","<type>", "true", "Write library in binary format, which enables quicker search.", false, false);
      registerStringOption_("create:writeDtaFiles","<type>", "false", "Write all library spectra as .dta files.", false, false);
      registerStringOption_("create:writeMgfFiles","<type>", "false", "Write all library spectra as one .mgf file.", false, false);

      addEmptyLine_();
      addText_("PEPXML IMPORT OPTIONS (Applicable with .pepXML file input) ");
      registerDoubleOption_("create:minimumProbabilityToInclude", "<prob>", 0.0, "Include all spectra identified with probability no less than <prob> in the library.", false);
      registerStringOption_("create:datasetName", "<name>", "", "Specify a dataset identifier for the file to be imported.", false);
      registerStringOption_("create:addMzXMLFileToDatasetName", "<type>", "false", "Add the originating mzXML file name to the dataset identifier. Good for keeping track of in which \nMS run the peptide is observed. (Turn off with -co!)", false, false);
      registerStringOption_("create:setDeamidatedNXST", "<type>", "false", "Set all asparagines (N) in the motif NX(S/T) as deamidated (N[115]). Use for glycocaptured peptides. (Turn off with -cg!).", false, false);
      registerStringOption_("create:setFragmentation", "<frag>", "","Set the fragmentation type of all spectra, overriding existing information.");
      registerStringOption_("create:setDeamidatedNXST", "<type>", "false", "Set all asparagines (N) in the motif NX(S/T) as deamidated (N[115]), and all asparagines not in the motif NX(S/T) as unmodified. Use for glycocaptured peptides.", false, false);
      //registerFlag_("cI", "Set the instrument and acquisition settings of the spectra.");

      addText_("LIBRARY MANIPULATION OPTIONS (Applicable with .splib files)");

      registerStringOption_("create:filterCriteria", "<pred>", "", "<pred> should be in quotes in the form “<attr> <op> <value>”. <attr> can refer to any of the fields and any comment entries. <op> can be ==, !=, <, >, <=, >=, =~ and !~. Multiple predicates can be separated by either & (AND logic) or | (OR logic), but not both. Default is off. ", false);
      registerStringOption_("create:combineAction", "<test>", "cJU", "Combine action.", false, false);
      setValidStrings_("create:combineAction", StringList::create("cJU,cJI,CJS,cJH"));
      registerStringOption_("create:buildAction", "<test>", "", "Build action.", false, false);
      setValidStrings_("create:buildAction", StringList::create("cAB,cAC,cAQ,cAD,cAM"));

      registerInputFile_("create:refreshDatabase","<file>", "" ,"Refresh protein mappings against the database <file> in FASTA format.", false, false);
      registerFlag_("create:refreshDeleteUnmapped", "Delete entries whose peptide sequences do not map to any protein during refreshing with -cD option.", true);
      registerFlag_("create:refreshDeleteMultimapped", "Delete entries whose peptide sequences map to multiple proteins during refreshing with the -cD option.", true);


      addEmptyLine_();
      addText_("CONSENSUS SPECTRUM CREATION OPTIONS (Applicable with -cAC option) ");
      registerIntOption_("create:minimumNumReplicates", "<num>", 1, "Minimum number of replicates required for each library entry. Peptide ions with fewer than <num> replicates will be excluded from library when creating consensus library.", false);
      registerStringOption_("create:removeDissimilarReplicates","<type>", "false", "Remove dissimilar replicates before creating consensus spectrum.", false, true);
      registerDoubleOption_("create:peakQuorum", "<frac>", 0.6, "Specify peak quorum: the fraction of all replicates required to contain a certain peak. Peaks not present in enough replicates will be deleted.", false, true);
      registerIntOption_("create:maximumNumPeaksUsed", "<num>", 300, "Maximum number of peaks in each replicate to be considered in creating consensus. Only the most intense <num> peaks by intensity will be considered.", false, true);
      registerIntOption_("create:maximumNumReplicates", "<num>", 100, "Maximum number of replicates used to build consensus spectrum.", false, true);
      registerIntOption_("create:maximumNumPeaksKept", "<num>", 150, "De-noise single spectra by keeping only the most intense <num> peaks.", false, true);
      registerStringOption_("create:replicateWeight", "<score>", "c_WGTS", "Select the type of score to weigh and rank the replicates.", false, true);
      setValidStrings_("create:buildAction", StringList::create("c_WGTS,c_WGTX,c_WGTP"));

      addEmptyLine_();
      addText_("BEST REPLICATE SELECTION OPTIONS (Applicable with -cAB option) ");
      //registerIntOption_("create:minimumNumReplicates", "<num>", 1, "Minimum number of replicates required for each library entry. Peptide ions with fewer than <num> replicates will be excluded from library when creating consensus library.", false);
      registerStringOption_("create:removeDissimilarReplicates","<type>", "false", "Remove dissimilar replicates before creating consensus spectrum.", false, true);

      addEmptyLine_();
      addText_("QUALITY FILTER OPTIONS (Applicable with -cAQ option) ");
      registerIntOption_("create:qualityLevelRemove","<level>",0,"Specify the stringency of the quality filter.",false);
      setMinInt_("create:qualityLevelRemove", 0);
      setMaxInt_("create:qualityLevelRemove", 5);
      registerIntOption_("create:qualityLevelMark","<level>",0,"-cL specifies the level for removal, -cl specifies the level for marking.\n<level> = 0: No filter.\n<level> = 1: Remove/mark impure spectra.\n<level> = 2: Also remove/mark spectra with a spectrally similar counterpart in the library that is better.\n<level> = 3: Also remove/mark inquorate entries (defined with -cr) that share no peptide sub-sequences with any other entries in the library.\n<level> = 4: Also remove/mark all singleton entries.\n<level> = 5: Also remove/mark all inquorate entries (defined with -cr).",false);
      setMinInt_("create:qualityLevelMark", 0);
      setMaxInt_("create:qualityLevelMark", 5);
      registerStringOption_("create:qualityPenalizeSingletons", "<type>", "true", "Apply stricter thresholds to singleton spectra during quality filters.", false, false);
      registerDoubleOption_("create:qualityImmuneProbThreshold", "<tresh>", 0.999, "Specify a probability above which library spectra are immune to quality filters.", false, true);
                            //0.999, "Specify a probability above which library spectra are immune to quality filters.", true);
      registerStringOption_("create:qualityImmuneMultipleEngines", "<type>", "true", "Make spectra identified by multiple sequence search engines immune to quality filters.", false, true);

      addEmptyLine_();
      addText_("DECOY GENERATION OPTIONS (Applicable with -cAD option) ");
      registerStringOption_("create:decoyConcatenate", "<type>", "false", "Concatenate real and decoy libraries.", false, false);
      registerIntOption_("create:decoySizeRatio", "<num>", 1, "Specify the (decoy / real) size ratio.", false, true);

      addEmptyLine_();
      addText_("SEMI-EMPIRICAL SPECTRUM GENERATION OPTIONS (Applicable with -cAM option) ");
      registerStringOption_("create:allowableModTokens", "<type>", "false", "Specify the set(s) of modifications allowed in semi-empirical spectrum generation by -cAM option.", false, true);



      addEmptyLine_();

      registerTOPPSubsection_("search","Search Mode");
      addText_("(II) Search Mode\nUsage: spectrast [ options ] <SearchFileName1> [ <SearchFileName2> ... <SearchFileNameN> ]\nwhere: SearchFileNameX = Name(s) of file containing unknown spectra to be searched.\n\t\t(Extension specifies format of file. Supports .mzXML, .mzData, .dta and .msp)");

      registerInputFile_("search:in_mzml", "<file>", "", "MZ Input file ");
      setValidFormats_("search:in_mzml",StringList::create("mzML,mzXML,mzData,mgf,dta,msp"));
      registerInputFile_("search:in_splib", "<file>", "", "splib Input file ");
      // Hannes fragen...
      //setValidFormats_("search:in_splib",StringList::create("splib"));
      registerInputFile_("search:in_spdix", "<file>", "", "spdix input file");
      registerOutputFile_("search:out", "<file>", "", "SpectraST search mode outputfile", true, false);
      registerInputFile_("search:libraryFile", "<file>", "", "Specify library file. Mandatory unless specified in parameter file. <file> must have .splib extension.", true, false);
      registerInputFile_("search:databaseFile", "<file>", "", "Specify a sequence database file. This will not affect the search in any way, but this information will be included in the output for any downstream data processing. <file> must have .fasta extension. If not set, SpectraST will try to determine this from the preamble of the library.", true, false);
      registerStringOption_("search:indexCacheAll", "<type>", "false", "Cache all entries in RAM. Requires a lot of memory (the library will usually be loaded almost in its entirety), but speeds up search for unsorted queries.", false, false);

      addEmptyLine_();
      addText_("CANDIDATE SELECTION AND SCORING OPTIONS ");
      registerDoubleOption_("search:indexRetrievalMzTolerance", "<tol>", 3.0, "Specify precursor m/z tolerance in Th. Monoisotopic mass is assumed.", false);
      registerStringOption_("search:indexRetrievalUseAverage", "<type>", "false", "Use average mass instead of monoisotopic mass.", false, false);

      vector<String> bool_strings;
      bool_strings.push_back("true");
      bool_strings.push_back("false");
      setValidStrings_("create:annotatePeaks", bool_strings);
      setValidStrings_("create:binaryFormat", bool_strings);
      setValidStrings_("create:writeDtaFiles", bool_strings);
      setValidStrings_("create:writeMgfFiles", bool_strings);
      setValidStrings_("create:addMzXMLFileToDatasetName", bool_strings);
      setValidStrings_("create:setDeamidatedNXST", bool_strings);
      setValidStrings_("create:setDeamidatedNXST", bool_strings);
      setValidStrings_("create:removeDissimilarReplicates", bool_strings);
      setValidStrings_("create:qualityPenalizeSingletons", bool_strings);
      setValidStrings_("create:qualityImmuneMultipleEngines", bool_strings);
      setValidStrings_("create:decoyConcatenate", bool_strings);
      setValidStrings_("create:allowableModTokens", bool_strings);
      setValidStrings_("search:indexCacheAll", bool_strings);
      setValidStrings_("search:indexRetrievalUseAverage", bool_strings);

    }

    ExitCodes main_(int , const char**)
    {
      // path to the log file
      String logfile(getStringOption_("log"));
      String spectrast_executable(getStringOption_("spectrast_executable"));
      //String unique_name = QDir::toNativeSeparators(String(File::getTempDirectory() + "/" + File::getUniqueName()).toQString()); // body for the tmp files

      //-------------------------------------------------------------
      // parsing parameters
      //-------------------------------------------------------------

      // get version of SpectraST
      QProcess qp;



      //qp.setWorkingDirectory("/tmp/test");
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

      //QStringList pepXMLFiles;
      QStringList qparam;
      QStringList arguments;

      QString filename;
      QString pepXMLfile;
      QString pepXMLbase;
      int qp_pid = getpid();
      QString tempdir = "/tmp/" + QString::number(qp_pid) + "/";

      QDir().mkdir(tempdir);

      // create parameters, todo: unterscheidung create und search

      if (getStringOption_("mode") == "search")
      {
        QString mzml_file = getStringOption_("search:in_mzml").toQString();
        QFileInfo fi_mzml(mzml_file);
        if (fi_mzml.exists())
        {
          QFile abc(mzml_file);
          abc.copy(tempdir + fi_mzml.fileName());
        }

        QString splib_file = getStringOption_("search:in_splib").toQString();
        QFileInfo fi_splib(splib_file);
        if (fi_splib.exists())
        {
          QFile abc(splib_file);
          abc.copy(tempdir + fi_splib.fileName());
        }

        QString spdix_file = getStringOption_("search:in_spdix").toQString();
        QFileInfo fi_spidx(spdix_file);
        if (fi_spidx.exists())
        {
          QFile abc(spdix_file);
          abc.copy(tempdir + fi_spidx.fileName());
        }
        qparam << "-sL" + tempdir + fi_splib.fileName();
        qparam << tempdir + fi_mzml.fileName();
      }



      else if (getStringOption_("mode") == "create")
      {
        StringList file_names = getStringList_("create:in");

        if (file_names.size() < 2)
        {

          writeLog_("SpectraST needs at least 1 pepXML and 1 mzXML file. Aborting!");
          printUsage_();
          return ILLEGAL_PARAMETERS;
        }

        for (StringList::Iterator file_it = file_names.begin();
             file_it != file_names.end(); ++file_it)
        {
          //cout <<  *file_it << '\n';
          QString copy = "cp";
          arguments << file_it->toQString();
          arguments << tempdir.toAscii().data();
          arguments.join("\t");

          filename = file_it->toQString();

          if (filename.endsWith("pepXML"))
          {
            pepXMLfile = file_it->toQString();
            QFileInfo fi(pepXMLfile);
            pepXMLbase = fi.fileName();
            pepXMLbase = tempdir + pepXMLbase; // + ".pepXML";
          }
          QProcess::execute(copy, arguments);
        }
        qparam << "-cP" + QString::number(getDoubleOption_("create:minimumProbabilityToInclude"));
        qparam << pepXMLbase;
      }

      cout << (String)qparam.join("\t") << endl;


      Int status = QProcess::execute(spectrast_executable.toQString(), qparam);

      if (status != 0)
      {
        writeLog_("Error: SpectraST problem! (Details can be seen in the logfile: \"" + logfile + "\")");
        writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory");

        if (getIntOption_("debug") <= 1)
        {
          writeDebug_("Removing temporary files", 10);
//          QString command = "rm -Rf " +  tempdir;

//          if (system(command.toAscii().data()))
//            // todo show warning.
//            writeDebug_("Cannot remove temporary files", 10);
        }
        return EXTERNAL_PROGRAM_ERROR;
      }

      // read SpectraST output
      writeDebug_("Reading output of SpectraST", 10);


      //copy spectrast create-mode output files
      if (getStringOption_("mode") == "create")
      {
        StringList tmp_files = QDir(tempdir).entryList(QDir::Files | QDir::Hidden);

        QString file_pepidx;
        file_pepidx = tempdir + File::removeExtension(File::basename(tmp_files[0])).toQString() + ".pepidx";
        QFileInfo fi_pepidx(file_pepidx);
        if (fi_pepidx.exists())
        {
          QFile abc(file_pepidx);
          abc.copy(getStringOption_("create:out_pepidx").toQString());
        }

        QString file_sptxt;
        file_sptxt = tempdir + File::removeExtension(File::basename(tmp_files[0])).toQString() + ".pepidx";
        QFileInfo fi_sptxt(file_sptxt);
        if (fi_sptxt.exists())
        {
          QFile abc(file_sptxt);
          abc.copy(getStringOption_("create:out_sptxt").toQString());
        }

        QString file_splib;
        file_splib = tempdir  + File::removeExtension(File::basename(tmp_files[0])).toQString() + ".splib";
        QFileInfo fi_splib(file_splib);
        if (fi_splib.exists())
        {
          QFile abc(file_splib);
          abc.copy(getStringOption_("create:out_splib").toQString()); // fileextension missing
        }

        QString file_spidx;
        file_spidx = tempdir  + File::removeExtension(File::basename(tmp_files[0])).toQString() + ".splib";
        QFileInfo fi_spidx(file_spidx);
        if (fi_spidx.exists())
        {
          QFile abc(file_spidx);
          abc.copy(getStringOption_("create:out_spidx").toQString()); // fileextension missing
        }

      }
      else if (getStringOption_("mode") == "search") // CRASH!!!
      {
        StringList tmp_files = QDir(tempdir).entryList(QDir::Files | QDir::Hidden);

        QString file_pepXML;
        file_pepXML = QDir(tempdir).absoluteFilePath((*(tmp_files.searchSuffix("pep.xml"))).toQString());
        QFileInfo fi_pepXML(file_pepXML);
        if (fi_pepXML.exists())
        {
          QFile abc(file_pepXML);
          abc.copy(getStringOption_("search:out").toQString()); // fileextension missing
        }

      }

      // TODO: delete temporary files
//      if (getIntOption_("debug") <= 1)
//      {
//        writeDebug_("Removing temporary files", 10);
//        QString command = "rm -Rf " +  tempdir;

//        if (system(command.toAscii().data()))
//          // todo show warning.
//          writeDebug_("Cannot remove temporary files", 10);
//      }

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPSpectraSTAdapter tool;

  return tool.main(argc,argv);
}

/// @endcond
