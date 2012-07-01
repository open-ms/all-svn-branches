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
// $Maintainer: Matthias Seybold $
// $Authors: Matthias Seybold$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
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
@page TOPP_PepitomeAdapter PepitomeAdapter

@brief Identifies peptides in MS/MS spectra via Pepitome (Open Mass Spectrometry Search Algorithm).

<CENTER>
<table>
 <tr>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ PepitomeAdapter \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
 </tr>
 <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
 </tr>
</table>
</CENTER>

@em Pepitome must be installed on the system to be able to use the @em PepitomeAdapter. See the Pepitome site
for further information on how to download and install @em Pepitome on your system. You might find that the lacurrdir Pepitome version
does not run on your system (to currdir this, run @em Pepitome from command line). If you encounter
an error message, try another Pepitome version

This adapter supports relative database filenames, which (when not found in the current working directory) is looked up in
the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

This wrapper has been currdired successfully with Pepitome, version 1.x.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PepitomeAdapter.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPPepitomeAdapter
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
 TOPPPepitomeAdapter()
   : TOPPBase("PepitomeAdapter","Annotates MS/MS spectra using Pepitome.")
 {
 }

protected:

 struct PepitomeVersion
 {
   PepitomeVersion ()
     : Pepitome_major(0), Pepitome_minor(0)
   {}

   PepitomeVersion (Int maj, Int min)
     : Pepitome_major(maj), Pepitome_minor(min)
   {}

   Int Pepitome_major;
   Int Pepitome_minor;

   bool operator < (const PepitomeVersion& v) const
   {
     if (Pepitome_major > v.Pepitome_major) return false;
     else if (Pepitome_major < v.Pepitome_major) return true;
     else // ==
     {
       if (Pepitome_minor > v.Pepitome_minor) return false;
       else if (Pepitome_minor < v.Pepitome_minor) return true;
     }
     return false;
   }
 };

 bool getVersion_(const String& version, PepitomeVersion& Pepitome_version_i) const
 {
   // we expect two components separated by '('
   IntList nums = IntList::create(StringList::create(version,'.'));
   if (nums.size() != 2)
   {
     return false;
   }
   Pepitome_version_i.Pepitome_major = nums[0];
   Pepitome_version_i.Pepitome_minor = nums[1];
   return true;
 }

 //int convert

 void registerOptionsAndFlags_()
 {
   registerInputFile_("Pepitome_executable", "<executable>", "", "The 'Pepitome' executable of the Pepitome installation", true, false, StringList::create("skipexists"));
   registerInputFile_("SpectralLibrary", "<library>", "", "SPTXT formated spectral library", true, false);
   registerInputFile_("ProteinDatabase", "<fasta>","", "FASTA protein database", true, false);
   registerInputFile_("mzML", "<file>", "", "MS/MS data filepath in a supported file format", true, false);
   registerOutputFile_("out", "<file>", "", "idXML Output File", true, false);








 }

 ExitCodes main_(int , const char**)
 {
   // path to the log file
   String logfile(getStringOption_("log"));
   String Pepitome_executable(getStringOption_("Pepitome_executable"));
   String unique_name = QDir::toNativeSeparators(String(File::getTempDirectory() + "/" + File::getUniqueName()).toQString()); // body for the tmp files

   //-------------------------------------------------------------
   // parsing parameters
   //-------------------------------------------------------------

   // get version of Pepitome
   QProcess qp;
   qp.setProcessChannelMode(QProcess::MergedChannels);
   qp.start(Pepitome_executable.toQString());
   bool success = qp.waitForFinished(-1);
   String output(QString(qp.readAllStandardOutput()));
   String Pepitome_version;
   PepitomeVersion Pepitome_version_i;
   if (!success/* || qp.exitStatus() != 0 || qp.exitCode() != 0*/)
   {
     writeLog_("Warning: unable to determine the version of Pepitome - the process returned an error. Call string was: " + Pepitome_executable + ". Make sure that the path to the Pepitome executable is correct!");
     return ILLEGAL_PARAMETERS;
   }
   else
   {
     vector<String> version_split;
     output.split("Pepitome ", version_split);
     String tmp = version_split[1];
     tmp.split(" (", version_split);
     String version_line = version_split[0];

     //cout << version_line << endl;


     if (getVersion_(version_line, Pepitome_version_i))
     {
       Pepitome_version = version_line;
       writeDebug_("Setting Pepitome version to " + Pepitome_version, 1);
     }
     else
     {
       writeLog_("Warning: Pepitome version output (" + output + ") not formatted as expected!");
     }
   }

   QStringList qparam;

   qparam << "-SpectralLibrary" << getStringOption_("SpectralLibrary").toQString();
   qparam << "-ProteinDatabase" << getStringOption_("ProteinDatabase").toQString();
   qparam << getStringOption_("mzML").toQString();

   //TODO: parse remaining command line options

   cout << (String)qparam.join("\t") << endl;

   Int status = QProcess::execute(Pepitome_executable.toQString(), qparam);
   if (status != 0)
   {
     writeLog_("Error: Pepitome problem! (Details can be seen in the logfile: \"" + logfile + "\")");
     writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory");

     if (getIntOption_("debug") <= 1)
     {
       // TODO: if no debuggin is enabled, delete temporary Pepitome param file
       writeDebug_("Removing temporary files", 10);
       // QFile(bla.toQString()).remove();
     }
     return EXTERNAL_PROGRAM_ERROR;
   }

   // read Pepitome output
   writeDebug_("Reading output of Pepitome", 10);

   // TODO: delete temporary files
   if (getIntOption_("debug") <= 1)
   {
     writeDebug_("Removing temporary files", 10);
     // QFile(bla.toQString()).remove();
   }

   //-------------------------------------------------------------
   // store output
   //-------------------------------------------------------------

   vector<ProteinIdentification> protein_ids;
   vector<PeptideIdentification> peptide_ids;
   IdXMLFile id_xml;

   QDir currdir;
   QString tmp;
   QString file_pepxml;

   tmp = currdir.currentPath() + "/";
   file_pepxml = tmp + File::removeExtension(File::basename(getStringOption_("mzML"))).toQString() + ".pepXML";

   QFileInfo fi_pepxml(file_pepxml);
   if (fi_pepxml.exists())
   {
     PepXMLFile pep_xml;
     pep_xml.load(file_pepxml.toStdString(), protein_ids, peptide_ids);
     id_xml.store(getStringOption_("out"), protein_ids, peptide_ids);
   }
   else
   {
     writeLog_("Error: Pepitome problem! " + file_pepxml + "not found!");
   }



   //pep_xml.load(String(pepxml_filename.toStdString()), protein_ids, peptide_ids);
   //id_xml.store(getStringOption_("out"), protein_ids, peptide_ids);



   return EXECUTION_OK;
 }
};


int main( int argc, const char** argv )
{
TOPPPepitomeAdapter tool;

return tool.main(argc,argv);
}

/// @endcond
