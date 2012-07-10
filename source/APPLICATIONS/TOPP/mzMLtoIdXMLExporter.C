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
// $Maintainer: $
// $Authors: Matthias Seybold $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <QFileInfo>
#include <QtCore/QFile>
#include <QStringList>


using namespace OpenMS;
using namespace std;

/**
  @page TOPP_mzMLtoIdXMLExporter mzMLtoIdXMLExporter

  @brief Imports an mzML file to an %OpenMS database.

  @deprecated Deprecated in OpenMS 1.9

  [Tool is no longer supported and will be removed in OpenMS 2.0]

  Besides the file to import, only the connection data has to be given.
  The data can then be retrieved by the @ref TOPP_DBExporter.

  The @em init flag can be used to create a new %OpenMS database.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_mzMLtoIdXMLExporter.cli
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES

class TOPPmzMLtoIdXMLExporter
  : public TOPPBase
{
  public:
    TOPPmzMLtoIdXMLExporter()
      : TOPPBase("mzMLtoIdXMLExporter","Imports data to an idXML File.")
    {

    }

  protected:
    void registerOptionsAndFlags_()
    {
      registerInputFile_("in", "<file>", "", "Textfile with mzML Filenames");
      registerStringOption_("folder", "<folder>", "", "Folder of the mzml files", true, false);
      registerOutputFile_("out", "<file>", "", "bla", true, false);
      registerOutputFile_("out_mzml", "<file>", "", "blablubb", true, false);

    }

    ExitCodes main_(int , const char**)
    {

      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------

      QString folder = getStringOption_("folder").toQString();



      //-------------------------------------------------------------
      // reading input
      //-------------------------------------------------------------

      TextFile input;

      vector<ProteinIdentification> protein_ids;
      vector<PeptideIdentification> peptide_ids;
      MSExperiment<> exp;
      MSExperiment<> out_exp;

      input.load(getStringOption_("in"));
      QString old_filepath = "";
      for (TextFile::Iterator it = input.begin()  ; it != input.end(); it++)
      {
        QString line = it->toQString();
        QStringList data = line.split("\t");
        //neu
        String seq_string = data.at(0).toAscii().data();
        AASequence seq = seq_string.substitute("m", "M(Oxidation)");
        QString filename = data.at(1);
        Peak2D::CoordinateType rt = data.at(2).toDouble();
        Peak2D::CoordinateType mz = data.at(3).toDouble();
        //neu
        int charge = data.at(4).toInt();

        QString filepath;
        filename = File::removeExtension(File::basename(filename)).toQString() + ".mzML";
        filepath = folder + filename;

        QFileInfo fi_filepath(filepath);
        if (fi_filepath.exists())
        {
          // rt und mz position des MS2 wird hier in einem vector von Peak2D übergeben

          vector<Peak2D> pcs;  // TODO: den musst du mit den daten rt und mz aus der Textdatei füllen
          Peak2D tmp;
          tmp.setRT(rt);
          tmp.setMZ(mz);
          pcs.push_back(tmp);

          MzMLFile mzml_file;
          if (filepath != old_filepath)
          {
            mzml_file.load(filepath, exp);
            old_filepath = filepath;
          }
          // aus Experiment exp wird extrahiert und in out_exp gespeichert
          for ( Size i = 0; i != exp.size(); ++i )
          {
            if ( exp[i].getMSLevel() == 2 )
            {
              if (!exp[i].getPrecursors().empty())
              {
                DoubleReal pc_mz = exp[i].getPrecursors()[0].getMZ();
                DoubleReal ms2_rt = exp[i].getRT(); // use rt of MS2 as we can't be sure there are MS1 in the experiment
                //bool found = false;
                for (Size j = 0; j != pcs.size(); ++j)
                {
                  DoubleReal mz_low = pcs[j].getMZ() - 0.1;
                  DoubleReal rt_low = pcs[j].getRT() - 0.1;
                  DoubleReal mz_high = pcs[j].getMZ() + 0.1;
                  DoubleReal rt_high = pcs[j].getRT() + 0.1;

                   cout << mz_low << " " << mz_high << " " << rt_low << " " << rt_high << endl;
                   cout << mz << " " << rt << endl;

                  if ( ms2_rt > rt_low && ms2_rt < rt_high && pc_mz > mz_low && pc_mz < mz_high)
                  {
                    cout << "DADADADADADADADADA" << endl;
                    out_exp.push_back(exp[i]);  // add spectrum

                    // TODO: eintrag zu peptide identifications hinzufügen
                    PeptideIdentification pep_id;
                    pep_id.setMetaValue("MZ", mz);
                    pep_id.setMetaValue("RT", rt);
                    pep_id.setScoreType("dot product"); // evtl anderen score eintragen
                    pep_id.setHigherScoreBetter(true);
                    PeptideHit pep_hit(0, 0, charge, seq);
                    //if at least one peptide hit is found
                    pep_id.insertHit(pep_hit);
                    peptide_ids.push_back(pep_id);


                    vector < ProteinIdentification > protein_ids(1);
                    protein_ids[0].setDateTime(DateTime::now());
                    protein_ids[0].setSearchEngine("MySearchEngine");
                    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

                    break;
                  }
                }
              }
            }
          }
        }

        else
        {
          writeLog_("ERROR: file " + filepath.toStdString() + " does not exist");
        }




      }

      //todo: out_exp in mzML speichern
      MzMLFile mzml_file;
      mzml_file.store(getStringOption_("out_mzml"), out_exp);


      IdXMLFile idxml_file;
      idxml_file.store(getStringOption_("out"), protein_ids, peptide_ids);

      return EXECUTION_OK;
    }
};

/// @endcond

int main( int argc, const char** argv )
{
  TOPPmzMLtoIdXMLExporter tool;
  return tool.main(argc,argv);
}
