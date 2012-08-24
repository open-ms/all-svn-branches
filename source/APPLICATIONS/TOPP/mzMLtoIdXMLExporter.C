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
#include <OpenMS/KERNEL/StandardTypes.h>



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


struct IDInfo
{
  DoubleReal mz;
  DoubleReal rt; // seconds
  Int charge;
  AASequence seq;
  PeptideIdentification pep_id;
  Size scan_index; // scan index in out_exp

};

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
    registerOutputFile_("out_idxml", "<file>", "", "idXML output file", true, false);
    registerOutputFile_("out_mzml", "<file>", "", "mzML output file", true, false);
    registerOutputFile_("out_mzml_test", "<file>", "", "mzML test file", true, false);
    registerOutputFile_("out_idxml_test", "<file>", "", "idXML output file", true, false);
    registerOutputFile_("out_mzml_train", "<file>", "", "mzML train file", true, false);
    registerOutputFile_("out_idxml_train", "<file>", "", "idXML output file", true, false);

    registerIntOption_("threshold", "<num>", 50, "Percentage of peptide hits used for Test Data", false, false);

  }

  ExitCodes main_(int , const char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    QString folder = getStringOption_("folder").toQString();
    float threshold = getIntOption_("threshold");


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    TextFile input;

    vector<PeptideIdentification> peptide_ids;
    vector<PeptideIdentification> out_test_pep_ids;
    vector<PeptideIdentification> out_train_pep_ids;
    MSExperiment<> exp;
    MSExperiment<> out_exp = exp;
    out_exp.clear(false);
    MSExperiment<> out_test = exp;
    out_test.clear(false);
    MSExperiment<> out_train = exp;
    out_train.clear(false);

    input.load(getStringOption_("in"));
    QString old_filepath = "";
    map<String, vector<IDInfo> > map_seq2id;

    for (TextFile::Iterator it = input.begin()  ; it != input.end(); it++)
    {
      QString line = it->toQString();
      QStringList data = line.split("\t");
      //neu

      IDInfo id;
      String seq_string = data.at(0).toAscii().data();
      id.seq = seq_string.substitute("s", "S(Phospho)").substitute("t", "T(Phospho)").substitute("y", "Y(Phospho)").substitute("e", "E").substitute("d", "D")
          .substitute("k", "K").substitute("i", "I").substitute("l", "I").substitute("m", "M(Oxidation)").substitute("q", "P(Pro->pyro-Glu)"); //substitute("k","K(Carbamyl)").
      QString filename = data.at(1);
      id.rt = data.at(2).toDouble() * 60.0;   //to get rt in seconds
      id.mz = data.at(3).toDouble();
      id.charge = data.at(4).toInt();

      QString filepath;
      filename = File::removeExtension(File::basename(filename)).toQString() + ".mzML";
      filepath = folder + filename;

      QFileInfo fi_filepath(filepath);
      if (fi_filepath.exists())
      {
        // rt und mz position des MS2 wird hier in einem vector von Peak2D übergeben
        vector<Peak2D> pcs;
        Peak2D tmp;
        tmp.setRT(id.rt);
        tmp.setMZ(id.mz);
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
          if ( exp[i].getMSLevel() == 2 ) // is tandem spectrum
          {
            if (!exp[i].getPrecursors().empty())
            {
              DoubleReal pc_mz = exp[i].getPrecursors()[0].getMZ();


              DoubleReal ms2_rt_s = exp[i].getRT(); // use rt of MS2 as we can't be sure there are MS1 in the experiment
              //bool found = false;
              for (Size j = 0; j != pcs.size(); ++j)
              {
                DoubleReal mz_low = pcs[j].getMZ() - 0.1;
                DoubleReal rt_low = pcs[j].getRT() - 1.0; // our data is not very accurate on RT
                DoubleReal mz_high = pcs[j].getMZ() + 0.1;
                DoubleReal rt_high = pcs[j].getRT() + 1.0;

                if ( ms2_rt_s > rt_low && ms2_rt_s < rt_high && pc_mz > mz_low && pc_mz < mz_high) // todo: compare precursor charge to id charge
                {
                  //                    cout << pc_mz << " " << ms2_rt_s << endl;
                  //                    cout << mz_low << " " << mz_high << " " << rt_low << " " << rt_high << endl;

                  PeakSpectrum ps = exp[i];
                  String tmp = ps.getNativeID();
                  String tmp2;
                  tmp2 = tmp.substr( 0, 41);

                  id.scan_index = out_exp.size();
                  String tmp3 = tmp2 + String( id.scan_index + 1 );
                  ps.setNativeID( tmp3 );

                  out_exp.push_back(ps);  // add spectrum
                  DoubleReal score = 0;
                  uInt rank = 0;

                  PeptideIdentification pep_id;
                  pep_id.setMetaValue("MZ", pc_mz);
                  pep_id.setMetaValue("RT", ms2_rt_s);
                  pep_id.setScoreType("Mascot"); // evtl anderen score eintragen
                  pep_id.setHigherScoreBetter(true);
                  PeptideHit pep_hit(score, rank, id.charge, id.seq);
                  //if at least one peptide hit is found
                  pep_id.insertHit(pep_hit);
                  peptide_ids.push_back(pep_id);
                  id.pep_id = pep_id;
                  map_seq2id[id.seq.toString()].push_back(id);
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


    for (map<String, vector<IDInfo> >::const_iterator it = map_seq2id.begin(); it != map_seq2id.end(); ++it)
    {

      UInt n_test = (UInt)((DoubleReal)threshold / 100.0 * (DoubleReal)it->second.size());

//      cout << "seq: " << it->first << " entries: " << it->second.size() << " n_test:" << n_test << endl;
      for (Size j = 0; j != n_test; ++j)
      {
//        cout << "train set: " << it->second[j].scan_index << " " << out_exp.size() << endl;
        out_train.push_back(out_exp[it->second[j].scan_index]);
        out_train_pep_ids.push_back(it->second[j].pep_id);
      }
      for (Size j = n_test; j < it->second.size(); ++j)
      {
//        cout << "test set: " << it->second[j].scan_index << " " << out_exp.size() << endl;
        out_test.push_back(out_exp[it->second[j].scan_index]);
        // TIMO out_test_pep_ids.push_back(it->second[j].pep_id);
        out_test_pep_ids.push_back(it->second[j].pep_id);

      }
    }






    MSExperiment<> new_out_train;
    Size index = 1;
    for ( Size i = 0; i != out_train.size(); ++i )
    {
      PeakSpectrum ps = out_train[i];
      String tmp = ps.getNativeID();
      String tmp2;
      tmp2 = tmp.substr( 0, 41);
      String tmp3 = tmp2 + String( index );
      ps.setNativeID( tmp3 );
      new_out_train.push_back(ps);
      index++;

    }


    MSExperiment<> new_out_test;
    index = 1;
    for ( Size i = 0; i != out_test.size(); ++i )
    {
      PeakSpectrum ps = out_test[i];
      String tmp = ps.getNativeID();
      String tmp2;
      tmp2 = tmp.substr( 0, 41);
      String tmp3 = tmp2 + String( index );
      ps.setNativeID( tmp3 );
      new_out_test.push_back(ps);
      index++;

    }







    cout << "size (total/train/test) sets: " << out_exp.size() << " " << out_train.size() << " " << out_test.size() << endl;

    MzMLFile mzmlfile_test;
    mzmlfile_test.store(getStringOption_("out_mzml_test"), new_out_test);

    MzMLFile mzmlfile_train;
    mzmlfile_train.store(getStringOption_("out_mzml_train"), new_out_train);

    MzMLFile mzmlfile;
    mzmlfile.store(getStringOption_("out_mzml"), out_exp);

    vector < ProteinIdentification > protein_ids(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("MySearchEngine");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    IdXMLFile idxml_file;
    idxml_file.store(getStringOption_("out_idxml"), protein_ids, peptide_ids);

    IdXMLFile idxml_file_test;
    idxml_file_test.store(getStringOption_("out_idxml_test"), protein_ids, out_test_pep_ids);

    IdXMLFile idxml_file_train;
    idxml_file_train.store(getStringOption_("out_idxml_train"), protein_ids, out_train_pep_ids);

    return EXECUTION_OK;
  }
};

/// @endcond

int main( int argc, const char** argv )
{
  TOPPmzMLtoIdXMLExporter tool;
  return tool.main(argc,argv);
}
