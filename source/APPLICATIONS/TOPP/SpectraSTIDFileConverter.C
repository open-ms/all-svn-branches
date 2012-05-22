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
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PepXMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**


*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPSpectraSTIDFileConverter : public TOPPBase
{
public:
  TOPPSpectraSTIDFileConverter() :
    TOPPBase("SpectraSTIDFileConverter", "Converts identification engine file formats.", true)
  {
  }

protected:
  void
  registerOptionsAndFlags_()
  {
    registerInputFile_("in_idXML", "<file>", "", "bla", true, false);
    registerInputFile_("in_mzml", "<file>", "", "bla", true);
  //  registerInputFile_("in_mzxml", "<file>", "", "bla", true);
    registerOutputFile_("out_pepXML", "<file>", "", "Output IdXML", true);
    //registerOutputFile_("out_mzXML", "<file>", "", "Out mzXML", false);
    //String formats("pepXML");
    //setValidFormats_("out_pepXML", StringList::create(formats));
  }

  ExitCodes
  main_(int, const char**)
  {
    String in_mzml = getStringOption_("in_mzml");
    String in_idXML = getStringOption_("in_idXML");
    String out_pepXML = getStringOption_("out_pepXML");
    //String out_mzXML = getStringOption_("out_mzXML");


    // load experiment data
    MSExperiment<> exp;
    MzMLFile mzml_file;
    mzml_file.load( in_mzml, exp );

    // load identification results
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids;
    IdXMLFile idxml_file;
    idxml_file.load(in_idXML, protein_ids, peptide_ids);

    for ( Size i = 0; i != exp.size(); ++i )
    {
      if ( exp[i].getMSLevel() == 2 ) // is tandem spectrum?
      {
        if (!exp[i].getPrecursors().empty())
        {
          DoubleReal ms2_rt = exp[i].getRT();

          for (vector<PeptideIdentification>::iterator cit = peptide_ids.begin(); cit != peptide_ids.end(); ++cit)
          {
            DoubleReal id_rt = cit->getMetaValue("RT");
            DoubleReal rt_low = id_rt - 0.00001;
            DoubleReal rt_high = id_rt + 0.00001;

            if ( ms2_rt > rt_low && ms2_rt < rt_high)
            {
              cit->setMetaValue("RT_index", i);
            }
          }
        }
      }
    }

    //idxml_file.store(out, protein_ids, peptide_ids);
    PepXMLFile pep_file;
    pep_file.store(out_pepXML, protein_ids, peptide_ids);

  //  String in_mzxml = getStringOption_("in_mzxml");
  //  MzXMLFile mzxml_file;
  //  mzxml_file.load( in_mzxml, exp );
  //  String out_mzXML = out_pepXML.substitute("pepXML", "mzXML");
  //  mzxml_file.store( out_mzXML, exp );


  return EXECUTION_OK;
}
};

int main(int argc, const char** argv)
{
  TOPPSpectraSTIDFileConverter tool;
  return tool.main(argc, argv);
}

///@endcond
