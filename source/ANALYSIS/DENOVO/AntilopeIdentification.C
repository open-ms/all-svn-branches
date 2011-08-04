// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringBayes.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeAlgorithm.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePostScoring.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePMCorrect.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include<OpenMS/FORMAT/IdXMLFile.h>
#include<OpenMS/ANALYSIS/ID/IDMapper.h>


using namespace OpenMS;

class DeNovoSequencerLagrange : public TOPPBase
{
  public:
    DeNovoSequencerLagrange():TOPPBase("DeNovoSequencerLagrange", "De novo sequencing of the input spectrum using ...", false)
    {
    }

  protected:

    //typedefs
    typedef std::vector<DoubleReal>DVec;
    typedef std::vector<UInt>UIntVec;
    typedef std::vector<bool>BoolVec;
    typedef std::vector<String>StringVec;
    typedef SpectrumGraphSeqan::VertexDescriptor VertexDescriptor;

    void registerOptionsAndFlags_()
    {
      std::vector<String>tru_fal(2);
      tru_fal.push_back("true");
      tru_fal.push_back("false");

      registerInputFile_("mz_data_file","<file>","","input spectra file mzData format",true);
      registerInputFile_("scoring_file","<file>","MBtan.bin","scoring function parameter file",false);
      registerDoubleOption_("delta","e.g 0.5",0.5,"allowed mz error for every peak",false);
      registerDoubleOption_("precision","e.g 0.1",0.1,"mz precision",false);
      registerStringOption_("tryptic","true/false","true","Flag whether spectra come from trypsin digested peptides",false);
      setValidStrings_("tryptic",tru_fal);
    }


    ExitCodes main_(int, const char**)
    {
      //parameter handling
      String mz_data_file = getStringOption_("mz_data_file");
      String scoring_param_file = getStringOption_("scoring_file");
      DoubleReal delta = getDoubleOption_("delta");
      DoubleReal precision = getDoubleOption_("precision");
      String tryptic_flag = getStringOption_("tryptic");


      //checking validity of files

      //spectrum data type
      FileHandler fh;
      FileTypes::Type mz_type = fh.getTypeByFileName(mz_data_file);

      if(mz_type!=FileTypes::MZDATA)
      {
        writeLog_("wrong spectrum file format");
        return INCOMPATIBLE_INPUT_DATA;
      }

      //scoring parameter file
      if (!File::exists(scoring_param_file))
      {
        writeLog_("missing scoring parameter file");
        return INPUT_FILE_NOT_FOUND;
      }

      if (File::empty(scoring_param_file))
      {
        writeLog_("empty scoring parameter file");
        return INPUT_FILE_EMPTY;
      }


      //if all files exist and are non-empty continue
      PeakMap spec_map;
      MzDataFile().load(mz_data_file,spec_map);


      // ----Begin

      //for evaluation purposes   !!! TO BE REMOVED LATER !!
      //the idxml file containing the true annotations
      //String idxml_file_name="../data/spectra/ISB_spectra.idXml";
      String idxml_file_name="../../data/spectra/Benchmark_test2.idXML";
      //String idxml_file_name="../data/spectra/SoCe.idXml";

      std::vector<PeptideIdentification> pep_id_vec;
      std::vector<ProteinIdentification> prot_id_vec;
      String tmp_str;

      IdXMLFile().load(idxml_file_name,prot_id_vec, pep_id_vec,tmp_str);
      IDMapper().annotate(spec_map,pep_id_vec,prot_id_vec);

      // ----End

      //compute the necessary data structures
      IdSetup::BoolVec A;
      IdSetup::UIntVec A_len;
      IdSetup::AASeqVecMap A_map;

      //generateDataStructures(A, A_len, A_map, delta, precision);
      //IdSetup setup(delta,precision);
      IdSetup setup; //@TODO do not use default contructor. get remaining parameters from command line
      setup.create_vector_A(A, A_len, A_map);

      UInt counted_spectra=0;  
      for(PeakMap::iterator it = spec_map.begin(); it != spec_map.end(); ++it)
      {
        DoubleReal parent_mass = it->getPrecursors()[0].getPosition()[0];
        DoubleReal charge = it->getPrecursors()[0].getCharge();
        parent_mass = (parent_mass * charge) - (charge - 1);

        //check validity of precursor mass
        it->sortByPosition();
        if(parent_mass < it->back().getMZ())
        {
          //throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "incorrect precurser ion m/z");
          std::cerr<<"ignoring spectrum nr: "<<it-spec_map.begin()<<std::endl;
          //continue;
        }

        ++counted_spectra;

        //perform parent mass correction
        ParentMassCorrection p_corr;
        DoubleReal corrected_mass = p_corr(*it,2,0.1,0.1);

        AASequence true_seq_tmp = it->getPeptideIdentifications()[0].getHits()[0].getSequence();
        DoubleReal true_parent_mass = true_seq_tmp.getMonoWeight(Residue::Full)+1;
        std::cout << "true mass:  " << true_parent_mass << "  computed mass:  " << p_corr.getCorrectedMass() * charge - (charge - 1) << std::endl;


        //it->getPrecursors()[0].setPosition(corrected_mass);
        it->getPrecursors()[0].setPosition((true_parent_mass + charge - 1) / charge);

        const PeakSpectrum& input_spec = *it;

        //create the spectrum graph for the input MSSpectrum
        SpectrumGraphSeqan spec_graph;

        //parameter setting
        OpenMS::Param spec_params;
        spec_params.setValue("delta",delta);
        spec_params.setValue("precision",precision);
        spec_params.setValue("tryptic",tryptic_flag);
        spec_graph.setParameters(spec_params);


        //create node set
        spec_graph.createNodeSet(input_spec, A, A_len);

        //load the scoringfunction with a given parameter file
        BayesScoring scoring_func;

        scoring_func.setDelta(delta);
        scoring_func.loadFile(scoring_param_file);
        spec_graph.scoreNodes(scoring_func,input_spec, true, 1.0, -2.0);

        std::cout<<"NODE SUMMARY"<<std::endl;
        for(Size node = 1; node < spec_graph.size(); ++node)
        {
          std::cout<<"id: "<<node<<" mass: "<<spec_graph.getRealMass(node)<<"  score:  "<<spec_graph.getScore(node)<<std::endl;

          SpectrumGraphSeqan::InterpretVec interprets;
          spec_graph.getInterpretations(node, interprets);
          for(Size j=0; j<interprets.size(); ++j)
          {
            std::cout<<IdSetup::type_to_str((IdSetup::ion_type)interprets[j].type_id)<<"   "<<interprets[j].peak_nr<<std::endl;
          }
        }

        //create the set of directed and undirected edges
        spec_graph.createUdirEdgeSet(input_spec);
        spec_graph.createEdgeSet(A, A_len);

        //SOLVING THE SHORTEST PATH PROBLEM using Lagrangian relaxation
        YenAlgorithm yen(&spec_graph);

        clock_t tstart = clock();
        yen.computeLongestPaths(20);
        std::cout<<"True time: "<<(clock()-tstart)/(DoubleReal)CLOCKS_PER_SEC<<std::endl;

        cout<<"after compute_longest_path"<<endl;
        vector<vector<VertexDescriptor> >paths = yen.get_longest_paths();
        DVec result_scores= yen.get_path_scores();

        vector<UIntVec>result_masses(paths.size());

        for(Size i = 0; i<paths.size();++i)
        {
          for(Size j = 0; j<paths[i].size()-1; ++j)
          {
            result_masses[i].push_back(spec_graph.getIntEdgeMass(paths[i][j], paths[i][j+1]));
          }
        }

        AASequence true_seq=it->getPeptideIdentifications()[0].getHits()[0].getSequence();

        //PROCESSING THE RESULTS
        std::cout<<"TRUE SEQUENCE:  "<<true_seq<<std::endl;
        for (UInt i = 0; i < result_masses.size(); ++i)
        {
          IdEval::IdAnnot tmp_annot;
          DoubleReal prefix=0;

          IdEval::generateSampleAnnotation(A_map, precision, result_masses[i], tmp_annot);
          for(UInt k=0; k<result_masses[i].size();++k)
          {
            prefix+=result_masses[i][k];
            std::cout<<prefix<<" -- ";
          }
          std::cout<<std::endl;

          std::cout << i << ":    ";
          tmp_annot.print();
          std::cout << std::endl;
          std::cout<<"SCORE: "<<result_scores[i]<<std::endl;
          std::cout<<"PATH:"<<std::endl;
          for(Size l=0; l<paths[i].size(); ++l)
          {
            std::cout<<paths[i][l]<<std::endl;
          }

          UInt correct_aa=0;
          DoubleReal correct_percents=0.0;
          IdEval::predictionQuality(result_masses[i], precision, it->getPeptideIdentifications()[0].getHits()[0].getSequence() ,correct_aa,correct_percents);
        }
      }
      return EXECUTION_OK;
    }
};


//-----------------------------------------------------------------------------------------------
//------------------------start of main----------------------------------------------------------
//-----------------------------------------------------------------------------------------------
int main(int argc, const char** argv)
{
  DeNovoSequencerLagrange tool;
  return tool.main(argc,argv);
}
