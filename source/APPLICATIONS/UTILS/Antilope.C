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

#define LAGRANGE
#undef USECPLEX

//#define MATRIX
#define ONLYPATHS

#undef Debug
//#define Debug

//#include <ilcplex/ilocplex.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
//#include <Scoring/ScoringFunctionMB.h>
//#include <Scoring/BayesScoringSmile.h>
//#include <Scoring/BayesScoring.h>
//#include <Scoring/ScoringFunctionPepNovo.h>
//#include <Graph/SpectrumGraph.h>
//#include <Graph/SpectrumGraph_boost.h>
//#include <Solver/ILP/ILPSolver.h>
//#include <Solver/ILP/ILPSolverCoinOr.h>
//#include <Preprocessing/IdSetup.h>


#ifdef USECPLEX
#include <ilcplex/ilocplex.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeILP.h>
#endif

#include <OpenMS/ANALYSIS/DENOVO/AntilopeAlgorithm.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeSpectrumGraph.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePMCorrect.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePostScoring.h>




//#include <Preprocessing/ParentMassCorrection.h>
//#include <Preprocessing/SpectrumPreprocessorBase.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>

#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>

#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>

#include <QtCore/QObject>
#include <QtCore/QProcess>
#include <QtCore/QByteArray>

//#include <smile/smile.h>
//#include <smile/smilearn.h>


using namespace OpenMS;

class DeNovoSequencer : public TOPPBase
{
  public:
    DeNovoSequencer():TOPPBase("DeNovoSequencer", "De novo sequencing of the input spectrum using ...", false)
    {
    }

  protected:

    //typedefs
    typedef std::vector<UInt>UIntVec;
    typedef std::vector<bool>BoolVec;
    typedef std::vector<String>StringVec;

    void registerOptionsAndFlags_()
    {
      std::vector<String>tru_fal;
      tru_fal.push_back("true");
      tru_fal.push_back("false");

      registerInputFile_("mz_data_file","<file>","","input spectra file mzData format",true);
      registerInputFile_("scoring_file","<file>","MBtan.bin","scoring function parameter file",false);
      //registerInputFile_("rank_scoring_file","<file>","rankscores","scoring function parameter file",false);
      registerDoubleOption_("delta","e.g 0.5",0.5,"allowed mz error for every peak",false);
      registerDoubleOption_("precision","e.g 0.1",0.1,"mz precision",false);      
      registerStringOption_("tryptic","true/false","true","Flag whether spectra come from trypsin digested peptides",false);
      setValidStrings_("tryptic",tru_fal);
      registerFlag_("local_hits", "set to 1 if local hits shall be considered", 0);      
    }

    ExitCodes main_(int, const char**)
    {


     AASequence seq1("RAGGLKL");
     AASequence seq2("RA[114.1]LKL");
     DoubleReal corr_percc;
     UInt corr_aaa;
     IdEval::predictionQuality(IdEval::IdAnnot(seq1), IdEval::IdAnnot(seq2), corr_aaa, corr_percc);
     std::cerr<<"CHECK QUALITY  "<<corr_aaa<<"::"<<corr_percc<<"  "<<seq2.toUnmodifiedString().size()<<std::endl;

      //parameter handling
      String mz_data_file = getStringOption_("mz_data_file");
      String scoring_param_file = getStringOption_("scoring_file");
      //String rankscoring_file = getStringOption_("rank_scoring_file");
      DoubleReal delta = getDoubleOption_("delta");
      DoubleReal precision = getDoubleOption_("precision");
      String tryptic_flag = getStringOption_("tryptic");
      bool local_hits = getFlag_("local_hits");
      local_hits=false;


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

//      //scoring parameter file
//      if (!File::exists(rankscoring_file))
//      {
//        writeLog_("missing scoring parameter file");
//        return INPUT_FILE_NOT_FOUND;
//      }

//      if (File::empty(rankscoring_file))
//      {
//        writeLog_("empty scoring parameter file");
//        return INPUT_FILE_EMPTY;
//      }


      //if all files exist and are non-empty continue
      PeakMap spec_map;
      MzDataFile().load(mz_data_file,spec_map);


      // ----Begin

      //for evaluation purposes   !!! TO BE REMOVED LATER !!
      //the idxml file containing the true annotations
      //String idxml_file_name="../data/spectra/ISB_spectra.idXml";
      String idxml_file_name="./data/Benchmark_test2.idXML";
      //String idxml_file_name="ISB_short_100.idXml";

      //String idxml_file_name="../data/spectra/SoCe.idXml";

      std::vector<PeptideIdentification> pep_id_vec;
      std::vector<ProteinIdentification> prot_id_vec;
      String tmp_str;

      OpenMS::Param mapper_param;
      mapper_param.setValue("rt_tolerance",0.001);
      mapper_param.setValue("mz_tolerance",0.001);
      mapper_param.setValue("mz_measure","Da");
      IDMapper mapper;
      mapper.setParameters(mapper_param);
      IdXMLFile().load(idxml_file_name,prot_id_vec, pep_id_vec,tmp_str);
      mapper.annotate(spec_map,pep_id_vec,prot_id_vec);

//      //DEBUG:
//      for(Size s=0; s<spec_map.size();++s)
//      {
//        std::cout<<"SeqComp: "<<pep_id_vec[s].getHits()[0].getSequence()<<" :::  "<<spec_map[s].getPeptideIdentifications()[0].getHits()[0].getSequence()<<std::endl;
//      }

      // ----End

      //compute the necessary data structures
      IdSetup::BoolVec A;
      IdSetup::UIntVec A_len;
      IdSetup::AASeqVecMap A_map;

      IdSetup setup; //TODO do not use default contructor. get remaining parameters from command line
      setup.create_vector_A(A, A_len, A_map);     

      UInt total_correct_residues = 0, total_possible_correct_residues=0, total_residues =0, total_quick_residues=0, total_residues_precision=0;
      UInt counted_spectra=0;

      //BayesScoringSmile scoring_func;
      BayesScoring scoring_func;
      //OpenMS::Param scr_params;
      //scr_params.setValue("delta",delta);
      //scoring_func.setParameters(scr_params);
      scoring_func.loadFile(scoring_param_file);

      //Scored types required by psm_score function in IdEval
      std::set<IdSetup::ion_type>scored_types;
      scored_types.insert(IdSetup::BIon);
      scored_types.insert(IdSetup::YIon);
      scored_types.insert(IdSetup::BIon2);
      scored_types.insert(IdSetup::YIon2);
      scored_types.insert(IdSetup::AIon);
      scored_types.insert(IdSetup::BIon_h2o);
      scored_types.insert(IdSetup::BIon_nh3);
      scored_types.insert(IdSetup::YIon_h2o);
      scored_types.insert(IdSetup::YIon_nh3);

      DoubleReal time_solver_summed = 0;
      DoubleReal time_solver_netto_summed = 0;
      for(PeakMap::iterator it = spec_map.begin(); it != spec_map.end(); ++it)
      //for(PeakMap::iterator it=spec_map.begin(); it!=spec_map.begin()+2;++it)
      {
        //std::cout<<"Iteration: "<<it-spec_map.begin()<<"  Sequence: "<<it->getPeptideIdentifications()[0].getHits()[0].getSequence()<<std::endl;
        std::cout << "Iteration: " << it-spec_map.begin() << std::endl;
        DoubleReal parent_mass = it->getPrecursors()[0].getPosition()[0];
        DoubleReal charge = it->getPrecursors()[0].getCharge();
        parent_mass = (parent_mass * charge) - (charge-1);

//        if (parent_mass>1400)
//        {
//          std::cout<<"dropping sequence: "<<it->getPeptideIdentifications()[0].getHits()[0].getSequence()<<" having weight: "<<parent_mass<<std::endl;
//          continue;
//        }
//        Size tmp_size=it->getPeptideIdentifications()[0].getHits()[0].getSequence().size();
//        if(it->getPeptideIdentifications()[0].getHits()[0].getSequence()[tmp_size-1]!='K' && it->getPeptideIdentifications()[0].getHits()[0].getSequence()[tmp_size-1]!='R')
//        {
//          std::cout<<"dropping sequence: "<<it->getPeptideIdentifications()[0].getHits()[0].getSequence()<<" due to non tryptic: "<<std::endl;
//          continue;
//        }
//        if (counted_spectra>250)
//        {
//            break;
//        }

        //check validity of precursor mass
        it->sortByPosition();

        //perform parent mass correction
        ParentMassCorrection p_corr;
        DoubleReal corrected_mass = p_corr(*it,2,0.1,0.1);

        AASequence true_seq_tmp=it->getPeptideIdentifications()[0].getHits()[0].getSequence();
        DoubleReal true_parent_mass=true_seq_tmp.getMonoWeight(Residue::Full, 1);
        std::cout<<"true mass:  "<<true_parent_mass<<"  computed mass:  "<<corrected_mass*charge -(charge-1)<<"  original mass: "<<parent_mass<<std::endl;

        it->getPrecursors()[0].setPosition(corrected_mass);
//        //it->getPrecursors()[0].setPosition((true_parent_mass+charge-1)/charge);
//        if(fabs(((true_parent_mass+charge-1)/charge)-corrected_mass)>=0.5)
//        {
//          std::cout<<"dropping sequence: "<<it->getPeptideIdentifications()[0].getHits()[0].getSequence()<<" having high error"<<std::endl;
//          continue;
//        }

        DoubleReal time_solver = clock();
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
        //spec_graph.scoreNodes(scoring_func,input_spec, true, 0.5, -0.5);
#ifdef Debug
        std::cout<<"NODE SUMMARY"<<std::endl;
        for(UInt node=0; node<spec_graph.size(); ++node)
        {
         std::cout << "id: " << node << " mass: " << spec_graph.getRealMass(node) << "  score:  " << spec_graph.getScore(node)<<std::endl;
        }
#endif

        //create the set of directed and undirected edges
        spec_graph.createUdirEdgeSet();
        spec_graph.createEdgeSet(A, A_len);
        spec_graph.setEdgeWeights();
        //spec_graph.createSpanningEdges();
        
#ifdef MATRIX
        clock_t mat_start = clock();
        MatrixSpectrumGraph mat_graph;
        mat_graph.createMatrixGraph2(spec_graph);
        
//        MatrixSpectrumGraph mat_tmp;
//        mat_tmp.createMatrixGraph2(spec_graph);
#endif

        //TODO: DO this during construction of spec graph!
/*
        //make local solutions possible by creating more edges
        DoubleReal upper_mass_bound = spec_graph.getRealMass(spec_graph.size()-1) - 400.0;
        if (local_hits)
        {
          for(UInt node=0; node<spec_graph.size(); ++node)
          {
            DoubleReal real_mass = spec_graph.getRealMass(node);
            if(!spec_graph.findDirectedEdge(0,node) && real_mass>400)
            {
              spec_graph.addDirEdgeToGraph(0,node,1);
            }
            if(!spec_graph.hasDirectedEdge(node,spec_graph.size()-1) && real_mass<upper_mass_bound)
            {
              spec_graph.addDirEdgeToGraph(node,spec_graph.size()-1,1);
            }
          }
        }

*/
        //SOLVING THE SHORTEST PATH PROBLEM

        //use the solver to compute the longest antisymmetric paths
        //de_novo_ILP ilp;

        std::vector<std::vector<UInt> > result_masses;
        std::vector<DoubleReal> scores;

        DoubleReal time_solver_netto = clock();
#ifdef LAGRANGE

#ifdef MATRIX
        YenAlgorithm solver(&mat_graph);
#else
        YenAlgorithm solver(&spec_graph);
#endif
        
        vector<vector<SpectrumGraphSeqan::VertexDescriptor> > paths;
        
        solver.computeLongestPaths(50, paths, scores);
//        std::cout << "TOTAL Matrix TIME: "<< (clock()-mat_start)/(DoubleReal)CLOCKS_PER_SEC << std::endl;
        //scores = solver.get_path_scores();
        //std::vector< std::vector<SpectrumGraphSeqan::VertexDescriptor> > paths = solver.get_longest_paths();

        result_masses.clear();
#ifdef MATRIX
        mat_graph.getMassesForPaths(spec_graph, paths, result_masses);
#else
        spec_graph.getMassesForPaths(paths, result_masses);
#endif

#endif
#ifdef USECPLEX

        de_novo_ILP ilp;
        Param ilp_par;
        ilp_par.setValue("suboptimals", 20);
        ilp.setParameters(ilp_par);


        //for each path the corresponding score
        std::vector<DoubleReal> result_scores;
        //solve the ILP
        clock_t tstart = clock();
        try
        {
          ilp.computeCandidates(spec_graph, result_masses, scores);
          std::cout<<"True time: "<<(clock()-tstart)/(DoubleReal)CLOCKS_PER_SEC<<std::endl;
        }
        catch(IloException &e)
        {
          std::cerr<<e.getMessage()<<std::endl;
          return INTERNAL_ERROR;
        }
        catch(...)
        {
          std::cerr<<"unbekannte exception"<<std::endl;
        }
#endif
        std::cout << "TOTAL TIME SOLVER " << it-spec_map.begin() << ": " << (clock()-time_solver)/CLOCKS_PER_SEC << std::endl;
        std::cout << "TOTAL TIME SOLVER NETTO" << it-spec_map.begin() << ": " << (clock()-time_solver_netto)/CLOCKS_PER_SEC << std::endl;
        time_solver_summed += (clock()-time_solver)/CLOCKS_PER_SEC;
        time_solver_netto_summed += (clock()-time_solver_netto)/CLOCKS_PER_SEC;
        
        std::cout << "TOTAL TIME SOLVER SUMMED" << it-spec_map.begin() << ": " << time_solver_summed << std::endl;
        std::cout << "TOTAL TIME SOLVER NETTO SUMMED" << it-spec_map.begin() << ": " << time_solver_netto_summed << std::endl;
        
        
#ifdef ONLYPATHS
        continue;
#endif

        //POSTPROCESSING OF THE RESULTS
        //the next steps are to resolve the abuigious multi-edges by rescoring the possible AA sequence combinations (all permutations) for each edge.
        //the best results for each multi-edge are then combined to the ultimate candidate. We keep one ultimate candidate for each generated path.
        //Finally we perform a Zhang Similarty scoring of each ultimate candidate to the spectrum and perform the final reranking according to the zhang scores.

        std::multimap<DoubleReal,IdEval::IdAnnot> ranked_candidates, final_cands;
        std::vector<IdEval::IdAnnot> resolved_candidates;

        UInt correct_aa=0;
        DoubleReal correct_percents = 0.0;
        DoubleReal best_correct_percents = 0.0;
        UInt best_correct_residues = 0;
        UInt best_possible_correct_residues = 0;
        DoubleReal best_possible_correct_percents = 0.0;
        UInt best_quick_residues=0;



        DoubleReal time1=clock();
        //IdEval::brute_rescoring(A_map, precision, result_masses, ranked_candidates, *it);
        //IdEval::smart_rescoring(A_map, precision, result_masses, ranked_candidates, *it);
        std::set<std::pair<DoubleReal,IdEval::IdAnnot> >candidate_set;
        PeakSpectrum disc_spec=*it;
        PeakSpectrum norm_spec=*it;
        SqrtMower().filterPeakSpectrum(norm_spec);
        Normalizer().filterPeakSpectrum(norm_spec);
        SqrtMower().filterPeakSpectrum(disc_spec);
        Normalizer().filterPeakSpectrum(disc_spec);
        //SpectrumPreprocessorInspect()(norm_spec);

        Normalizer normalizer;
        Param n_param(normalizer.getParameters());
        n_param.setValue("method", "to_one");
        normalizer.setParameters(n_param);
        normalizer.filterSpectrum(norm_spec);

        //scoring_func.normalizeIntensity(disc_spec);
        //SpectrumPreprocessorInspect()(disc_spec);

        IdEval::smart_rescoring2(A_map, precision, result_masses, candidate_set, disc_spec, scoring_func);

        std::set<AASequence>unique_check;
        for(std::set<std::pair<DoubleReal,IdEval::IdAnnot> >::reverse_iterator setit= candidate_set.rbegin(); setit!=candidate_set.rend(); ++setit)
        {
          std::cout<<setit->second.sequence<<"   ::  "<<setit->first<<std::endl;
          if(unique_check.insert(setit->second.sequence).second)
          {
            resolved_candidates.push_back(setit->second);
          }
        }
        std::cerr<<"rescoring time: "<<(clock()-time1)/CLOCKS_PER_SEC<<std::endl;

        resolved_candidates.resize(std::min(resolved_candidates.size(),Size(400)));


        //--------------------------------------------------------------------------------
        //----------------START: Rescoring Using OpenMS SpectralAlignment-----------------
        //--------------------------------------------------------------------------------
        //with the top 400 now perform spectral alignment for final score:

        TheoreticalSpectrumGenerator t_gen;
        Param params;
        params.setValue("add_isotopes", "true");
        params.setValue("add_losses", "true");
        params.setValue("a_intensity", 0.1);
        params.setValue("x_intensity", 0.2);
        params.setValue("c_intensity", 0.2);
        params.setValue("z_intensity", 0.2);
        params.setValue("relative_loss_intensity", 0.1);


        for(std::vector<IdEval::IdAnnot>::const_iterator res_cand_it=resolved_candidates.begin(); res_cand_it!=resolved_candidates.end(); ++res_cand_it)
        {
          /*
          std::cerr<<res_cand_it->sequence<<std::endl;
          UInt aa_index_mass=(UInt)floor (res_cand_it->sequence[(Size)0].getMonoWeight (Residue::Internal)/precision+0.5);
          if(A_map[aa_index_mass].size()>1)
          {
            params.setValue("add_first_prefix_ion","true");
          }
          else
          {
            params.setValue("add_first_prefix_ion","false");
          }

          t_gen.setParameters(params);
          RichPeakSpectrum tmp;
          t_gen.getSpectrum(tmp, res_cand_it->sequence, 1);
          params.setValue("y_intensity", 0.5);
          params.setValue("b_intensity", 0.5);
          params.setValue("add_losses","false");
          params.setValue("add_isotopes", "true");
          t_gen.setParameters(params);
          t_gen.addPeaks(tmp, res_cand_it->sequence, Residue::AIon,1);
          //t_gen.addPeaks(tmp, res_cand_it->sequence, Residue::CIon,1);
          //t_gen.addPeaks(tmp, res_cand_it->sequence, Residue::ZIon,1);
          //t_gen.addPeaks(tmp, res_cand_it->sequence, Residue::XIon,1);
          t_gen.addPeaks(tmp, res_cand_it->sequence, Residue::YIon, 2);
          //t_gen.addPeaks(tmp, res_cand_it->sequence, Residue::BIon, 2);


          PeakSpectrum tmp_peak_spec;          

          //transform the RichPeakSpectrum into a PeakSpectrum
          for (RichPeakSpectrum::iterator spec_it = tmp.begin(); spec_it != tmp.end(); ++spec_it)
          {
            Peak1D tmp_peak;
            tmp_peak.setIntensity(spec_it->getIntensity());
            tmp_peak.setMZ(spec_it->getMZ());
            tmp_peak_spec.push_back(tmp_peak);
          }

          //IdEval::calibrateSpectrum(norm_spec,res_cand_it->sequence, delta);

          SpectrumAlignmentScore sac;
          Param params_sac;
          params_sac.setValue("use_gaussian_factor","true");
          params_sac.setValue("tolerance",0.3);
          params_sac.setValue("use_linear_factor", "false");
          sac.setParameters(params_sac);
          DoubleReal score = sac(tmp_peak_spec, norm_spec);
          */
          //score=score/res_cand_it->sequence.toUnmodifiedString().length();
          
          DoubleReal score = IdEval::scorePSM(norm_spec, res_cand_it->sequence, scored_types, delta);

          //ZhangSimilarityScore zhang;
          //zhang.setParameters(params_sac);
          //DoubleReal score = zhang(tmp_peak_spec, norm_spec);

          std::cout<<"rescored: "<<res_cand_it->sequence<<"  "<<score<<std::endl;
          final_cands.insert(std::pair<DoubleReal, IdEval::IdAnnot>(score,*res_cand_it));
        }

        //DEBUG stuff
        std::cout<<"fasta output"<<std::endl;
        for(std::vector<IdEval::IdAnnot>::const_iterator res_cand_it=resolved_candidates.begin(); res_cand_it!=resolved_candidates.end(); ++res_cand_it)
        {
          std::cout<<">> "<<res_cand_it-resolved_candidates.begin()<<std::endl;
          std::cout<<res_cand_it->sequence<<std::endl;
        }
        //DEBUG stuff end

        resolved_candidates.clear();
        for(std::multimap<DoubleReal,IdEval::IdAnnot>::reverse_iterator ranked_it=final_cands.rbegin(); ranked_it!=final_cands.rend(); ++ranked_it)
        {
          resolved_candidates.push_back(ranked_it->second);
        }
        //--------------------------------------------------------------------------------
        //----------------END: Rescoring Using OpenMS SpectralAlignment-------------------
        //--------------------------------------------------------------------------------


        //Evaluation:
        //compute the quality of the resolved candidates and compare to the quality of the best possible candidate

        AASequence true_seq =it->getPeptideIdentifications()[0].getHits()[0].getSequence();
        AASequence best_sequence;

        std::cout<<"+++ Results (True Sequence: "<<true_seq<<" )"<<std::endl;
        for(std::vector<IdEval::IdAnnot>::const_iterator res_cand_it=resolved_candidates.begin(); res_cand_it!=resolved_candidates.end(); ++res_cand_it)
        {
          IdEval::predictionQuality(*res_cand_it, true_seq ,correct_aa,correct_percents);
          std::cout<<" ++ "<<res_cand_it->sequence<<"  correct aa's: "<<correct_aa<<"  correct \% : "<<correct_percents<<"     "<<std::endl; //sample_candidates[res_cand_it-resolved_candidates.begin()]<<std::endl;
          if(correct_percents>best_correct_percents)
          {
            best_correct_percents=correct_percents;
            best_correct_residues=correct_aa;
            best_sequence=res_cand_it->sequence;
          }
        }
        total_correct_residues+=best_correct_residues;
        total_residues+=true_seq.size();
        total_residues_precision+=best_sequence.size();
        std::cout<<">> Best Hit Report  "<<counted_spectra<<std::endl;
        //std::cout<<">>> "<<best_correct_percents<<" "<<best_id_rank<<" "<<std::endl;
        std::cout<<">>> "<<best_correct_percents<<std::endl;


        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<">>>>Total performance until here: "<<(DoubleReal)total_correct_residues/total_residues<<std::endl;
        std::cout<<">>>>Total precision until here: "<<(DoubleReal)total_correct_residues/total_residues_precision<<std::endl;
        std::cout<<">++ Average Peptide Length: "<<total_residues/(DoubleReal)(counted_spectra+1)<<std::endl;


        //calculate what would have been possible
        IdEval::IdAnnot tmp_annot(it->getPeptideIdentifications()[0].getHits()[0].getSequence());
        //total_residues+=tmp_annot.sequence.size();

        //UInt quick_computed_corr_res=0;
        time1=clock();
        for (UInt i = 0; i < result_masses.size(); ++i)
        {
          best_quick_residues=std::max(best_quick_residues, IdEval::possiblePredictionQuality(A_map, result_masses[i], precision, true_seq));
          std::cout<<"quickRes: "<<best_quick_residues<<std::endl;

//          std::vector<IdEval::IdAnnot>resolved_annots;
//          IdEval::generateAllAnnotations(A_map, precision, result_masses[i], resolved_annots);
//          for(Size res_ann_id=0; res_ann_id<resolved_annots.size(); ++ res_ann_id)
//          {
//            UInt correct_aa=0;
//            DoubleReal correct_percents=0.0;
//            IdEval::predictionQuality(resolved_annots[res_ann_id], tmp_annot ,correct_aa,correct_percents);
//            if(correct_percents>best_possible_correct_percents)
//            {
//              best_possible_correct_residues = correct_aa;
//              best_possible_correct_percents=correct_percents;
//            }
//          }
        }

        total_quick_residues+=best_quick_residues;
        total_possible_correct_residues+=best_possible_correct_residues;
        std::cerr<<"possible quality computetion: "<<(clock()-time1)/CLOCKS_PER_SEC<<std::endl;

        std::cout<<">> Best possible Hit Report"<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;
        //std::cout<<">>>>Total possible performance until here: "<<(DoubleReal)total_possible_correct_residues/total_residues<<std::endl;
        std::cout<<">>>>Total possible performance until here (QUICK): "<<(DoubleReal)total_quick_residues/total_residues<<std::endl;
        //std::cout<<"-->>comparison local: "<<best_correct_percents<<"  "<<best_possible_correct_percents<<std::endl;
        std::cout<<"-->>comparison local(QUICK): "<<best_correct_percents<<"  "<<best_quick_residues/total_residues<<std::endl;

        ++counted_spectra;

      }
      return EXECUTION_OK;
    }
};


//-----------------------------------------------------------------------------------------------
//------------------------start of main----------------------------------------------------------
//-----------------------------------------------------------------------------------------------
int main(int argc, const char** argv)
{
  DeNovoSequencer tool;
  return tool.main(argc,argv);
}
