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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/TARGETED/PSProteinInference.h>
#include <OpenMS/ANALYSIS/TARGETED/PSLPFormulation.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/SYSTEM/StopWatch.h>

namespace OpenMS
{


  PSLPFormulation::PSLPFormulation():DefaultParamHandler("PSLPFormulation"), solver_(LPWrapper::SOLVER_GLPK)
{
  //model_ = new LPWrapper();
  defaults_.setValue("rt:min_rt",960.,"Minimal rt in seconds.");
	defaults_.setMinFloat("rt:min_rt",0.);
	defaults_.setValue("rt:max_rt",3840.,"Maximal rt in seconds.");
	defaults_.setMinFloat("rt:max_rt",0.);
	defaults_.setValue("rt:rt_step_size",30.,"rt step size in seconds.");
	defaults_.setMinFloat("rt:rt_step_size",1.);
  
  defaults_.setValue("thresholds:min_protein_probability",0.2,"Minimal protein probability for a protein to be considered in the ILP");
  defaults_.setValue("thresholds:min_protein_id_probability",0.95,"Minimal protein probability for a protein to be considered identified.");
	defaults_.setValue("thresholds:min_pt_weight",0.5,"Minimal pt weight of a precursor");
  defaults_.setValue("thresholds:min_pred_pep_prob",0.5,"Minimal predicted peptide probability of a precursor");
	defaults_.setMinFloat("thresholds:min_pred_pep_prob",0.);
	defaults_.setMaxFloat("thresholds:min_pred_pep_prob",1.);

	defaults_.setMinFloat("thresholds:min_pt_weight",0.);
	defaults_.setMaxFloat("thresholds:min_pt_weight",1.);
	defaults_.setValue("thresholds:min_rt_weight",0.5,"Minimal rt weight of a precursor");
	defaults_.setMinFloat("thresholds:min_rt_weight",0.);
	defaults_.setMaxFloat("thresholds:min_rt_weight",1.);
	defaults_.setValue("mz_tolerance",25.,"Allowed precursor mass error tolerance in ppm.");
	defaults_.setValue("rt_tolerance",300,"Allowed precursor rt deviation in seconds.");

  defaults_.setValue("combined_ilp:k1",0.2,"combined ilp: weight for z_i");
  defaults_.setMinFloat("combined_ilp:k1",0.);
  //	defaults_.setMaxFloat("combined_ilp:k1",1.);
  defaults_.setValue("combined_ilp:k2",0.2,"combined ilp: weight for x_j,s*int_j,s");
  defaults_.setMinFloat("combined_ilp:k2",0.);
  //	defaults_.setMaxFloat("combined_ilp:k1",1.);
  defaults_.setValue("combined_ilp:k3",0.4,"combined ilp: weight for -x_j,s*w_j,s");
  defaults_.setMinFloat("combined_ilp:k3",0.);
  //	defaults_.setMaxFloat("combined_ilp:k1",1.);
	
	defaultsToParam_();
  
}

PSLPFormulation::~PSLPFormulation()
{
  //delete model_;
}

void PSLPFormulation::createAndSolveCombinedLPFeatureBased_(const FeatureMap<>& features,std::vector<std::vector<DoubleReal> >& intensity_weights,
                                                    std::set<Int>& charges_set,std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
                                                    std::vector<IndexTriple>& variable_indices,std::vector<Int>& solution_indices,
                                                    UInt ms2_spectra_per_rt_bin,Size number_of_scans, Size step_size)
{
  DoubleReal k2 = param_.getValue("combined_ilp:k2");
  std::cout << "k2: "<<k2 <<std::endl;

	//#define DEBUG_OPS
  model_ = new LPWrapper();
  model_->setSolver(solver_);
	Int counter = 0;

	std::cout << "Feature Based: Build model: first objective"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	// add objective function
	///////////////////////////////////////////////////////////////////////
	model_->setObjectiveSense(LPWrapper::MAX); // maximize
	// max \sum_j x_jk * signal_jk
	//                    column_index, feature_index,scan
	
	//

  // first get maximimal intensity, as we use it for normalization
  DoubleReal max_int =0.;
	for(Size i = 0; i < features.size(); ++i)
		{
			// first check if charge state is allowed
			// charge not in "ChargeFilter" list
#ifdef DEBUG_OPS
			std::cout << "feat: "<<i <<" charge "<<features[i].getCharge() << std::endl;
#endif
			if (charges_set.count(features[i].getCharge())<1) continue;
      if((DoubleReal)features[i].getMetaValue("msms_score") > max_int) max_int = features[i].getMetaValue("msms_score");
    }
	std::cout << features.size() << " features.\n";
	for(Size i = 0; i < features.size(); ++i)
		{
			// first check if charge state is allowed
			// charge not in "ChargeFilter" list
#ifdef DEBUG_OPS
			std::cout << "feat: "<<i <<" charge "<<features[i].getCharge() << std::endl;
#endif
			if (charges_set.count(features[i].getCharge())<1) continue;
			if(mass_ranges[i].size()==0)
				{
					std::cout  << "No mass ranges for "<<features[i].getRT() << " "<<features[i].getMZ()<<std::endl;
				}
#ifdef DEBUG_OPS
			if(mass_ranges[i].size() > 0)
				{
					std::cout << "start_scan "<< mass_ranges[i][0].first << " ?= "<<features[i].getQuality(0)
										<< " stop scan "<< (mass_ranges[i].end()-1)->first<< " ?= "	<< features[i].getQuality(1)-1<<std::endl;
				}
#endif

			DoubleReal msms_score = features[i].getMetaValue("msms_score");

			Size c = 0;
			// go through all rts of the current feature
			for(Size s_idx = 0; s_idx < mass_ranges[i].size();s_idx+=2) // walk in size two steps
 			{
					Size s = mass_ranges[i][s_idx].first;


					////////////////////////////////////////////////////////////////////////

						

#ifdef DEBUG_OPS
					std::cout << "add column "<<counter << std::endl;
#endif
					IndexTriple triple;
					triple.feature = i;
					triple.scan = s;
          Size index = model_->addColumn();
					triple.variable = index;
					variable_indices.push_back(triple);
          model_->setColumnBounds(index,0,1,LPWrapper::DOUBLE_BOUNDED);
					model_->setColumnType(index,LPWrapper::BINARY); // binary variable
					model_->setColumnName(index,(String("x_")+i+","+s).c_str());					
          //#ifdef DEBUG_OPS	
          std::cout << "feat "<<i << " scan "<< s<< " "
                    << intensity_weights[i][c] <<" msms_score "
                    << msms_score <<" max_int "<<max_int
                    << " obj: "<<intensity_weights[i][c]*msms_score
                    << " anderes obj.: "<< intensity_weights[i][c]*msms_score/max_int
                    <<std::endl;
          //           intensity_weights[i][c]*msms_score/max_int * (*x_j_s_)[counter];
          //#endif
          model_->setObjective(index,intensity_weights[i][c]* (DoubleReal)features[i].getMetaValue("msms_score"));
					//cmodel_->setObjective(counter,intensity_weights[i][c]);
					++counter;
					//					if(intensity_weights[i][c]*msms_score > max_int) max_int = intensity_weights[i][c]*msms_score;
					if(msms_score > max_int) max_int = msms_score;
					//if(intensity_weights[i][c] > max_int) max_int = intensity_weights[i][c];
					++c;
					//std::cout <<"max_int "<<max_int<<std::endl;
				}
			
		}
	std::cout << "now normalize the objective values by "<<max_int<<"\n";
	//	normalize and invert objectives
	for(Int i = 0; i < counter; ++i)
		{
			//			std::cout <<i <<" "<< cmodel_->objective(i)<<"\t";
			model_->setObjective(i,k2*model_->getObjective(i)/max_int);
			//			std::cout <<cmodel_->objective(i)<<"\n";
		}
	
	
	///////////////////////////////////////////////////////////////////////
	// add constraints
	///////////////////////////////////////////////////////////////////////
#ifdef DEBUG_OPS	
	std::cout << "and now the constraints:"<<std::endl;
#endif
	///////////////////////////////////////////////////////////////////////
	// 1: ensure that each precursor is acquired maximally once
	///////////////////////////////////////////////////////////////////////
	addPrecursorAcquisitionNumberConstraint_(variable_indices,features.size(),1);
	///////////////////////////////////////////////////////////////////////
	// 2: do not exceed rt bin capacity
	///////////////////////////////////////////////////////////////////////
	addRTBinCapacityConstraint_(variable_indices,number_of_scans,ms2_spectra_per_rt_bin);
	///////////////////////////////////////////////////////////////////////
	// 3: add step size constraint
	///////////////////////////////////////////////////////////////////////
  if(step_size > 0) addStepSizeConstraint_(variable_indices,step_size);
	solveILP(solution_indices);  
}

  
void PSLPFormulation::createAndSolveILPFeatureBased_(const FeatureMap<>& features,std::vector<std::vector<DoubleReal> >& intensity_weights,
                                                     std::set<Int>& charges_set,std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
                                                     std::vector<IndexTriple>& variable_indices,std::vector<int>& solution_indices,
                                                     UInt ms2_spectra_per_rt_bin,
                                                     Size number_of_scans)
{
	Int counter = 0;
  model_ = new LPWrapper();
  model_->setSolver(solver_);
  //#define DEBUG_OPS	
	std::cout << "Feature Based: Build model: first objective"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	// add objective function
	///////////////////////////////////////////////////////////////////////
	model_->setObjectiveSense(LPWrapper::MAX); // maximize
	// max \sum_j x_jk * signal_jk
	//                    column_index, feature_index,scan
	
	//

	for(Size i = 0; i < features.size(); ++i)
		{
			// first check if charge state is allowed
			// charge not in "ChargeFilter" list
#ifdef DEBUG_OPS
			std::cout << "feat: "<<i <<" charge "<<features[i].getCharge() << std::endl;
#endif
			if (charges_set.count(features[i].getCharge())<1) continue;
			if(mass_ranges[i].size()==0) continue;
#ifdef DEBUG_OPS
			if(mass_ranges[i].size() > 0)
				{
					std::cout << "start_scan "<< mass_ranges[i][0].first << " ?= "<<features[i].getQuality(0)
										<< " stop scan "<< (mass_ranges[i].end()-1)->first<< " ?= "	<< features[i].getQuality(1)-1<<std::endl;
				}
#endif

		
			Size c = 0;
			// go through all rts of the current feature
			for(Size s_idx = 0; s_idx < mass_ranges[i].size();s_idx+=2) // walk in size two steps
				{
					Size s = mass_ranges[i][s_idx].first;


					////////////////////////////////////////////////////////////////////////

						
					
#ifdef DEBUG_OPS
					std::cout << "add column "<<counter << std::endl;
#endif
					IndexTriple triple;
					triple.feature = i;
					triple.scan = s;
          Int index = model_->addColumn();
          triple.variable = index;
					variable_indices.push_back(triple);
          
          std::cout << index << " variable index"<<std::endl;
          model_->setColumnBounds(index,0,1,LPWrapper::DOUBLE_BOUNDED);
          model_->setColumnType(index,LPWrapper::BINARY); // binary variable
					model_->setColumnName(index,(String("x_")+i+","+s));
          //#ifdef DEBUG_OPS	
					std::cout << "feat "<<i << " scan "<< s << " intensity_weight "
										<< intensity_weights[i][c] <<std::endl;
          //#endif
					model_->setObjective(index,intensity_weights[i][c]);
					++counter;
					++c;
				}
		}
	
	///////////////////////////////////////////////////////////////////////
	// add constraints
	///////////////////////////////////////////////////////////////////////
#ifdef DEBUG_OPS	
	std::cout << "and now the constraints:"<<std::endl;
#endif
	///////////////////////////////////////////////////////////////////////
	// 1: ensure that each precursor is acquired maximally once
	///////////////////////////////////////////////////////////////////////
#ifdef DEBUG_OPS	
	std::cout << "first the number of times a precursors is acquired"<<std::endl;
#endif
  addPrecursorAcquisitionNumberConstraint_(variable_indices,features.size(),1);

	///////////////////////////////////////////////////////////////////////
	// 2: do not exceed rt bin capacity
	///////////////////////////////////////////////////////////////////////
#ifdef DEBUG_OPS	
	std::cout << "and now the rt bin capacity"<<std::endl;
	std::cout << ms2_spectra_per_rt_bin << " rt bin capacity"<<std::endl;
#endif
  addRTBinCapacityConstraint_(variable_indices,number_of_scans,ms2_spectra_per_rt_bin);

#ifdef DEBUG_OPS	
	model_->writeProblem("/home/zerck/data/tmp/test_pis_problem.mps","MPS");
#endif
	
	solveILP(solution_indices);
}

void PSLPFormulation::addPrecursorAcquisitionNumberConstraint_(std::vector<IndexTriple>& variable_indices,
                                                               Size number_of_features,UInt number_of_msms_per_precursor)
{
	///////////////////////////////////////////////////////////////////////
	// ensure that each precursor is acquired maximally once
	///////////////////////////////////////////////////////////////////////
  //	std::cout << "now the number of times a precursors is acquired"<<std::endl;
	Size j = 0;
	for(Size i = 0; i < number_of_features;++i)
		{
			Size start = j;
			while(j < variable_indices.size() && variable_indices[j].feature == i)
				{
#ifdef DEBUG_OPS
					std::cout << j << " "<<variable_indices[j].variable << " "
										<< variable_indices[j].feature << " "
										<< variable_indices[j].scan<<std::endl;
#endif
					++j;
				}

			Size stop = j;
      std::vector<DoubleReal> entries(stop-start);
			std::vector<Int> indices(stop-start);
#ifdef DEBUG_OPS
			std::cout << "feature "<<i <<" ";//<<features[i].getMZ() <<" "<<features[i].getRT()<<" ";
			std::cout << stop-start<<"variables in equation\n";
#endif
			Size c = 0;
			for(Size k = start; k < stop; ++k)
				{
					entries[c] = 1.;
					indices[c] = variable_indices[k].variable;
					//					std::cout << j<<" "<<indices[j]<<std::endl;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			String name = "PREC_ACQU_LIMIT_" + String(i);
			if(stop-start > 0)
        {
          model_->addRow(indices,entries,name,0,(Int)number_of_msms_per_precursor,LPWrapper::UPPER_BOUND_ONLY);
        }
#ifdef DEBUG_OPS
			std::cout << stop-start << " "<<name<<std::endl;
			std::cout << "added row"<<std::endl;
#endif
			
		}
}

void PSLPFormulation::addRTBinCapacityConstraint_(std::vector<IndexTriple>& variable_indices,
                                                  Size max_rt_index,UInt ms2_spectra_per_rt_bin)
{
	//////////////////////////////////////////////////////////////////////////
	// constraint : rt_bin_capacity
	//////////////////////////////////////////////////////////////////////////
	// sort variable_indices according to their scan number
	sort(variable_indices.begin(),variable_indices.end(),ScanLess());
	Size j = 0;
	for(Size i = 0; i < max_rt_index;++i)
		{
			// first determine number of indices:
			Size start = j;
			while(j < variable_indices.size() && variable_indices[j].scan == i)
				{
					++j;
				}
			// no feature occuring in this scan
			if(start == j) continue;

			Size stop = j;
			Size c = 0;			
      std::vector<double> entries(stop-start);
      std::vector<int> indices(stop-start);
			for(Size s = start; s < stop; ++s)
				{
					entries[c] = 1.;
					indices[c] = variable_indices[s].variable;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
      model_->addRow(indices,entries,(String("RT_CAP")+i),0,ms2_spectra_per_rt_bin,LPWrapper::UPPER_BOUND_ONLY);// only upper bounded problem -> lower bound is ignored
#ifdef DEBUG_OPS
			std::cout << "added row"<<std::endl;
#endif
			
		}

}

void PSLPFormulation::addStepSizeConstraint_(std::vector<IndexTriple>& variable_indices,UInt step_size)
{
 	//////////////////////////////////////////////////////////////////////////
	// constraint : step size
	//////////////////////////////////////////////////////////////////////////
	std::cout << "and now the step size"<<std::endl;
	std::cout << step_size << " step size"<<std::endl;

  std::vector<DoubleReal> entries(variable_indices.size(),1.);
  std::vector<Int> indices(variable_indices.size());
	for(Size i = 0; i < variable_indices.size();++i)
		{
			indices[i] = i;
    }      
#ifdef DEBUG_OPS
  std::cout << "\nadd row "<<std::endl;
#endif
  model_->addRow(indices,entries,"step_size",0,step_size,LPWrapper::UPPER_BOUND_ONLY);// only upper bounded problem -> lower bound is ignored
#ifdef DEBUG_OPS
  std::cout << "added row"<<std::endl;
#endif
		
}

void PSLPFormulation::updateStepSizeConstraint(Size num_realized_precursors,UInt step_size)
{
  Int row_index =  model_->getRowIndex("step_size");
  std::cout << "index "<<row_index<<std::endl;
  std::cout << "step_size:row_upper before "<<model_->getRowUpperBound(row_index);
  model_->setRowBounds(row_index,0,(DoubleReal)(step_size+num_realized_precursors),LPWrapper::UPPER_BOUND_ONLY);
  std::cout << " and after "<<model_->getRowUpperBound(row_index)<<std::endl;
}


void PSLPFormulation::solveILP(std::vector<int>& solution_indices,Int /*iteration*/,bool /*targeted*/,Size /*rt_bin_capacity*/)
{
  solution_indices.clear();
  LPWrapper::SolverParam param;
//   model_->writeProblem(String("test_prob_")+String(solver_)+String(".mps"),"MPS");
  model_->solve(param);

  for (Int column = 0; column < model_->getNumberOfColumns(); ++column)
  {
    double value = model_->getColumnValue(column);
    std::cout << value << " "<< model_->getColumnType(column) << std::endl;
    if ((fabs(value) > 0.5 && model_->getColumnType(column) == LPWrapper::BINARY) ||
      (fabs(value) > 0.5 && model_->getColumnType(column) == LPWrapper::INTEGER))
    {
          //          std::cout << value << " "<< model_->getColumnType(column) << std::endl;      
#ifdef DEBUG_OPS	
      std::cout << model_->getColumnName(column) << " is in optimal solution" << std::endl;
#endif
      solution_indices.push_back((int) column);
    }
  }
}
  

DoubleReal PSLPFormulation::getObjectiveValue(Int index)
{
  return model_->getObjective(index);
}


void PSLPFormulation::updateFeatureILPVariables(FeatureMap<>& new_features,
                                                std::vector<IndexTriple>& variable_indices,
                                                std::map<Size,std::vector<String> >& feature_constraints_map)
{
	std::cout << "update feature ilp variables"<<std::endl;
	DoubleReal min_rt = param_.getValue("rt:min_rt");
	DoubleReal max_rt = param_.getValue("rt:max_rt");
	DoubleReal rt_step_size = param_.getValue("rt:rt_step_size");

	Size max_index = (Size)ceil((max_rt-min_rt) / rt_step_size);
	for(Size f = 0; f < new_features.size();++f)
		{
			Size f_index = new_features[f].getMetaValue("feature_index");
			// find corresponding variables
			Size f_v_idx = 0;
			while(f_v_idx < variable_indices.size() && f_index != variable_indices[f_v_idx].feature)
				{
					++f_v_idx;
				}
			if(f_v_idx == variable_indices.size()) std::cout << "This should not happen!"<<std::endl;
			else
				{
					// now set the corresponding variables in the ILP to 1 as this peptide is already acquired
					// find the right spectrum
					Size rt_index = std::min((Size)std::max(0.,ceil((new_features[f].getRT()-min_rt)/rt_step_size)) ,max_index);
					bool existing = false;
					while(f_v_idx < variable_indices.size() && variable_indices[f_v_idx].feature == f_index)
						{
							if(variable_indices[f_v_idx].scan == rt_index)
								{
									existing = true;
                  model_->setColumnBounds(variable_indices[f_v_idx].variable,1.,
                                          model_->getColumnUpperBound(variable_indices[f_v_idx].variable),
                                          LPWrapper::FIXED);
									std::cout << "set column lower von "<<model_->getColumnName(variable_indices[f_v_idx].variable)
														<< " auf "<<model_->getColumnLowerBound(variable_indices[f_v_idx].variable) <<std::endl;
									//									cmodel_->setObjective(variable_indices[f].variable,score);
									break;
 								}
							++f_v_idx;
						}
				}
			std::map<Size,std::vector<String> >::iterator c_iter = feature_constraints_map.find(f);
			if(c_iter!=feature_constraints_map.end())
				{
					for(Size c = 0; c < c_iter->second.size();++c)
						{
							//std::cout <<"constraint name: "<< c_iter->second[c]<<std::endl;
							Int row = model_->getRowIndex(c_iter->second[c]);
							if(row != -1) // if constraint is still existing
								{
									//									std::cout << "delete row "<<row<<std::endl;
									model_->deleteRow(row);
									//std::cout << "deleted row.\n";
								}
						}
				}
			
		}// 	for(Size f = 0; f < new_features.size();++f)
  std::cout << "update feature ilp variables --->done"<<std::endl;
}


void PSLPFormulation::updateCombinedILP(FeatureMap<>& features,
                                        PrecursorIonSelectionPreprocessing& preprocessed_db,
                                        std::vector<IndexTriple>& variable_indices,
                                        std::vector<String>& new_protein_accs,
                                        std::vector<String>& protein_accs,
                                        PSProteinInference& prot_inference,
                                        Size& variable_counter,
                                        std::map<String,std::vector<Size> >& protein_feature_map,
                                        Feature& new_feature, std::map<String,Size>& protein_variable_index_map)
{
  std::cout << "Update Combined ILP."<<std::endl;
	DoubleReal min_prot_coverage = param_.getValue("thresholds:min_protein_id_probability");
	DoubleReal min_pt = param_.getValue("thresholds:min_pt_weight");
	DoubleReal min_rt_weight = param_.getValue("thresholds:min_rt_weight");
  DoubleReal min_pred_pep_weight = param_.getValue("thresholds:min_pred_pep_prob");
  DoubleReal mz_tolerance = param_.getValue("mz_tolerance");
  bool use_detectability = true;//param_.getValue("use_detectability") == "true" ? true : false;
	DoubleReal min_protein_probability = param_.getValue("thresholds:min_protein_probability");
  DoubleReal k1 =  param_.getValue("combined_ilp:k1");
  std::cout << "k1 "<<k1<<std::endl;
	std::cout << "min_protein_probability = "<<min_protein_probability<< " min_prot_coverage "
            << min_prot_coverage <<std::endl;
	std::cout << "parsed all parameters"<<std::endl;
	std::cout << "log(1.-min_prot_coverage) = "<<log(1.-min_prot_coverage)<<std::endl;

	sort(variable_indices.begin(),variable_indices.end(),FeatureIndexLess());
  StopWatch timer;
  timer.start();
	if(new_feature.getPeptideIdentifications().size() > 0
		 && new_feature.getPeptideIdentifications()[0].getHits().size()>0)
		{
			// if a selected feature yielded a peptide id, the peptide probability needs to be considered in the protein constraint
			DoubleReal pep_score = new_feature.getPeptideIdentifications()[0].getHits()[0].getScore();
			Size index = new_feature.getMetaValue("variable_index");
			const std::vector<String>& accs = new_feature.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions();
			// check all proteins that were already detected (only there we need to update a constraint)
			for(Size pa = 0; pa < protein_accs.size();++pa)
				{
					if(find(accs.begin(),accs.end(),protein_accs[pa]) == accs.end())  continue;
					DoubleReal weight= 1.;
					Int row= model_->getRowIndex((String("PROT_COV_")+protein_accs[pa]).c_str());
          std::cout << protein_accs[pa] << " index "<<row << " "<<model_->getElement(row,index)<<std::endl;
					if(model_->getElement(row,index) != 0.)
						{
// 							std::cout << "getElement("<<protein_accs[pa]<<","<<index<<")="
// 												<< cmodel_->getElement(row,index) << "\t";
							if(fabs(pep_score*weight - 1.)< 0.000001)
								{
									model_->setElement(row,index,-log( 0.000001)/log(1.-min_prot_coverage)); // pseudocount to prevent entering inf
									std::cout << protein_accs[pa] <<" setElement("<<row<<","<<model_->getColumnName(index)<<"=  "
														<< model_->getElement(row,index)<<std::endl;
								}
							else
								{
									model_->setElement(row,index,-log(1.-pep_score*weight)/log(1.-min_prot_coverage));
									std::cout << protein_accs[pa]<<" setElement("<<row<<","<<model_->getColumnName(index)<<"=  "
														<< model_->getElement(row,index)<<std::endl;
								}
// 							std::cout << "getElement("<<protein_accs[pa]<<","<<index<<")="
// 												<< cmodel_->getElement(row,index) << "\n";
						}
// 					else std::cout << "getElement("<<protein_accs[pa]<<","<<index<<") ="
// 												 << cmodel_->getElement(row,index) << " ==? 0.\n";
				}//for(Size pa = 0; pa < protein_accs.size();++pa)
			
			std::cout << "Now search for matching features for "<<accs.size() <<" proteins."<<std::endl;
			for(Size prot = 0; prot < accs.size();++prot)
				{
					// first enter the new feature to the corresponding proteins
					std::map<String,std::vector<Size> >::iterator prot_var_iter = protein_feature_map.find(accs[prot]);
					std::cout << "now enter "<<new_feature.getRT() << "  "<<new_feature.getMZ()
										<< " with variable index "
										<< index
										<< " to prot_variable_map "<<accs[prot] << std::endl;
					if(prot_var_iter == protein_feature_map.end())
						{
							std::vector<Size> vec;
							vec.push_back((Size)new_feature.getMetaValue("feature_index"));
							protein_feature_map.insert(std::make_pair(accs[prot],vec));
						}
					else prot_var_iter->second.push_back((Size)new_feature.getMetaValue("feature_index"));
					// if protein prob exceeds min threshold
					if((prot_inference.getProteinProbability(accs[prot]) > min_protein_probability)
             && (prot_inference.isProteinInMinimalList(accs[prot])))
						{
							std::map<String,std::vector<DoubleReal> >::const_iterator map_iter = preprocessed_db.getProteinPTMap().find(accs[prot]);
							if(map_iter !=  preprocessed_db.getProteinPTMap().end())
								{
									std::vector<String>::iterator prot_acc_iter = find(protein_accs.begin(),protein_accs.end(),accs[prot]);
									if(prot_acc_iter == protein_accs.end())
										{
                      // enter new variable y_i
                      if(prot_inference.getProteinProbability(accs[prot]) >= min_prot_coverage )
                        {
                          new_protein_accs.push_back(accs[prot]);
                          updateObjFunction_(accs[prot],features,preprocessed_db,variable_indices);
                        }
											protein_accs.push_back(accs[prot]);
											// insert protein to ILP
											// first add penalty
											// insert penalty variable
                      Int index = model_->addColumn();
                      model_->setColumnName(index,(String("y_")+map_iter->first).c_str());
                      model_->setColumnBounds(index,0.,1.,LPWrapper::DOUBLE_BOUNDED);
                      //cmodel_->setColumnIsInteger(variable_counter,true); //try without integer constraint
											model_->setObjective(index,k1);
											protein_variable_index_map.insert(make_pair(map_iter->first,index));
                      std::cout << "entered protein variable "<< (String("y_")+map_iter->first)
                                << "\tit consists of:"
                                << std::endl;
											//	std::cout << variable_counter <<" penalty for "<<map_iter->first << std::endl;
											++variable_counter;

											std::vector<Int> indices;
											std::vector<DoubleReal> entries;


											// now go through protein_feature_map and enter them to
											std::map<String,std::vector<Size> >::iterator prot_var_iter = protein_feature_map.find(accs[prot]);
											if(prot_var_iter != protein_feature_map.end())
												{
													// TODO: problem: they need not all correspond to the current feature
													for(Size var = 0; var < prot_var_iter->second.size();++var)
														{
															Int variable = features[prot_var_iter->second[var]].getMetaValue("variable_index");
															DoubleReal score = features[prot_var_iter->second[var]].getPeptideIdentifications()[0].
																getHits()[0].getScore();
															bool found_index = false;
															for(Size idx = 0; idx < indices.size();++idx)
																{
																	// TODO: really compare here with the variable index?
																	// isn't the feature index more interesting here? think about it!
																	if(indices[idx] == variable)
																		{
																			DoubleReal weight = 1.;
																			if(fabs(score*weight- 1.) < 0.000001) // pseudocount to prevent entering -inf
																				{
																					entries[idx] = -log( 0.000001) / log(1.-min_prot_coverage);
																				}
																			// if the current entry has a better peptide id, take it instead of the old one 
																			else if(entries[idx] > -log(1.-score*weight))
																				{
																					entries[idx] = -log(1.-score*weight) / log(1.-min_prot_coverage);
																				}
                                      std::cout << entries[idx]<<" "<<score*weight<<"\n";
																			found_index = true;
																			break;
																		}
																}
															if(!found_index)
																{
																	std::cout << "feature that was already measured is entered into prot_cov "
																						<< prot_var_iter->first << " matches variable no. "
																						<< variable << std::endl;
																	indices.push_back(variable);
																	DoubleReal score = features[prot_var_iter->second[var]].getPeptideIdentifications()[0].
                                    getHits()[0].getScore();
																	DoubleReal weight = 1.;
                                  std::cout <<features[prot_var_iter->second[var]].getMZ()<<" "<<features[prot_var_iter->second[var]].getRT()<<std::endl;
																	if(fabs(weight * score- 1.) < 0.000001)
																		{
																			entries.push_back(-log( 0.000001) / log(1.-min_prot_coverage));
																			std::cout <<"pseudocount needed: "<< score
																								<<" * "<<weight << " = "<<weight*score<<" -> "
																								<<-log( 0.000001) << std::endl;
																		}
																	else
																		{
																			entries.push_back(-log(1.-weight*pep_score) / log(1.-min_prot_coverage));
																			std::cout << score<<" * "<<weight << " = "<<weight*score<<" -> "
																								<<-log(1.-weight*score) << std::endl;
																		}
                                  std::cout << entries.back()<<" "<<score*weight<<"\n";
																}
														}
												}//	if(prot_var_iter != protein_feature_map.end())

											// then go through all (tryptic) peptides of this protein
											std::cout <<"protein "<<map_iter->first << " with "<< map_iter->second.size()
																<< " entries in pt_map"<<std::endl;
											const std::vector<DoubleReal>& masses = preprocessed_db.getMasses(map_iter->first);
											for(Size p = 0;p < map_iter->second.size();++p)
												{
// 													std::cout <<preprocessed_db.getRT(map_iter->first,p)
// 																		<<" "<< masses[p]<< " "<<map_iter->second[p]
// 																		<<std::endl;
//													std::cout << map_iter->second[p] << " > "<<min_pt<<std::endl;
													if(map_iter->second[p] > min_pt)
														{
															std::vector<Int> matching_features;
															// go through all features
															for(Size f = 0; f < features.size();++f)
																{
																	if(features[f].getMetaValue("fragmented")=="true") continue;
																	// 			std::cout << masses[p] << " - "<< features[f].getMZ() << " -> "
																	// 																		<< fabs(masses[p] - features[f].getMZ())/masses[p] * 1e06<<" < "
																	// 																		<< mz_tolerance<<std::endl;
																	// and check if they match the current peptide
																	if(fabs(masses[p] - features[f].getMZ())/masses[p] * 1e06 < mz_tolerance)
																		{
																			// we need to compare the rt of each scan separately, not of the whole feature
 																			DoubleReal rt_weight = preprocessed_db.getRTProbability(map_iter->first,p,features[f]);
                                      //                                    std::cout << "rt_weight : "<<rt_weight <<" >? "<< min_rt_weight
                                      //        << std::endl;
                                      
																			if(rt_weight > min_rt_weight)
																				{
																					
 																					std::cout << features[f].getMZ() << " "<<features[f].getRT()<<" "
 																										<< rt_weight <<" is matching peptide "<<p <<"\n";
                                          
																			
																					// store all matching features to this peptide
																					matching_features.push_back((Int)f);
																					// if yes: add as variable to the protein coverage constraint
																					// find corresponding variables
																					Size f_v_idx = 0;
																					while(f_v_idx < variable_indices.size() && f != variable_indices[f_v_idx].feature)
																						{
																							++f_v_idx;
																						}
																					if(f_v_idx == variable_indices.size())
																						{
																							std::cout << features[f].getMZ() << " "<< features[f].getRT()<<" "
																												<< /*rt_weight <<*/" is matching peptide "<< p
																												<< ", but not existing in variable indices???"
																												<< "--->This should not happen!"<<std::endl;
																						}
																					else
																						{
																							// TODO: enter here really all scans for this feature?
																							// or only the ones with rt-weight > min_rt_weight
																							while(f_v_idx < variable_indices.size() && f == variable_indices[f_v_idx].feature)
																								{
																									bool found_index = false;
																									for(Size idx = 0; idx < indices.size();++idx)
																										{
																											if(indices[idx] == (Int)f_v_idx)
																												{
																													found_index = true;
																													break;
																												}
																										}
																									if(!found_index)
																										{
                                                      // weight is detectability * rt_weight
                                                      DoubleReal dt;
                                                      if(use_detectability)
                                                        {
                                                          dt = map_iter->second[p];
                                                        }
                                                      else
                                                        {
                                                          dt = 1.;
                                                        }
                                                      DoubleReal weight = dt * rt_weight;//preprocessed_db.getRTDTProbability(dt*curr_rt_weight);
                                                          
                                                      std::cout << dt << " * " << rt_weight
                                                                << " = "<< dt*rt_weight
                                                                << " -> "
                                                                << weight
                                                                << " -> "<< log(1.-weight) / log(1.-min_prot_coverage)
                                                                << " obj: "<<model_->getObjective(f_v_idx)
                                                                << " columnname: "<<model_->getColumnName(f_v_idx)
                                                                << std::endl;
                                                      
                                                      std::cout << weight << " <? "<< min_pred_pep_weight<<"? "
                                                                << (weight < min_pred_pep_weight) <<"\n";
                                                      if(weight < min_pred_pep_weight)
                                                        {
                                                          ++f_v_idx;
                                                          continue;
                                                        }
                                                      indices.push_back((Int)f_v_idx);
                                                      if(fabs(1.-weight) <  0.000001)
                                                        {
                                                          entries.push_back(-log(0.000001)/ log(1.-min_prot_coverage));
                                                        }
                                                      else entries.push_back(-log(1.-weight)/ log(1.-min_prot_coverage));
                                                      if(features[f].metaValueExists("pred_pep_prob"))
                                                        {
                                                          std::cout << "old pred_pep_prob "<<features[f].getMetaValue("pred_pep_prob")
                                                                    << " vs. new "<<weight<<std::endl;
                                                        }
                                                      features[f].setMetaValue("pred_pep_prob",weight);
                                                      // enter this constraint for the this feature into the feature_constraint_map
                                                      //																													std::map<Size,std::vector<String> >::iterator f_c_iter = feature_constraints_map.find(f_v_idx);
                                                    
																										}
																									++f_v_idx;
																								}// while(f_v_idx < variable_indices.size() && f == variable_indices[f_v_idx].feature)
																						}//else
                                        }//if(rt_weight > min_rt_weight)
																		}//if(fabs(masses[p] - features[f].getMZ())/masses[p] * 1e06 < mz_tolerance)
																}//for(Size f = 0; f < features.size();++f)
														}//if(map_iter->second[p] > min_pt)
												}//for(Size p = 0;p < map_iter->second.size();++p)

											// now enter protein variable y_i
											indices.push_back(protein_variable_index_map[map_iter->first]);
											entries.push_back(1.);
											//Int i = distance(preprocessed_db.getProteinPTMap().begin(),map_iter);			
#ifdef DEBUG_OPS
											std::cout << "\nadd row ";
											std::cout << (String("PROT_COV_")+map_iter->first) <<"\t"<<(String("PROT_COV_")+i).c_str();
											std::cout << " "<<indices.size()<< " entries "<<(indices[0])<< " "<<(entries[0])<<std::endl;
#endif
											
											// at the moment we want a certain protein coverage
											model_->addRow(indices,entries,String("PROT_COV_")+map_iter->first,0,0.,LPWrapper::UPPER_BOUND_ONLY);
#ifdef DEBUG_ILP
                      // print constraint:
                      std::cout <<" the constraint: \n";
                      for(Size idx = 0; idx < indices.size();++idx)
                        {
                          std::cout <<indices[idx]<<"\t";
                        }
                      std::cout << "\n";
                      for(Size idx = 0; idx < entries.size();++idx)
                        {
                          std::cout <<entries[idx]<<"\t";
                        }
                      std::cout << "\n";
#endif
											//											std::cout << "added row"<<std::endl;
										}//if(prot_acc_iter == protein_accs.end())
								}
							else
								{
									std::cout << "Protein not present in preprocessed db: "<<accs[prot]<<std::endl;
								}
						}
				}
		}
// 	else if(!not_solve)
  else 
    {
			// if a selected feature didn't yield any peptide id, all concerned
			// protein constraints need to be updated (feature has to be deleted from this constraint)
			Size index = new_feature.getMetaValue("variable_index");
			Size f_index = new_feature.getMetaValue("feature_index");
			
			std::cout << "new feature had no peptide id"<<std::endl;
			Size f_v_idx = 0;
			while(f_v_idx < variable_indices.size() && variable_indices[f_v_idx].feature != f_index)
				{
					++f_v_idx;
				}

			// TODO: maybe we need to update old protein constraints, as the selected feature didn't yield any
			// id
      std::set<String> updated_constraints;
			for(Size pa = 0; pa < protein_accs.size();++pa)
				{
					Int row= model_->getRowIndex(String("PROT_COV_")+protein_accs[pa]);
										
					Size f_v_idx2 = f_v_idx;
					while(f_v_idx2 < variable_indices.size() && f_index == variable_indices[f_v_idx2].feature)
						{
							index = variable_indices[f_v_idx2].variable;
							if(model_->getElement(row,index) != 0.)
								{
                  updated_constraints.insert(protein_accs[pa]);
									std::cout << "getElement("<<protein_accs[pa]<<","<<index<<")="
														<< model_->getElement(row,index) << "\t";
									model_->setElement(row,index,0.);
									std::cout << "getElement("<<protein_accs[pa]<<","<<index<<")="
														<< model_->getElement(row,index) << "\n";
								}
							++f_v_idx2;
						}
				}
//       // show updated constraints:
//       std::set<String>::iterator citer = updated_constraints.begin();
//       Size num_feat_variables = variable_indices.size();
//       for(;citer != updated_constraints.end();++citer)
//         {
//           std::cout << *citer<<":\n";
//           Int row= cmodel_->row((String("PROT_COV_")+protein_accs[pa]).c_str());
//           for(Size idx = 0; idx < num_feat_variables;++idx)
//             {
//               if(cmodel_->getElement(row,idx) != 0.)
//                 {
//                   std::cout << cmodel_->getElement(row,idx) << " "
//                 }
//             }
            
//         }
      
//       cmodel_->writeMps((String("/home/zerck/data_project_maldi/results/ilp_ips/mit_prot_prophet/brukerdata/20100414/mps/iteration_") + String(iteration) + String(".mps")).c_str());
		}
  std::cout << "updated Combined ILP\n";
  timer.stop();
  std::cout <<timer.getClockTime()<<" seconds needed to update combined ILP.\n";

}


void PSLPFormulation::updateObjFunction_(String acc,FeatureMap<>& features,
                                         PrecursorIonSelectionPreprocessing& preprocessed_db,
                                         std::vector<IndexTriple>& variable_indices)
{
  std::cout << "Update Obj. function of combined ILP."<<std::endl;
	DoubleReal min_pt = param_.getValue("thresholds:min_pt_weight");
	DoubleReal min_rt_weight = param_.getValue("thresholds:min_rt_weight");
	DoubleReal mz_tolerance = param_.getValue("mz_tolerance");
  DoubleReal log_weight = param_.getValue("combined_ilp:k3");
  bool use_detectability = true;//param_.getValue("use_detectability") == "true" ? true : false;
  std::cout << "k3: "<<log_weight<<std::endl;
	std::cout << "parsed all parameters"<<std::endl;
  std::vector<IndexTriple> variable_indices_copy = variable_indices;
	sort(variable_indices_copy.begin(),variable_indices_copy.end(),VariableIndexLess());

	
  // if protein prob exceeds min threshold (TODO: is this necessary? we only use this ILP if either a protein is
  // identified or no more features can be found for this certain protein)
  // 					if(graph.getProteinProbability(accs[prot]) > min_protein_probability)
  // 						{
  std::map<String,std::vector<DoubleReal> >::const_iterator map_iter = preprocessed_db.getProteinPTMap().
    find(acc);
  if(map_iter !=  preprocessed_db.getProteinPTMap().end())
    {
      // then go through all (tryptic) peptides of this protein
      std::cout <<"protein "<<map_iter->first << " with "<< map_iter->second.size()
                << " entries in pt_map"<<std::endl;
      const std::vector<DoubleReal>& masses = preprocessed_db.getMasses(map_iter->first);
      for(Size p = 0;p < map_iter->second.size();++p)
        {
          // 													std::cout <<preprocessed_db.getRT(map_iter->first,p)
          // 																		<<" "<< masses[p]<< " "<<map_iter->second[p]
          // 																		<<std::endl;
          //std::cout << map_iter->second[p] << " > "<<min_pt<<std::endl;
          if(map_iter->second[p] > min_pt)
            {
              // go through all features
              for(Size f = 0; f < features.size();++f)
                {
                  if(features[f].getMetaValue("fragmented")=="true") continue;
                  // 			std::cout << masses[p] << " - "<< features[f].getMZ() << " -> "
                  // 																		<< fabs(masses[p] - features[f].getMZ())/masses[p] * 1e06<<" < "
                  // 																		<< mz_tolerance<<std::endl;
                  // and check if they match the current peptide
                  if(fabs(masses[p] - features[f].getMZ())/masses[p] * 1e06 < mz_tolerance)
                    {
                      DoubleReal rt_weight = preprocessed_db.getRTProbability(map_iter->first,p,features[f]);
                      if(rt_weight > min_rt_weight)
                        {
                          
                          std::cout << features[f].getMZ() << " "<<features[f].getRT()<<" "
                                    << rt_weight <<" is matching peptide "<<p <<"\n";
                          
													
                          // if yes: add as variable to the protein coverage constraint
                          // find corresponding variables
                          Size f_v_idx = 0;
                          while(f_v_idx < variable_indices_copy.size() && f != variable_indices_copy[f_v_idx].feature)
                            {
                              ++f_v_idx;
                            }
                          if(f_v_idx == variable_indices_copy.size())
                            {
                              std::cout << features[f].getMZ() << " "<< features[f].getRT()<<" "
                                        << /*rt_weight <<*/" is matching peptide "<< p
                                        << ", but not existing in variable indices???"
                                        << "--->This should not happen!"<<std::endl;
                            }
                          else
                            {
                              // TODO: enter here really all scans for this feature?
                              // or only the ones with rt-weight > min_rt_weight
                              while(f_v_idx < variable_indices_copy.size() && f == variable_indices_copy[f_v_idx].feature)
                                {
                                  if(model_->getObjective(f_v_idx) < 0.00000001)
                                    {
                                      ++f_v_idx;
                                      continue;
                                    }
                                  DoubleReal dt;
                                  if(use_detectability)
                                    {
                                      dt = map_iter->second[p];
                                    }
                                  else dt = 1.;
                                  // weight is detectability * rt_weight
                                  DoubleReal weight = dt * rt_weight;//preprocessed_db.getRTDTProbability(dt*curr_rt_weight);
                                  DoubleReal obj = model_->getObjective(f_v_idx);
                                  
                                  std::cout << features[f].getMZ() << " "
                                            << features[f].getRT()<<" "
                                            << rt_weight <<" is matching peptide "<<p <<" with score: ";
                                  std::cout << dt << " * " << rt_weight
                                            << " = "<< dt*rt_weight
                                            << " -> "<< log(1.-weight)<<"*"<<log_weight
                                            << " obj: "<<model_->getObjective(f_v_idx);
                                  // if(fabs(1.-weight)<0.000001)
                                  //                                         {
                                  //                                           if(obj+log_weight * log(0.000001) < 0.)
                                  //                                             {
                                  //                                               cmodel_->setObjective(f_v_idx,0.001);
                                  //                                             }
                                  //                                           else cmodel_->setObjective(f_v_idx,obj+log_weight * log(0.000001));
                                  //                                           //    std::cout <<"setting obj to "<<cmodel_->objective(f_v_idx)<<std::endl;
                                  //                                         }
                                  //                                       else if(obj+log_weight * log(1.-weight) < 0.)
                                  //                                         {
                                  //                                           cmodel_->setObjective(f_v_idx,0.001);
                                  //                                         }
                                  //                                       else cmodel_->setObjective(f_v_idx,obj+log_weight*(log(1.-weight)));
                                  if(log_weight * weight >  obj && obj > 0.)
                                    {
                                      model_->setObjective(f_v_idx,0.001);
                                    }
                                  else model_->setObjective(f_v_idx, obj - weight*log_weight);
                                  
                                  std::cout << " -> "<<model_->getObjective(f_v_idx)
                                            << " columnname: "<<model_->getColumnName(f_v_idx)
                                            << std::endl;
                                  
                                  ++f_v_idx;
                                }// while(f_v_idx < variable_indices_copy.size() && f == variable_indices_copy[f_v_idx].feature)
                            }// else if(f_v_idx == variable_indices.size())
                        }// if(rt_weight > min_rt_weight)
                    }// if(fabs(masses[p] - features[f].getMZ())/masses[p] * 1e06 < mz_tolerance)
                }//for(Size f = 0; f < features.size();++f)
            }//if(map_iter->second[p] > min_pt)
        }//for(Size p = 0;p < map_iter->second.size();++p)
    }//if(map_iter !=  preprocessed_db.getProteinPTMap().end())
  std::cout << "Update Obj. function of combined ILP -->done"<<std::endl;	
}

  

} // namespace OpenMS
