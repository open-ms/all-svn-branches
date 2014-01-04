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


#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>

#include <OpenMS/ANALYSIS/DENOVO/AntilopeSpectrumGraph.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringBayes.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringRank.h>

#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#undef Debug
//#define Debug

namespace OpenMS
{
  const DoubleReal SpectrumGraphSeqan::INFINITYdist = seqan::_getInfinityDistance<DoubleReal>();
  
  
  
  
  //-------------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------constructors----------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  
  // default constructor
  SpectrumGraphSeqan::SpectrumGraphSeqan() :
  DefaultParamHandler("SpectrumGraphSeqan"),
  graph()
  {
    defaults_.setValue("tryptic", "true", "whether the query spectrum is from peptides digested with trypsin");
    std::vector<String>valid_tryptic;
    valid_tryptic.push_back("true");
    valid_tryptic.push_back("false");
    defaults_.setValidStrings("tryptic", valid_tryptic);
    
    defaults_.setValue("precision", 0.1, "The m/z precision");
    defaults_.setMaxFloat("precision", 1.0);
    defaults_.setMinFloat("precision", 0.01);
    defaults_.setValue("delta", 0.5, "the allowed mass error for each peak");
    defaults_.setMaxFloat("delta", 1.0);
    defaults_.setMinFloat("delta", 0.1);
    
    defaultsToParam_();
  }
  
  //copy constructor
  SpectrumGraphSeqan::SpectrumGraphSeqan(const SpectrumGraphSeqan & in) :
  DefaultParamHandler(in),
  graph(in.graph),
  conflicting_nodes(in.conflicting_nodes),
  pairwise_conflicting(in.pairwise_conflicting),
  conflict_clusters(in.conflict_clusters),
  conflict_clusters_inv(in.conflict_clusters_inv),
  vertices_(in.vertices_),
  edges_(in.edges_),
  ordering_(in.ordering_),
  edge_weights_(in.edge_weights_)
  {
  }
  
  //assignment operator
  SpectrumGraphSeqan& SpectrumGraphSeqan::operator=(const SpectrumGraphSeqan& in)
  {
    if (this != &in)
    {
      DefaultParamHandler::operator=(in);
      conflicting_nodes = in.conflicting_nodes;
      pairwise_conflicting = in.pairwise_conflicting;
      conflict_clusters = in.conflict_clusters;
      conflict_clusters_inv = in.conflict_clusters_inv;
      vertices_ = in.vertices_;
      edges_ = in.edges_;
      ordering_ = in.ordering_;
      graph = in.graph;
      edge_weights_ = in.edge_weights_;
    }
    return *this;
  }
  
  
  //-------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------edge_weight--------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  // return weight of edge (i,j)
  DoubleReal SpectrumGraphSeqan::getEdgeWeight(VertexDescriptor idx_l, VertexDescriptor idx_r) const
  {
    EdgeDescriptor ed = seqan::findEdge(graph, idx_l, idx_r);
    if (!ed)
    {
      return INFINITYdist;
    }
    else
    {
      return seqan::getProperty(edge_weights_, ed);
    }
  }
  
  // return weight of edge (i,j)
  DoubleReal SpectrumGraphSeqan::getEdgeWeight(EdgeDescriptor ed) const
  {
      return seqan::getProperty(edge_weights_, ed);
  }
  
  
  //-------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------getEdgeWeights-----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  
  /// return all edge weights (i,j)
  void SpectrumGraphSeqan::getEdgeWeights(seqan::String<DoubleReal> &weights) const
  {
    weights = edge_weights_;
  }
  
  //-------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------setEdgeWeights-----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  
  /// set all edge weights (i,j)
  void SpectrumGraphSeqan::setEdgeWeights()
  {
    seqan::resizeEdgeMap(graph, edge_weights_);
    EdgeIterator it(graph);
    
    while(!seqan::atEnd(it))
    {
      VertexDescriptor source = seqan::sourceVertex(graph, *it);
      DoubleReal factor = 1.0;
      if (seqan::getProperty(edges_, *it).length >= 2)
      {
        factor = 0.5;
      }
      DoubleReal score = (-1) * seqan::getProperty(vertices_, source).score * factor;
      
      seqan::property(edge_weights_, *it) = score;
      ++it;
    }
  }
  
  
  //-------------------------------------------------------------------------------------------------------------------------
  //-----------------------------------------------has_directed_edge---------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  // return whether to nodes i and j are connected via a directed edge
  SpectrumGraphSeqan::EdgeDescriptor SpectrumGraphSeqan::findDirectedEdge(VertexDescriptor i, VertexDescriptor j) const
  {
    return seqan::findEdge(graph, i, j);
  }
  
  //-------------------------------------------------------------------------------------------------------------------------
  //----------------------------------------------has_undirected_edge--------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  // return wheter to nodes i and j are connected via a undirected edge i.e are complementary
  bool SpectrumGraphSeqan::hasUndirectedEdge(VertexDescriptor i, VertexDescriptor j) const
  {
    if(i < j)
      return pairwise_conflicting[i][j];
    else
      return pairwise_conflicting[j][i];
  }
  
  //-------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------create_node_set----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  /// creating the node set by creating a specGraph G and adding the nodes. actually for each peak two nodes are created
  // a b and a y ion. later on masses  appearing twice are joined to one node.
  UInt SpectrumGraphSeqan::createNodeSet(const PeakSpectrum & spectrum, const IdSetup::BoolVec &A, const IdSetup::UIntVec &A_len)
  {
    
    //parameter extraction
    DoubleReal delta = param_.getValue("delta");
    DoubleReal precision = param_.getValue("precision");
    String tryptic_flag = param_.getValue("tryptic");
    bool tryptic = tryptic_flag == "true" ? true : false;
    
    const DoubleReal charge = spectrum.getPrecursors()[0].getCharge();
    const DoubleReal parent_mass = (spectrum.getPrecursors()[0].getPosition()[0]*charge) - (charge-1);
    const DoubleReal parent_peptide_mass = parent_mass - 19;
    
    std::cout << "parentmass: " << parent_mass << std::endl;
    
    //Filter the Spectrum. Do not create nodes for every peak, only for those that are likely to correspond
    //to a BION or a YION
    
    //local copy of the spectrum
    PeakSpectrum tmp_spectrum(spectrum);
    
    //create one node for each considered ion type (BION and YION)
    std::vector<Node> nodes;
    
    nodes.push_back(Node(0.0, 0, 0, -1, 1)); // artificial source node
    
    UInt integer_MZ;
    DoubleReal real_MZ;
    
    //preprocess spectrum
    IdSetup::filterSpectrumInspect(tmp_spectrum);// preproc;
    Normalizer().filterPeakSpectrum(tmp_spectrum);  //TODO remove
    
    //compute the rank scores
    DVector b_rank_scores, y_rank_scores;
    RankScoringFunction rnk_scf;
    
    rnk_scf.loadModel();
    rnk_scf.getRankScores(b_rank_scores,y_rank_scores,tmp_spectrum);
    
    DoubleReal max_node_mass = parent_peptide_mass - 56.; //TODO exact mass of minimum amino acid
    
    UInt index = 0;
    for (PeakSpectrum::ConstIterator it = tmp_spectrum.begin(); it != tmp_spectrum.end(); ++it)
    {
      if(it->getMZ() < parent_peptide_mass)
      {
        //add node for B-Ions
        //but no B1 ion peak
        real_MZ = IdSetup::getPrefixMass(IdSetup::BIon, it->getMZ(), parent_mass);
        integer_MZ = (UInt) floor(real_MZ / precision + 0.5);
        
        if (real_MZ < max_node_mass && (integer_MZ >= A.size() || (A[integer_MZ] && (A_len[integer_MZ] > 1))))
        {
          nodes.push_back(Node(real_MZ, integer_MZ, b_rank_scores[index], index, IdSetup::BIon, it->getIntensity()));//b-ion
          //nodes.push_back(Node(real_MZ, integer_MZ, it->getIntensity(), index, IdSetup::BIon, it->getIntensity()));//b-ion
        }
        
        //add node for Y-Ions
        real_MZ = IdSetup::getPrefixMass(IdSetup::YIon, it->getMZ(), parent_mass);
        integer_MZ = (UInt) floor(real_MZ / precision + 0.5);
        if (real_MZ < max_node_mass && (integer_MZ>=A.size() || A[integer_MZ]))
        {
//          nodes.push_back(Node(real_MZ, integer_MZ, y_rank_scores[index], index, IdSetup::YIon, it->getIntensity()));//y-ion
          nodes.push_back(Node(real_MZ, integer_MZ, it->getIntensity(), index, IdSetup::YIon, it->getIntensity()));//y-ion
        }
        
        //add node for Y2-Ions
        real_MZ = IdSetup::getPrefixMass(IdSetup::AIon, it->getMZ(), parent_mass);
        integer_MZ = (UInt) floor(real_MZ / precision + 0.5);
        
        if (real_MZ < max_node_mass && (integer_MZ >= A.size() || (A[integer_MZ] && (A_len[integer_MZ] > 1))))
        {
          nodes.push_back(Node(real_MZ, integer_MZ, b_rank_scores[index], index, IdSetup::BIon, it->getIntensity()));//b-ion
//          nodes.push_back(Node(real_MZ, integer_MZ, it->getIntensity(), index, IdSetup::AIon, it->getIntensity()));//b-ion
        }
        
        //add node for B-NH3-Ions
        real_MZ = IdSetup::getPrefixMass(IdSetup::BIon_nh3, it->getMZ(), parent_mass);
        integer_MZ = (UInt) floor(real_MZ / precision + 0.5);
        
        if (real_MZ < max_node_mass && (integer_MZ >= A.size() || (A[integer_MZ] && (A_len[integer_MZ] > 1))))
        {
          //          nodes.push_back(Node(real_MZ, integer_MZ, b_rank_scores[index], index, IdSetup::BIon, it->getIntensity()));//b-ion
          nodes.push_back(Node(real_MZ, integer_MZ, it->getIntensity(), index, IdSetup::BIon_nh3, it->getIntensity()));//b-ion
        }

        //add node for Y-H2O-Ions
        real_MZ = IdSetup::getPrefixMass(IdSetup::YIon_h2o, it->getMZ(), parent_mass);
        integer_MZ = (UInt) floor(real_MZ / precision + 0.5);
        
        if (real_MZ < max_node_mass && (integer_MZ >= A.size() || (A[integer_MZ])))
        {
          //          nodes.push_back(Node(real_MZ, integer_MZ, b_rank_scores[index], index, IdSetup::BIon, it->getIntensity()));//b-ion
          nodes.push_back(Node(real_MZ, integer_MZ, it->getIntensity(), index, IdSetup::YIon_h2o, it->getIntensity()));//b-ion
        }

      }
      ++index;
    }
    
    //if the data are tryptic digested peptides, the last amino acid is supposed to be a K or R
    //hence the corresponding nodes are added
    if (tryptic)
    {
      real_MZ = parent_peptide_mass - 128.09496;
      integer_MZ = floor(real_MZ / precision + 0.5);
      nodes.push_back(Node(real_MZ, integer_MZ, 0, -2, IdSetup::BIon));//K
      
      real_MZ = parent_peptide_mass - 156.10111;
      integer_MZ = floor(real_MZ / precision + 0.5);
      
      nodes.push_back(Node(real_MZ, integer_MZ, 0, -2, IdSetup::BIon));//R
    }
    
    // artificial sink node
    const UInt parent_peptide_mass_int = (UInt) floor(parent_peptide_mass / precision + 0.5);
    nodes.push_back(Node(parent_peptide_mass, parent_peptide_mass_int, 0, -1, 0));
    
    //sort the array of nodes by mass
    std::sort(nodes.begin(), nodes.end());
    
#ifdef Debug
    //test output
    for (UInt i = 0; i < nodes.size(); ++i)
    {
      std::cout << "0 -- node " << i << "  mass:  " << nodes[i].int_mass <<"  "<<nodes[i].real_mass<< "  interprets: ";
      std::set<Size>::const_iterator it;
      for (it = nodes[i].generating_peaks.begin(); it != nodes[i].generating_peaks.end(); ++it)
      {
        std::cout << " " << *it;
      }
      std::cout << "  has intensity:  " << nodes[i].score << std::endl;
    }
#endif
    
    std::cout << "numV before merging: " << nodes.size() << std::endl;
    
    //---------------------------------------------------------------
    //--------------------NODE MERGING-------------------------------
    
    //vector delete[i] true if node i deleted
    std::vector<bool> deleted(nodes.size(), false);
    
    DoubleReal min_dist = 0;
    
    DoubleReal last_mass, mass_diff;
    
    UInt ln, rn, last_node;
    
    //merging nodes that are within the allowed M/Z-delta
    while (min_dist < delta)
    {
      min_dist = (UInt) (delta + 1);
      last_node = 0;
      last_mass = 0;
      ln = 0;
      rn = 0;
      
      for (UInt i = 2; i < nodes.size(); ++i)
      {
        if (!deleted[i])
        {
          mass_diff = nodes[i].real_mass - last_mass;
          //if mass difference is smaller than all previous remember nodes
          if (mass_diff < min_dist)
          {
            min_dist = mass_diff;
            ln = last_node;
            rn = i;
          }
          
          last_node = i;
          last_mass = nodes[i].real_mass;
        }
      }
      
      //if the minimal distance is smaller than delta the nodes are joined
      if (min_dist < delta)
      {
        //debug info
#ifdef Debug
        std::cout << "joined nodes " << ln << ": " << nodes[ln].real_mass << " and " << rn << " : "
        << nodes[rn].real_mass << std::endl;
#endif
        
        nodes[ln].int_mass = (UInt) floor(0.5 * (nodes[ln].int_mass + nodes[rn].int_mass) + 0.5);
        nodes[ln].real_mass = 0.5 * (nodes[ln].real_mass + nodes[rn].real_mass);
        
        nodes[ln].score = (std::max(0.0,nodes[ln].score) + std::max(0.0,nodes[rn].score));
        
        nodes[ln].generating_peaks.insert(nodes[rn].generating_peaks.begin(), nodes[rn].generating_peaks.end());
        
        deleted[rn] = true;
      }
    }
#ifdef Debug
    //test output
    for (UInt i = 0; i < nodes.size(); ++i)
    {
      std::cout << "0 -- node " << i << "  mass:  " << nodes[i].int_mass <<"  "<<nodes[i].real_mass<< "  interprets: ";
      std::set<Size>::const_iterator it;
      for (it = nodes[i].generating_peaks.begin(); it != nodes[i].generating_peaks.end(); ++it)
      {
        std::cout << " " << *it;
      }
      std::cout << "  has intensity:  " << nodes[i].score << std::endl;
    }
#endif
    
    //in case of tryptically digested peptides remove all nodes with masses greater than of the two only possible last nodes
    if (tryptic)
    {
      UInt i = nodes.size() - 2;
      
      while (nodes[i].real_mass > parent_peptide_mass - 156.1 - delta)
      {
        if(fabs(nodes[i].real_mass - (parent_peptide_mass - 156.1) ) > delta &&
           fabs(nodes[i].real_mass - (parent_peptide_mass - 128.1) ) > delta)
        {
          deleted[i] = true;
          std::cout<<"delete node "<< i <<std::endl;
        }
        --i;
      }
    }//end of tryptic specific actions
    
    
    //remove deleted nodes from vector
    std::vector<Node> tmp_nodes(nodes);
    nodes.clear();
    
    Size pos = 0;
    for (std::vector<Node>::iterator it = tmp_nodes.begin(); it != tmp_nodes.end(); ++it, ++pos)
    {
      if (!deleted[pos])
      {
        nodes.push_back(*it);
      }
    }
    
    tmp_nodes.clear();
    
    vertices_.clear();
    vertices_ = nodes;
    
    seqan::resize(ordering_, nodes.size());
    for (Size i = 0; i < nodes.size(); ++i)
      ordering_[i] = i;
    
    
    std::cout << "numV after merging: " << vertices_.size() << std::endl;
    
#ifdef Debug
    //test output
    for (UInt i = 0; i < nodes.size(); ++i)
    {
      std::cout << "0 -- node " << i << "  mass:  " << nodes[i].int_mass <<"  "<<nodes[i].real_mass<< "  interprets: ";
      std::set<Size>::const_iterator it;
      for (it = nodes[i].generating_peaks.begin(); it != nodes[i].generating_peaks.end(); ++it)
      {
        std::cout << " " << *it;
      }
      std::cout << "  has intensity:  " << nodes[i].score << std::endl;
    }
#endif
    
    return 0;
  }
  
  //-------------------------------------------------------------------------------------------------------------------------
  //-----------------------------------------------scoreNodes-----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  UInt SpectrumGraphSeqan::scoreNodes(BayesScoring &scoring_func, const PeakSpectrum &spectrum, bool add, DoubleReal weight, DoubleReal threshold)
  {
    //this copy is used for the intensity score. flatened by square-routing before normalizing to one
    PeakSpectrum tmp_spectrum(spectrum);
    
    String tryptic_flag = param_.getValue("tryptic");
    bool tryptic = tryptic_flag == "true" ? true : false;
    
    DoubleReal parent_mass = tmp_spectrum.getPrecursors()[0].getPosition()[0];
    DoubleReal charge = tmp_spectrum.getPrecursors()[0].getCharge();
    parent_mass = (parent_mass*charge) - (charge-1);
    
    
    PeakSpectrum disc_spectrum(tmp_spectrum);
    scoring_func.normalizeIntensity(disc_spectrum);
    
    
    //the conditional probability score
    for (Size i = 1; i < vertices_.size()-1; ++i)
    {
      std::vector<IdSetup::interpretation> interprets;
      
      if(add)
      {
        DoubleReal tmp_score_old = vertices_[i].score;
        
        vertices_[i].score += weight * scoring_func.getLogScore(disc_spectrum, vertices_[i].real_mass, interprets);
#ifdef Debug
        std::cout<<"Node: "<<i<<" prev score: "<<tmp_score_old<<" New score: "<< vertices_[i].score << std::endl;
#endif
      }
      else
      {
        vertices_[i].score = scoring_func.getLogScore(disc_spectrum, vertices_[i].real_mass, interprets);
      }
      
      //recompute the peak indices of the original spectrum
      for(std::vector<IdSetup::interpretation>::iterator peak_iter = interprets.begin(); peak_iter != interprets.end(); ++peak_iter)
      {
        peak_iter->peak_nr = spectrum.findNearest(disc_spectrum[peak_iter->peak_nr].getMZ());
      }
      
      vertices_[i].peak_interprets = interprets;
    }
    
    //filter nodes below threshold
    std::vector<Node>tmp_vertices = vertices_;
    vertices_.clear();
    vertices_.push_back(tmp_vertices[0]);
    
    for(Size i = 1; i < tmp_vertices.size(); ++i)
    {
      if(tmp_vertices[i].score >= threshold ||
         (tryptic && i >= tmp_vertices.size() - 3))
      {
        vertices_.push_back(tmp_vertices[i]);
        std::cout << "v " << vertices_.size() << ": " << tmp_vertices[i].score << " " << tmp_vertices[i].real_mass << std::endl;
      }
    }
    
    seqan::resize(ordering_, vertices_.size());
    
    return 0;
  }
  
  
  //-------------------------------------------------------------------------------------------------------------------------
  //-----------------------------------------------create_edge_set-----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  ///create the directed edges connecting nodes with a mass difference equal +/- precision to the mass of a single or several amino acids
  UInt SpectrumGraphSeqan::createEdgeSet(const IdSetup::BoolVec &A, const IdSetup::UIntVec &A_len)
  {
    //TODO as parameter!
    const UInt max_edge_len = 3;
    const Size number_of_nodes = vertices_.size();
    
    //generate vertex order -- required for efficient dag-shortest-path computation in seqan
    seqan::clearVertices(graph);
    seqan::String<VertexDescriptor> vertex_order;
    seqan::resize(vertex_order, number_of_nodes);
    
    for(Size i = 0; i < number_of_nodes; ++i)
    {
      seqan::addVertex(graph);
      seqan::assignValue(vertex_order, i, i);
    }
    
    
    //list of all nodes that have a mass difference to the current node which equals one or several amino acid masses. is a candidate for a new edge if not redundant
    std::vector<VertexDescriptor> target_candidates;
    
    String tryptic_flag = param_.getValue("tryptic");
    bool tryptic = tryptic_flag == "true" ? true : false;
    
    UInt mass_diff_index;
    
    seqan::String <VertexDescriptor> pred_map;
    seqan::String <DoubleReal> dist_map;
    seqan::String <DoubleReal> weight_map;
    
    Size source_node, target_node, nodes_end;
    nodes_end = vertices_.size();
    
    //if tryptic the edges from the last two nodes to the sink must be added now
    if (tryptic)
    {
      Node &last_v = seqan::property(vertices_, vertices_.size() - 1);
      Node &tryp_v = seqan::property(vertices_, vertices_.size() - 2);
      Node &tryp_v2 = seqan::property(vertices_, vertices_.size() - 3);
      
      Edge tmp_edge_a = {(last_v.real_mass - tryp_v.real_mass), last_v.int_mass - tryp_v.int_mass, 1};
      edges_.push_back(tmp_edge_a);
      seqan::addEdge(graph, vertices_.size() - 2, vertices_.size() - 1);
      seqan::appendValue(weight_map, 1);
      
      Edge tmp_edge_b = {(last_v.real_mass - tryp_v2.real_mass), last_v.int_mass - tryp_v2.int_mass, 1};
      edges_.push_back(tmp_edge_b);
      seqan::addEdge(graph, vertices_.size() - 3, vertices_.size() - 1);
      seqan::appendValue(weight_map, 1);
      
      --nodes_end;
    }
    
    // create edge for each pair of nodes fullfilling the constraints, please read Chen paper pp 392 for understanding the role of the array A_
    for (UInt len = 1; len <= max_edge_len; len++)
    {
      for (source_node = 0; source_node < nodes_end; ++source_node)
      {
        for(target_node = source_node + 1; target_node < nodes_end; ++target_node)
        {
          //massDiffIndex is the mass difference between nodes j and i (integer mass difference)
          mass_diff_index = seqan::property(vertices_, target_node).int_mass - seqan::property(vertices_, source_node).int_mass;
          
          if(mass_diff_index >= A.size())
          {
            //if mass_diff exceeds range we can terminate loop
            break;
          }
          
          if (A[mass_diff_index] && A_len[mass_diff_index] == len)
          {
            target_candidates.push_back(target_node);
          }
        }
        
        Node& source = seqan::property(vertices_, source_node);
        
        if (len == 1)
        {
          for (std::vector<VertexDescriptor>::iterator target_cand_it = target_candidates.begin(); target_cand_it
               != target_candidates.end(); ++target_cand_it)
          {
            if(!pairwise_conflicting[source_node][*target_cand_it])
            {
              Node& target = seqan::property(vertices_, *target_cand_it);
              Edge tmp_edge;
              tmp_edge.length = len;
              tmp_edge.real_mass = target.real_mass - source.real_mass;
              tmp_edge.int_mass = target.int_mass - source.int_mass;
              edges_.push_back(tmp_edge);
              seqan::addEdge(graph, (VertexDescriptor)source_node, *target_cand_it);
              seqan::appendValue(weight_map, 1);
#ifdef Debug
              std::cout << "added single AA edge  " << source_node << " --> " << *target_cand_it << std::endl;
#endif
            }
          }
        }
        else if(target_candidates.size() > 0)
        {
          //compute shortest path from node i to the last candidate and simlultaniously also to all other preceding candidates
          //for nodes that are already reachable no edge is added as it would be redundant
#ifdef Debug
          std::cerr<<"weight map legnth:  " << seqan::length(weight_map)<<std::endl;
#endif
          dagShortestPathST(graph, (VertexDescriptor)source_node, target_candidates.back(), weight_map, pred_map, dist_map, vertex_order);
          
          //for all nodes not reached we can add the edge to the graph. it is not redundant
          for (std::vector<VertexDescriptor>::iterator target_cand_it = target_candidates.begin(); target_cand_it
               != target_candidates.end(); ++target_cand_it)
          {
            if (pred_map[*target_cand_it] == seqan::getNil<VertexDescriptor>() && !pairwise_conflicting[source_node][*target_cand_it])
            {
              Node& target = seqan::property(vertices_, *target_cand_it);
              Edge tmp_edge;
              tmp_edge.length = len;
              tmp_edge.real_mass = target.real_mass - source.real_mass;
              tmp_edge.int_mass = target.int_mass - source.int_mass;
              edges_.push_back(tmp_edge);
              seqan::addEdge(graph, (VertexDescriptor)source_node, *target_cand_it);
              seqan::appendValue(weight_map, 1);
#ifdef Debug
              std::cout << "added " << len << " AA edge  " << source_node << " --> " << *target_cand_it << std::endl;
#endif
            }
#ifdef Debug
            //test_output
            /*
             else
             {
             std::cout << "REDUNDANT: " << source_node << " --> " << *target_cand_it << std::endl;
             std::cout << *target_cand_it << "  " <<pred_map[*target_cand_it] <<std::endl;
             }
             */
#endif
          }
        }
        target_candidates.clear();
      }
    }
    
    seqan::resize(vertex_active_, vertices_.size(), true);
    seqan::resize(edge_active_, edges_.size(), true);
    return 0;
  }
  
  
  //-------------------------------------------------------------------------------------------------------------------------
  //-------------------------------------------create_undirected_edges-------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  
  //Rather important is the inciden list for each node, so that in the lagrangian relaxation we can use the index
  //of each undirected edge as identifier for the subgradients.
  //together with the vector of all undirected edges with the information: terminal nodes + id all useful information is accessible
  UInt SpectrumGraphSeqan::createUdirEdgeSet()
  {
    //create the set of conflict clusters (hyperedges)
    std::vector<Size> mapping;
    Size nextId = 0;
    
    for (Size i = 1; i < vertices_.size() - 1; ++i)
    {
      std::set<Size>::iterator gen_it;
      for (gen_it = vertices_[i].generating_peaks.begin(); gen_it != vertices_[i].generating_peaks.end(); ++gen_it)
      {
        if (*gen_it >= mapping.capacity())
          mapping.reserve(2 * *gen_it);
        if (*gen_it >= mapping.size())
          mapping.resize(*gen_it + 1, -1u);

        if (mapping[*gen_it] == -1u)
        {
          mapping[*gen_it] = nextId;
          ++nextId;
        }


        Size id = mapping[*gen_it];
        assert(id <= conflict_clusters.size());

        if (id == conflict_clusters.size())
          conflict_clusters.push_back(std::vector<VertexDescriptor>(1,i));
        else
          conflict_clusters[id].push_back(i);
      }
    }
   
    //remove singleton or duplicated clusters
    Size new_size = 0;
    for (Size i = 0; i < conflict_clusters.size(); ++i)
    {
      if (conflict_clusters[i].size() >= 2)
      {
        conflict_clusters[i].swap(conflict_clusters[new_size]);
        ++new_size;
      }
    }
    conflict_clusters.resize(new_size);
    std::sort(conflict_clusters.begin(), conflict_clusters.end());
    conflict_clusters.erase(std::unique(conflict_clusters.begin(), conflict_clusters.end()),conflict_clusters.end());
    
    
    //create the inverse table --> for each node all clusters it appears in
    conflict_clusters_inv.assign(vertices_.size(), std::vector<Size>());
    conflicting_nodes.assign(vertices_.size(), std::vector<VertexDescriptor>());
    pairwise_conflicting.assign(vertices_.size(), std::vector<bool>(vertices_.size(), false));
    for (Size i = 0; i < conflict_clusters.size(); ++i)
    {
      for (Size j = 0; j < conflict_clusters[i].size(); ++j)
      {
        conflict_clusters_inv[conflict_clusters[i][j]].push_back(i);
        for (Size k = j+1; k < conflict_clusters[i].size(); ++k)
        {
          conflicting_nodes[conflict_clusters[i][j]].push_back(conflict_clusters[i][k]);
          conflicting_nodes[conflict_clusters[i][k]].push_back(conflict_clusters[i][j]);
          
          pairwise_conflicting[conflict_clusters[i][j]][conflict_clusters[i][k]] = true;
          pairwise_conflicting[conflict_clusters[i][k]][conflict_clusters[i][j]] = true;
        }
      }
    }
    
    //remove duplication entries in conflicting_nodes
    for (Size i = 0; i < vertices_.size(); ++i)
    {
      std::vector<VertexDescriptor> &cni = conflicting_nodes[i];
      std::sort(cni.begin(), cni.end());
      cni.erase(std::unique(cni.begin(), cni.end()), cni.end());
    }
    
//    DEBUG OUTPUT
//    for (Size i = 0; i < conflict_clusters.size(); ++i)
//    {
//      std::cout << "cluster " << i << ": ";
//      for (Size j = 0; j < conflict_clusters[i].size(); ++j)
//      {
//        std::cout << conflict_clusters[i][j] << '\t';
//      }
//      std::cout << std::endl;
//    }
//    
//    for (Size i = 0; i < conflict_clusters_inv.size(); ++i)
//    {
//      std::cout << "vertex " << i << "contained in: ";
//      for (Size j = 0; j < conflict_clusters_inv[i].size(); ++j)
//      {
//        std::cout << conflict_clusters_inv[i][j] << '\t';
//      }
//      std::cout << std::endl;
//    }
    
    return 0;
  }
  
  
  
  
  
  
  
  
  
  //    //clear the previous entries
  //    for (Size i = 0; i < conflict_edge_ids_.size(); ++i)
  //    {
  //      conflict_edge_ids_.clear();
  //    }
  //
  //    //create an index of conflicting interpretations
  //    std::vector<std::vector<Size> >conflict_index(spectrum.size());
  //
  //    for (Size i = 1; i < vertices_.size() - 1; ++i)
  //    {
  //      std::set<Size>::iterator gen_it;
  //      for (gen_it = seqan::property(vertices_, i).generating_peaks.begin(); gen_it != seqan::getProperty(vertices_, i).generating_peaks.end(); ++gen_it)
  //      {
  //        conflict_index[*gen_it].push_back(i);
  //      }
  //    }
  //
  //    UInt edge_id = 0;
  //    for (Size i = 0; i < conflict_index.size(); ++i)
  //    {
  //      std::vector<Size>::iterator it1, it2;
  //
  //      for(it1 = conflict_index[i].begin(); it1 != conflict_index[i].end(); ++it1)
  //      {
  //        for(it2 = it1 + 1; it2 < conflict_index[i].end(); ++it2)
  //        {
  //          if(!conflict_edge_matrix_[*it1][*it2] && (seqan::property(vertices_, *it2).real_mass - seqan::property(vertices_, *it1).real_mass > 57.0) )
  //          {
  //            ConflictEdge c_edge;
  //            c_edge.left_node = *it1;
  //            c_edge.right_node = *it2;
  //            conflict_edges_.push_back(c_edge);
  //
  //            conflict_edge_ids_[*it1].push_back(edge_id);
  //            conflict_edge_ids_[*it2].push_back(edge_id);
  //
  //            conflicting_nodes_[*it1].push_back(*it2);
  //            conflicting_nodes_[*it2].push_back(*it1);
  //            conflict_edge_matrix_[*it1][*it2] = true;
  //            conflict_edge_matrix_[*it2][*it1] = true;
  //#ifdef Debug
  //            std::cout<<"add undirected edge: "<<*it1<<"  <--> "<<*it2<<std::endl;
  //#endif
  //
  //            ++edge_id;
  //          }
  //        }
  //      }
  //    }
  //
  //    return 0;
  //  }
  
  //  void SpectrumGraphSeqan::createSpanningEdges()
  //  {
  //    EdgeIterator edge_it(graph);
  //    spanning_edges_.assign(size(), std::vector<EdgeDescriptor>());
  //    
  //    while(!seqan::atEnd(edge_it))
  //    {
  //      VertexDescriptor source = seqan::sourceVertex(graph, *edge_it);
  //      VertexDescriptor target = seqan::targetVertex(graph, *edge_it);
  //      for(Size j = source + 1; j < target; ++j)
  //      {
  //        spanning_edges_[j].push_back(*edge_it);
  //      }
  //      ++edge_it;
  //    }
  //    return;
  //  }
  
  }//namespace OpenMS
  
