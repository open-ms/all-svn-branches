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

#ifndef SPECTRUMGRAPH_SEQAN_H
#define SPECTRUMGRAPH_SEQAN_H

// general stuff
#include<vector>
#include<map>
#include<list>

//denovo includes
#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeIonScoringBayes.h>
//#include <Common/CommonTypes.h>

//OpenMS includes
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>


#include <seqan/graph_algorithms.h>
#include <seqan/graph_types.h>


namespace OpenMS
{

class SpectrumGraphSeqan: public DefaultParamHandler
{

public:

  typedef std::vector<DoubleReal> DVector;

  typedef IdSetup::interpretation interpretation;
  typedef std::vector<interpretation>InterpretVec;

  typedef seqan::Graph<seqan::Directed<> > TGraph;
  typedef seqan::VertexDescriptor<TGraph>::Type VertexDescriptor;
  typedef seqan::EdgeDescriptor<TGraph>::Type EdgeDescriptor;

  typedef seqan::Iterator<TGraph, seqan::VertexIterator>::Type VertexIterator;
  typedef seqan::Iterator<TGraph, seqan::EdgeIterator>::Type EdgeIterator;
  typedef seqan::Iterator<TGraph, seqan::OutEdgeIterator>::Type OutEdgeIterator;

  typedef seqan::Size<TGraph>::Type GSize;
  const static DoubleReal INFINITYdist;


  struct Edge
  {
    DoubleReal real_mass;
    Size int_mass;
    Size length;
  };

  struct ConflictEdge
  {
    VertexDescriptor left_node;
    VertexDescriptor right_node;
  };


  class Node
  {

  public:

    /// floating point mass of the node
    DoubleReal real_mass;

    /// integer mass
    UInt int_mass;

    ///the peak interpretations of this node
    std::vector<interpretation> peak_interprets;

    ///the generating peaks of this nodes (more than one if merged)
    std::set<Size> generating_peaks;

    /// the score of the node representing its reliability -- needed for edge weights
    DoubleReal score;

    ///the intensity of the peaks interpreted by the node
    DoubleReal intensity;

    bool operator<(const Node& other) const
    {
      return int_mass < other.int_mass;
    }

    ///constructor
    Node(DoubleReal rmass, UInt imass, DoubleReal inscore = -1, Int peak_id = -1, UInt type = 0, DoubleReal intens = 0) :
      real_mass(rmass),
      int_mass(imass),
			peak_interprets(std::vector<interpretation>(0)),
      score(inscore),
      intensity(intens)      
    {
      if(peak_id >= 0)
      generating_peaks.insert(peak_id);
      std::cout << "created node: " << rmass << "  " << imass << "  Score-> " << inscore << " " << peak_id << " "
          << type << std::endl;
    }

    ///default Constructor
    Node() :
      real_mass(0),
      int_mass(0),
      peak_interprets(std::vector<interpretation>(0)),
      generating_peaks(std::set<Size>()),
      score(0),
      intensity(0)
    {
    }

    ///destructor
    ~Node()
    {      
    }

    ///copy Constructor
    Node(const Node & node) :
      real_mass(node.real_mass),
      int_mass(node.int_mass),
      peak_interprets(node.peak_interprets),
      generating_peaks(node.generating_peaks),
      score(node.score),
      intensity(node.intensity)
    {
    }

    ///assignment operator
    Node &operator=(const Node & node)
    {
      if (this == &node)
        return *this;
      else
      {
        real_mass = node.real_mass;
        int_mass = node.int_mass;
        peak_interprets = node.peak_interprets;
        generating_peaks = node.generating_peaks;
        score = node.score;
        intensity = node.intensity;
      }
      return *this;
    }
  }; //end class Node


//public members
public:

  //The Seqan Graph
  TGraph graph;

  /// default constructor
  SpectrumGraphSeqan();

  /// copy constructor
  SpectrumGraphSeqan(const SpectrumGraphSeqan & in);

  /// assignment operator
  SpectrumGraphSeqan& operator=(const SpectrumGraphSeqan& in);

  /// return score of edge (i,j)
  DoubleReal getEdgeWeight(VertexDescriptor, VertexDescriptor) const;

  /// return all edge weights (i,j)
  void getEdgeWeights(seqan::String<DoubleReal> &weights) const;

  /// return whether two nodes i and j are connected via undirected edge --> are complementary
  bool hasUndirectedEdge(VertexDescriptor, VertexDescriptor) const;

  /// return whether two nodes i and j are connected via directed edge
  EdgeDescriptor findDirectedEdge(VertexDescriptor, VertexDescriptor) const;

  /// adding a node
  UInt addNodeToGraph(DoubleReal, int, DoubleReal, DoubleReal = -1);

  /// method to construct the spectrum graph (generate set of edges by editing adjacency lists)
  UInt createEdgeSet(const IdSetup::BoolVec &A, const IdSetup::UIntVec &A_len);

  /// construct the set of undirected edges
  UInt createUdirEdgeSet(const PeakSpectrum &spectrum, UInt mode = 1);

  /// method to construct the node set
  UInt createNodeSet(const PeakSpectrum & spectrum, const IdSetup::BoolVec &A, const IdSetup::UIntVec &A_len);  

  /// returns the number of nodes in the graph
  Size size() const
  {
    return vertices_.size();
  }

  /// return the real mass of a node
  DoubleReal getRealMass(VertexDescriptor vertex_idx)
  {
    return seqan::property(vertices_, vertex_idx).real_mass;
  }

  /// return the score of a node
  DoubleReal getScore(VertexDescriptor vertex_idx)
  {
    return seqan::property(vertices_, vertex_idx).score;
  }

  /// return the interpretations of a node
  void getInterpretations(VertexDescriptor vertex_idx, InterpretVec &interprets)
  {
    interprets = seqan::property(vertices_, vertex_idx).peak_interprets;
  }

  /// return the real mass of edge from node i to j
  DoubleReal getRealEdgeMass(VertexDescriptor idx_l, VertexDescriptor idx_r) const
  {
    EdgeDescriptor ed = seqan::findEdge(graph, idx_l, idx_r);
    if (!ed)
      return -1;
    else
      return seqan::property(edges_, ed).real_mass;
  }

  /// return the real mass of edge from node i to j
  DoubleReal getIntEdgeMass(VertexDescriptor i, VertexDescriptor j) const
  {
    EdgeDescriptor ed = seqan::findEdge(graph, i, j);
    if (!ed)
      return -1;
    else
      return seqan::property(edges_, ed).int_mass;
  }

  /// return the number of clusters (in our case actually the number of undirected edges)
  UInt getNumberOfClusters() const
  {
    return conflict_edges_.size();
  }

  /// return the clusters, for each node the id of the clusters containing it
  const std::vector<Size>& getClusters(VertexDescriptor vertex_idx) const
  {    
    return conflict_edge_ids_[vertex_idx];
  }

  /// return a list of nodes contained in a certain cluster
  const std::vector<VertexDescriptor> & getConflictingNodes(VertexDescriptor vertex) const
  {
    return conflicting_nodes_[vertex];
  }

  /// return a list of nodes contained in a certain cluster
  void getConflictingNodesByCluster(std::vector<VertexDescriptor> &vertices, Size cluster_it) const
  {
    vertices.push_back(conflict_edges_[cluster_it].left_node);
    vertices.push_back(conflict_edges_[cluster_it].right_node);
  }

  const std::vector<VertexDescriptor>& getTopologicalOrdering()
  {
    return ordering_;
  }

  void activateVertex(VertexDescriptor v)
  {
    seqan::assignProperty(vertex_active_, v, true);
  }

  void deactivateVertex(VertexDescriptor v)
  {
    seqan::assignProperty(vertex_active_, v, false);
  }

  void activateEdge(EdgeDescriptor e)
  {
    seqan::assignProperty(edge_active_, e, true);
  }

  void deactivateEdge(EdgeDescriptor e)
  {
    seqan::assignProperty(edge_active_, e, false);
  }


  //TODO can also be removed as function and be done automatically in create node set
  /// score nodes using the learned parameters stored in file passed as parameter
  /// if add is true, then the scores accoring to scoring_func are added to previous node scores weighted by factor weight
  /// Nodes with a score below threshold are removed
  UInt scoreNodes(BayesScoring &scoring_func, const PeakSpectrum &spectrum, bool add = false,
      DoubleReal weight = 1.0, DoubleReal threshold = (-1) * std::numeric_limits<float>::max());


//protected members
protected:

  /// vector containing all undirected edges the contained nodes.
  /// their position is the conflict_edge_ids index stored fore every node
  std::vector<ConflictEdge> conflict_edges_;

  std::vector<std::vector<Size> > conflict_edge_ids_;

  /// adjacency matrix for the conflict edges
  std::vector<std::vector<bool> >conflict_edge_matrix_;


  std::vector<std::vector<VertexDescriptor> > conflicting_nodes_;

  /// vertices of the graph
  std::vector <Node> vertices_;

  /// edges of the graph
  std::vector<Edge> edges_;

  /// the ordering of vertices (is trivial but required for efficient call of seqan dag shortest path)
  std::vector<VertexDescriptor> ordering_;

  /// vector of bools for each vertex (true --> vertex is active (default), false --> vertex is inactive)
  seqan::String<bool> vertex_active_;

  /// vector of bools for each edge (true --> edge is active (default), false --> edge is inactive)
  seqan::String<bool> edge_active_;


};//SpectrumGraphSeqan


}//namespace



#endif // SPECTRUMGRAPH_SEQAN_H
