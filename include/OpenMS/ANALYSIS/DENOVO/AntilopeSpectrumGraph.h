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
  class MatrixSpectrumGraph;

  class SpectrumGraphSeqan: public DefaultParamHandler
  {
    friend class MatrixSpectrumGraph;

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
    
    //access for every cluster C all nodes in C
    typedef std::vector<std::vector<VertexDescriptor> > TConflictClusterList;

    //access for every node v all clusters containing v
    typedef std::vector<std::vector<Size> > TConflictClusterInvList;

    //access for every node all conflicting node ids
    typedef std::vector<std::vector<VertexDescriptor> > TConflictingNodesList;


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

    enum VERTEXSTATUS
    {
      FREE,
      ACTIVE,
      INACTIVE
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
#ifdef Debug
        std::cout << "created node: " << rmass << "  " << imass << "  Score-> " << inscore << " " << peak_id << " "
                  << type << std::endl;
#endif
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
    
    /// return score of edge
    DoubleReal getEdgeWeight(EdgeDescriptor) const;

    /// return all edge weights (i,j)
    void getEdgeWeights(seqan::String<DoubleReal> &weights) const;

    /// return whether two nodes i and j are connected via undirected edge --> are complementary
    bool hasUndirectedEdge(VertexDescriptor, VertexDescriptor) const; // TODO: get rid off

    /// return whether two nodes i and j are connected via directed edge
    EdgeDescriptor findDirectedEdge(VertexDescriptor, VertexDescriptor) const;

    /// adding a node
    //UInt addNodeToGraph(DoubleReal, int, DoubleReal, DoubleReal = -1);

    /// method to construct the spectrum graph (generate set of edges by editing adjacency lists)
    UInt createEdgeSet(const IdSetup::BoolVec &A, const IdSetup::UIntVec &A_len);

    /// construct the set of undirected edges
    UInt createUdirEdgeSet();

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
      return conflict_clusters.size();
    }

    /// return the clusters, for each node the id of the clusters containing it
    const std::vector<Size>& getClusters(VertexDescriptor vertex_idx) const
    {
      return conflict_clusters_inv[vertex_idx];
    }

    /// return a list of nodes contained in a certain cluster
    const std::vector<VertexDescriptor>& getConflictingNodes(VertexDescriptor vertex_id) const
    {
      return conflicting_nodes[vertex_id];
    }

    /// return a list of nodes contained in a certain cluster
    const std::vector<VertexDescriptor>& getCluster(Size cluster_id) const
    {
      return conflict_clusters[cluster_id];
    }

    const seqan::String<VertexDescriptor>& getTopologicalOrdering() const
    {
      return ordering_;
    }
    
    void getMassesForPaths(const std::vector<std::vector<VertexDescriptor> > &paths, std::vector<std::vector<UInt> >&prms) const
    {
      prms.resize(paths.size());
      std::vector<SpectrumGraphSeqan::VertexDescriptor>::const_iterator left, right;
      for (Size i = 0; i < paths.size(); ++i)
      {
        prms[i].clear();
        prms[i].reserve(paths[i].size());
        right = paths[i].begin();
        left = right++;
        for (; right != paths[i].end(); ++left, ++right)
          prms[i].push_back(getIntEdgeMass(*left, *right));
      }
    }
    
    void setEdgeWeights();

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
    
    bool isPathFeasible(const std::vector<VertexDescriptor>& path)
    {
      std::set<Size> clusters;
      
      for (Size i = 0; i < path.size(); ++i)
      {
        for (Size j = 0; j < conflict_clusters_inv[path[i]].size(); ++j)
        {
          if (!clusters.insert(conflict_clusters_inv[path[i]][j]).second)
            return false;
        }
      }
      return true;
    }
    
//    const std::vector<EdgeDescriptor>& getSpanningEdges(VertexDescriptor v) const
//    {
//      return spanning_edges_[v];
//    }
    
//    void createSpanningEdges();

    template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
    static void dagShortestPathST(seqan::Graph<TSpec> const& g,
                                  TVertexDescriptor const source,
                                  TVertexDescriptor const target,
                                  TWeightMap const& weight,
                                  TPredecessorMap& predecessor,
                                  TDistanceMap& distance,
                                  seqan::String<TVertexDescriptor> const& order);

//    template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
//    static void dagShortestPathST(seqan::Graph<TSpec> const& g,
//                                  TVertexDescriptor const source,
//                                  TVertexDescriptor const target,
//                                  TWeightMap const& weight,
//                                  TPredecessorMap& predecessor,
//                                  TDistanceMap& distance);
    template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
    static void dagShortestPathST(seqan::Graph<TSpec> const& g,
                                  TVertexDescriptor const source,
                                  TVertexDescriptor const target,
                                  TWeightMap const& weight,
                                  TPredecessorMap& predecessor,
                                  TDistanceMap& distance,
                                  seqan::String<bool> const & is_vertex_inactive,
                                  seqan::String<bool> const & is_edge_inactive,
                                  std::set<TVertexDescriptor> const &forced_vertices,
                                  seqan::String<TVertexDescriptor> const& order
                                  );


    //TODO can also be removed as function and be done automatically in create node set
    /// score nodes using the learned parameters stored in file passed as parameter
    /// if add is true, then the scores accoring to scoring_func are added to previous node scores weighted by factor weight
    /// Nodes with a score below threshold are removed
    UInt scoreNodes(BayesScoring &scoring_func, const PeakSpectrum &spectrum, bool add = false,
                    DoubleReal weight = 1.0, DoubleReal threshold = (-1) * std::numeric_limits<float>::max());


    //protected members
  protected:
    
    std::vector<std::vector<bool> > pairwise_conflicting;
    TConflictClusterInvList conflict_clusters_inv;
    TConflictClusterList conflict_clusters;
    TConflictingNodesList conflicting_nodes;


    //std::vector<std::vector<VertexDescriptor> > conflicting_nodes_;
    
    //std::vector<std::vector<EdgeDescriptor> > spanning_edges_;

    /// vertices of the graph
    std::vector <Node> vertices_;

    /// edges of the graph
    std::vector<Edge> edges_;

    /// the ordering of vertices (is trivial but required for efficient call of seqan dag shortest path)
    //std::vector<VertexDescriptor> ordering_;

    /// vector of bools for each vertex (true --> vertex is active (default), false --> vertex is inactive)
    seqan::String<bool> vertex_active_;

    /// vector of bools for each edge (true --> edge is active (default), false --> edge is inactive)
    seqan::String<bool> edge_active_;
    
    /// edge weights of directed edges
    seqan::String<DoubleReal> edge_weights_;
    
    seqan::String<VertexDescriptor> ordering_;


  };//SpectrumGraphSeqan
  
  
  class MatrixSpectrumGraph : public SpectrumGraphSeqan
  {
    protected:
    
    
    bool hasUndirectedEdge(VertexDescriptor, VertexDescriptor) const
    {
      //this is not available for matrix spectrum graph -- todo refactor
      return false;
    };
    
    
    // last node in lower half of spectrum graph
    Size border_node;
    
    std::vector<std::vector<VertexDescriptor> > matrix_vertices;
    
    std::vector<std::pair<Size,Size> > matrix_vertices_inv;
    
  public:
    
    UInt createMatrixGraph(const SpectrumGraphSeqan &base_graph)
    {
      //identify border node
      Size base_last_id = base_graph.vertices_.size() - 1;
      DoubleReal pm = base_graph.vertices_.back().real_mass;
      DoubleReal pm_2 = pm / 2;
      matrix_vertices.assign(base_graph.size(), std::vector<VertexDescriptor>(base_graph.size(), -1u));
      
      //border node is first node >= than 0.5 parent mass
      for (border_node = 0; base_graph.vertices_[border_node].real_mass < pm_2; ++border_node);
      
      TGraph graph_tmp;
      
      std::vector<std::pair<Size,Size> > matrix_vertices_inv_tmp;
      matrix_vertices_inv_tmp.reserve(base_graph.size() * base_graph.size());
      
      //create vertex set
      for (Size i = 0; i < border_node; ++i)
      {
        for (Size j = base_last_id; j >= border_node; --j)
        {
          if ( !base_graph.pairwise_conflicting[i][j] &&
              fabs(base_graph.vertices_[i].real_mass + base_graph.vertices_[j].real_mass - pm) <= 400)
          {
            matrix_vertices[i][j] = seqan::addVertex(graph_tmp);
            matrix_vertices_inv_tmp.push_back(std::make_pair(i,j));
          }
        }
      }
      
      VertexDescriptor matrix_source = matrix_vertices[0][base_last_id];
      VertexDescriptor matrix_sink = seqan::addVertex(graph_tmp);
     

      DoubleReal max_dist_e = 0;

      std::vector<std::pair<EdgeDescriptor, DoubleReal> > edge_weights_tmp;
      edge_weights_tmp.reserve(seqan::numEdges(base_graph.graph) * base_graph.size());
      
      EdgeIterator base_it(base_graph.graph);
      while (!seqan::atEnd(base_it))
      {
        VertexDescriptor s = seqan::sourceVertex(base_it);
        VertexDescriptor t = seqan::targetVertex(base_it);
        max_dist_e = std::max(max_dist_e, fabs(base_graph.vertices_[s].real_mass - base_graph.vertices_[t].real_mass));

        //std::cerr << "s: " << s << " t: " << t << std::endl;
        
        if (t < border_node) // within lower half of spectrum graph
        {
          for (Size i = border_node; i < base_graph.size(); ++i)
          {
            if (matrix_vertices[s][i] == -1 || matrix_vertices[t][i] == -1)
              continue;
            
            //check whether target vertex is below diagonal
            if (base_graph.vertices_[t].real_mass + base_graph.vertices_[i].real_mass > pm)
            {
              //std::cerr << "(" << s << ',' << i << ") -> (" << t << ',' << i << ")" << std::endl;
              EdgeDescriptor edge = seqan::addEdge(graph_tmp, matrix_vertices[s][i], matrix_vertices[t][i]);
              edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(base_graph.edge_weights_, *base_it)));
            }
          }
        }
        else if (s >= border_node) // within upper half of spectrum graph
        {
          for (Size i = 0; i < border_node; ++i)
          {
            if (matrix_vertices[i][s] == -1 || matrix_vertices[i][t] == -1)
              continue;
            
            //check whether target vertex is above diagonal
            if (base_graph.vertices_[s].real_mass + base_graph.vertices_[i].real_mass <= pm)
            {
              //std::cerr << "(" << i << ',' << t << ") -> (" << i << ',' << s << ")" << std::endl;
              EdgeDescriptor edge = seqan::addEdge(graph_tmp, matrix_vertices[i][t], matrix_vertices[i][s]);
              edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(base_graph.edge_weights_, *base_it)));
            }
          }
        }
        else //possible accepting vertex in matrix spectrum graph
        {
          if (matrix_vertices[s][t] != -1)
          {
            EdgeDescriptor edge = seqan::addEdge(graph_tmp, matrix_vertices[s][t], matrix_sink);
            edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(base_graph.edge_weights_, *base_it)));
          }
        }
        seqan::goNext(base_it);
      }
      
      // edge weights
      seqan::String<DoubleReal> weights_tmp;
      seqan::resizeEdgeMap(graph_tmp, weights_tmp);
      for (Size i = 0; i < edge_weights_tmp.size(); ++i)
      {
        seqan::property(weights_tmp, edge_weights_tmp[i].first) = edge_weights_tmp[i].second;
      }

      // reduce graph by removing vertices not included in any source - sink path
      edge_weights_tmp.clear();
      
      std::cerr << "SpecGraph Size: " << seqan::numVertices(base_graph.graph) << " " << seqan::numEdges(base_graph.graph) << std::endl;
      std::cerr << "MatrixGraph Size: " << seqan::numVertices(graph_tmp) << " " << seqan::numEdges(graph_tmp) << std::endl;

      seqan::String<VertexDescriptor> pred;
      seqan::String<Size> distF, distR;
      seqan::breadthFirstSearch(graph_tmp, matrix_source, pred, distF);
      TGraph graph_rev;
      seqan::transpose(graph_tmp, graph_rev);
      
      seqan::breadthFirstSearch(graph_rev, matrix_sink, pred, distR);
      
      DoubleReal tt = clock();
      std::cerr << "1" << std::endl;
      
      std::vector<VertexDescriptor> new_v_ids(seqan::numVertices(graph_tmp), -1u);
      VertexIterator v_it(graph_tmp);
      
      while (!seqan::atEnd(v_it))
      {
        if (seqan::property(distF, *v_it) != _getInfinityDistance(distF) &&
            seqan::property(distR, *v_it) != _getInfinityDistance(distR))
        {
          //include into final graph
          if (*v_it != 0)
          {
            seqan::addVertex(graph); //this is used as entrance vertex
            matrix_vertices_inv.push_back(std::make_pair(-1u, -1u));
          }
          new_v_ids[*v_it] = seqan::addVertex(graph);
          if (*v_it != 0)
          {
            EdgeDescriptor edge = seqan::addEdge(graph, new_v_ids[*v_it]-1, new_v_ids[*v_it]);
            edge_weights_tmp.push_back(std::make_pair(edge, 0));
          }

          if (*v_it != matrix_sink)
          {
            matrix_vertices[matrix_vertices_inv_tmp[*v_it].first][matrix_vertices_inv_tmp[*v_it].second] = new_v_ids[*v_it];
            matrix_vertices_inv.push_back(matrix_vertices_inv_tmp[*v_it]);
          }
        }
        else
        {
          //exclude from final graph
          matrix_vertices[matrix_vertices_inv_tmp[*v_it].first][matrix_vertices_inv_tmp[*v_it].second] = -1u;
        }
        seqan::goNext(v_it);
      }
      
      matrix_source = new_v_ids[matrix_source];
      matrix_sink = new_v_ids[matrix_sink];
      
      std::vector<std::set<VertexDescriptor> > incoming_v(base_graph.size());
      
      EdgeIterator e_it(graph_tmp);
      while (!seqan::atEnd(e_it))
      {
        VertexDescriptor s = sourceVertex(e_it);
        VertexDescriptor t = targetVertex(e_it);
        
        if (new_v_ids[s] != -1u && new_v_ids[t] != -1u)
        {
          if (new_v_ids[t] == matrix_sink)
          {
            EdgeDescriptor edge = seqan::addEdge(graph, new_v_ids[s], new_v_ids[t]);
            edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(weights_tmp, *e_it)));
          }
          else
          {
            VertexDescriptor incoming = new_v_ids[t]-1;
            if (matrix_vertices_inv[new_v_ids[s]].first != matrix_vertices_inv[new_v_ids[t]].first)
            {
              //vertical edge
              incoming_v[matrix_vertices_inv[new_v_ids[t]].first].insert(incoming);
            }
            else
            {
              //horizontal edge
              incoming_v[matrix_vertices_inv[new_v_ids[t]].second].insert(incoming);
            }
            
            EdgeDescriptor edge1 = seqan::addEdge(graph, new_v_ids[s], incoming);
            edge_weights_tmp.push_back(std::make_pair(edge1, seqan::property(weights_tmp, *e_it)));
          }
        }
        seqan::goNext(e_it);
      }
      
      // edge weights
      seqan::resizeEdgeMap(graph, edge_weights_);
      for (Size i = 0; i < edge_weights_tmp.size(); ++i)
      {
        seqan::property(edge_weights_, edge_weights_tmp[i].first) = edge_weights_tmp[i].second;
      }
      
      std::cerr << "Took: " << (clock() - tt)/CLOCKS_PER_SEC << std::endl;
      std::cerr << "2" << std::endl;
      
      
      // conflict clusters
      conflict_clusters.clear();
      conflict_clusters_inv.clear();
      
      conflict_clusters.resize(base_graph.conflict_clusters.size());
      conflict_clusters_inv.resize(seqan::numVertices(graph));
      conflicting_nodes.resize(seqan::numVertices(graph));
      for (Size i = 0; i < base_graph.conflict_clusters.size(); ++i)
      {
        for (Size j = 0; j < base_graph.conflict_clusters[i].size(); ++j)
        {
          VertexDescriptor node_id = base_graph.conflict_clusters[i][j];
//          std::cout << "base info cluster: " << node_id << " " << base_graph.conflict_clusters_inv[node_id].size() << std::endl;
          conflict_clusters[i].insert(conflict_clusters[i].end(), incoming_v[node_id].begin(), incoming_v[node_id].end());
        }
        
        for (Size j = 0; j < conflict_clusters[i].size(); ++j)
        {
           conflict_clusters_inv[conflict_clusters[i][j]].push_back(i);
        }
      }
      

      for (Size i = 0; i < conflict_clusters.size(); ++i)
      {
        for (Size j = 0; j < conflict_clusters[i].size(); ++j)
        {
          for (Size k = j+1; k < conflict_clusters[i].size(); ++k)
          {
            conflicting_nodes[conflict_clusters[i][j]].push_back(conflict_clusters[i][k]);
            conflicting_nodes[conflict_clusters[i][k]].push_back(conflict_clusters[i][j]);
          }
        }
      }

      //remove duplicate entries in conflicting_nodes
      for (Size i = 0; i < vertices_.size(); ++i)
      {
        std::vector<VertexDescriptor> &cni = conflicting_nodes[i];
        std::sort(cni.begin(), cni.end());
        cni.erase(std::unique(cni.begin(), cni.end()),cni.end());
      }
      
      seqan::resize(ordering_, seqan::numVertices(graph));
      for (Size i = 0; i < seqan::numVertices(graph); ++i)
        ordering_[i] = i;
      
        

      
//      
//      EdgeIterator e_it(graph);
//      std::vector<VertexDescriptor> edges_tmp;
//      while (!seqan::atEnd(e_it))
//      {
//        VertexDescriptor s = sourceVertex(e_it);
//        VertexDescriptor t = targetVertex(e_it);
//        
//        if (seqan::property(distF, s) != _getInfinityDistance(distF) &&
//            seqan::property(distR, s) != _getInfinityDistance(distR) &&
//            seqan::property(distF, t) != _getInfinityDistance(distF) &&
//            seqan::property(distR, t) != _getInfinityDistance(distR))
//        {
//          edges_tmp.push_back(s);
//          edges_tmp.push_back(t);
//        }
//      }
//      
//      seqan::clear(graph);
//      seqan::addEdges(graph, edges_tmp, edges_tmp.size()/2);
//      seqan::resizeEdgeMap(graph, edge_weights_);
//
//      EdgeIterator e_it2(graph);
//      while (!seqan::atEnd(e_it2))
//      {
//        VertexDescriptor s = sourceVertex(e_it);
//        VertexDescriptor t = targetVertex(e_it);
//        
//        if (matrix_vertices_inv[s].first != matrix_vertices_inv[t].first)
//          seqan::property(edge_weights_, *e_it2) = edge_weights_tmp[i].second;
//        
//      }
//      for (Size i = 0; i < edge_weights_tmp.size(); ++i)
//      {
//        seqan::property(edge_weights_, edge_weights_tmp[i].first) = edge_weights_tmp[i].second;
//      }
//
//      
//      
//      
//      VertexIterator v_it(graph);
//      ordering_.clear();
//      ordering_.reserve(seqan::numVertices(graph));
//      
//      while (!seqan::atEnd(v_it))
//      {
//        VertexDescriptor v = *v_it;
//        //std::cerr << v << std::endl;
//        seqan::goNext(v_it);
//        
//        if (seqan::property(distF, v) == _getInfinityDistance(distF) || seqan::property(distR, v) == _getInfinityDistance(distR))
//        {
//          seqan::removeVertex(graph, v);
//          matrix_vertices[matrix_vertices_inv[v].first][matrix_vertices_inv[v].second] = -1u;
//        }
//      }
//      std::cerr << "Took: " << (clock() - tt)/CLOCKS_PER_SEC << std::endl;
//      std::cerr << "2" << std::endl;

//      DEBUG STUFF
//      seqan::goBegin(v_it);
//      DoubleReal max_dist = 0;
//      while (!seqan::atEnd(v_it))
//      {
//        if (*v_it == matrix_sink)
//        {
//          seqan::goNext(v_it);
//          continue;
//        }
//        
//        DoubleReal tmpd = fabs(base_graph.vertices_[matrix_vertices_inv[*v_it].first].real_mass + base_graph.vertices_[matrix_vertices_inv[*v_it].second].real_mass - pm);
//        
//        std::cerr << matrix_vertices_inv[*v_it].first << " : " << matrix_vertices_inv[*v_it].second << " " << tmpd << std::endl;
//        max_dist = std::max(max_dist, tmpd);
//        seqan::goNext(v_it);
//      }
//      std::cerr << "MAX DIST: " << max_dist << "  Edge: " << max_dist_e << std::endl;

      
//      
//      // conflict clusters
//      for (Size i = 0; i < base_graph.conflict_clusters.size(); ++i)
//      {
//        for (Size j = 0; j < base_graph.conflict_clusters[i].size(); ++j)
//        {
//          VertexDescriptor node_id = base_graph.conflict_clusters[i][j];
//          
//          if (base_graph.vertices_[node_id].real_mass < pm_2)
//          {
//            for (Size k = border_node; k < base_graph.size(); ++k)
//            {
//              VertexDescriptor matrix_node_id = matrix_vertices[node_id][k];
//              if (matrix_node_id != -1u)
//              {
//                conflict_clusters[i].push_back(matrix_node_id);
//                conflict_clusters_inv[matrix_node_id].push_back(i);
//              }
//            }
//          }
//          else
//          {
//            for (Size k = 0; k < border_node; ++k)
//            {
//              VertexDescriptor matrix_node_id = matrix_vertices[k][node_id];
//              if (matrix_node_id != -1u)
//              {
//                conflict_clusters[i].push_back(matrix_node_id);
//                conflict_clusters_inv[matrix_node_id].push_back(i);
//              }
//            }
//          }
//        }
//      }
//      for (Size i = 0; i < conflict_clusters.size(); ++i)
//      {
//        for (Size j = 0; j < conflict_clusters[i].size(); ++j)
//        {
//          conflict_clusters_inv[conflict_clusters[i][j]].push_back(i);
//          for (Size k = j+1; k < conflict_clusters[i].size(); ++k)
//          {
//            conflicting_nodes[conflict_clusters[i][j]].push_back(conflict_clusters[i][k]);
//            conflicting_nodes[conflict_clusters[i][k]].push_back(conflict_clusters[i][j]);
//            
//            //pairwise_conflicting[conflict_clusters[i][j]][conflict_clusters[i][k]] = true;
//            //pairwise_conflicting[conflict_clusters[i][k]][conflict_clusters[i][j]] = true;
//          }
//        }
//      }
//
//      //remove duplication entries in conflicting_nodes
//      for (Size i = 0; i < vertices_.size(); ++i)
//      {
//        std::vector<VertexDescriptor> &cni = conflicting_nodes[i];
//        std::sort(cni.begin(), cni.end());
//        cni.resize(std::distance(cni.begin(), std::unique(cni.begin(), cni.end() )));
//      }


//      //conflict_clusters.resize(base_graph.conflict_clusters.size());
//      std::set<std::pair<VertexDescriptor, VertexDescriptor> > conflict_check;
//      
//      conflict_clusters_inv.resize(seqan::numVertices(graph));
//      conflicting_nodes.resize(seqan::numVertices(graph));
//      
//      std::cerr << "Clusters: " << base_graph.conflict_clusters.size() << std::endl;
//      for (Size i = 0; i < base_graph.conflict_clusters.size(); ++i)
//      {
//        std::cerr << i << " " << base_graph.conflict_clusters.size() << std::endl;
//        for (Size j = 0; j < base_graph.conflict_clusters[i].size(); ++j)
//        {
//          VertexDescriptor node_id_j = base_graph.conflict_clusters[i][j];
//          DoubleReal mass_j = base_graph.vertices_[node_id_j].real_mass;
//          
//          for (Size k = j+1; k < base_graph.conflict_clusters[i].size(); ++k)
//          {
//            VertexDescriptor node_id_k = base_graph.conflict_clusters[i][k];
//            DoubleReal mass_k = base_graph.vertices_[node_id_k].real_mass;
//            
//            if (node_id_k < node_id_j)
//            {
//              //swap
//              std::swap(node_id_k, node_id_j);
//              std::swap(mass_k, mass_j);
//            }
//            
//            if (node_id_k < border_node && (mass_k - mass_j) < 57)
//              continue; //both in lower half and too close to appear in same path
//            
//            if (node_id_j >= border_node && (mass_k - mass_j) < 57)
//              continue; //both in upper half and too close to appear in same path
//
//            if (node_id_j < border_node && node_id_k >= border_node && fabs(mass_j + mass_k - pm) < 57)
//              continue; // in different halves but within the basic ion types mass corridor
//            
//            std::cerr << node_id_j << " " << mass_j << " : " << node_id_k << " " << mass_k << " cs: " << base_graph.conflict_clusters[i].size() << " PM: " << pm <<  std::endl;
//            std::cerr << mass_k - mass_j << "  " << mass_j + mass_k - pm << std::endl;
//            
//            if (node_id_k < border_node)
//            {
//              for (Size p1 = border_node; p1 < base_graph.size(); ++p1)
//              {
//                VertexDescriptor matrix_node_id_j = matrix_vertices[node_id_j][p1];
//                
//                if (matrix_node_id_j == -1u)
//                  continue;
//                
//                for (Size p2 = 0; p2 <= p1; ++p2)
//                {
//                  VertexDescriptor matrix_node_id_k = matrix_vertices[node_id_k][p2];
//                  if (matrix_node_id_k == -1u)
//                    continue;
//
//                  if (conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
//                  {
//                    conflict_clusters.resize(conflict_clusters.size()+1);
//                    conflict_clusters.back().push_back(matrix_node_id_j);
//                    conflict_clusters.back().push_back(matrix_node_id_k);
//                    
////                    if (!conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
////                    {
////                      std::cout << "redundant conflict1: " << matrix_node_id_j << ", " << matrix_node_id_k << std::endl;
////                      exit(1);
////                    }
//                    
//                    conflict_clusters_inv[matrix_node_id_j].push_back(conflict_clusters.size()-1);
//                    conflict_clusters_inv[matrix_node_id_k].push_back(conflict_clusters.size()-1);
//                  }
//                }
//              }
//            }
//            else if (node_id_j >= border_node)
//            {
//              for (Size p1 = 0; p1 < border_node; ++p1)
//              {
//                VertexDescriptor matrix_node_id_j = matrix_vertices[p1][node_id_j];
//                
//                if (matrix_node_id_j == -1u)
//                  continue;
//                
//                for (Size p2 = 0; p2 <= p1; ++p2)
//                {
//                  VertexDescriptor matrix_node_id_k = matrix_vertices[p2][node_id_k];
//                  if (matrix_node_id_k == -1u)
//                    continue;
//                  
//                  if (conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
//                  {
//                    conflict_clusters.resize(conflict_clusters.size()+1);
//                    conflict_clusters.back().push_back(matrix_node_id_j);
//                    conflict_clusters.back().push_back(matrix_node_id_k);
//                    
//  //                  if (conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
//  //                  {
//  //                    std::cout << "redundant conflict2: " << matrix_node_id_j << ", " << matrix_node_id_k << std::endl;
//  //                    exit(1);
//  //                  }
//                    
//                    conflict_clusters_inv[matrix_node_id_j].push_back(conflict_clusters.size()-1);
//                    conflict_clusters_inv[matrix_node_id_k].push_back(conflict_clusters.size()-1);
//                  }
//                }
//              }
//            }
//            else
//            {
//              for (Size p1 = border_node; p1 < node_id_k; ++p1)
//              {
//                VertexDescriptor matrix_node_id_j = matrix_vertices[node_id_j][p1];
//                
//                if (matrix_node_id_j == -1u)
//                  continue;
//                
//                for (Size p2 = 0; p2 < node_id_j; ++p2)
//                {
//                  VertexDescriptor matrix_node_id_k = matrix_vertices[p2][node_id_k];
//                  if (matrix_node_id_k == -1u)
//                    continue;
//                  
//                  if (conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
//                  {
//                    conflict_clusters.resize(conflict_clusters.size()+1);
//                    conflict_clusters.back().push_back(matrix_node_id_j);
//                    conflict_clusters.back().push_back(matrix_node_id_k);
//                    
////                    if (!conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
////                    {
////                      std::cout << "redundant conflict3: " << matrix_node_id_j << ", " << matrix_node_id_k << std::endl;
////                      exit(1);
////                    }
//                    
//                    conflict_clusters_inv[matrix_node_id_j].push_back(conflict_clusters.size()-1);
//                    conflict_clusters_inv[matrix_node_id_k].push_back(conflict_clusters.size()-1);
//                  }
//                }
//              }
//              for (Size p1 = node_id_k+1; p1 < base_graph.size(); ++p1)
//              {
//                VertexDescriptor matrix_node_id_j = matrix_vertices[node_id_j][p1];
//                
//                if (matrix_node_id_j == -1u)
//                  continue;
//                
//                for (Size p2 = node_id_j + 1; p2 < border_node; ++p2)
//                {
//                  VertexDescriptor matrix_node_id_k = matrix_vertices[p2][node_id_k];
//                  if (matrix_node_id_k == -1u)
//                    continue;
//
//                  if (conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
//                  {
//                    conflict_clusters.resize(conflict_clusters.size()+1);
//                    conflict_clusters.back().push_back(matrix_node_id_j);
//                    conflict_clusters.back().push_back(matrix_node_id_k);
//                    
////                    if (!conflict_check.insert(std::make_pair(matrix_node_id_j, matrix_node_id_k)).second)
////                    {
////                      std::cout << "redundant conflict4: " << matrix_node_id_j << ", " << matrix_node_id_k << std::endl;
////                      exit(1);
////                    }
//                    
//                    conflict_clusters_inv[matrix_node_id_j].push_back(conflict_clusters.size()-1);
//                    conflict_clusters_inv[matrix_node_id_k].push_back(conflict_clusters.size()-1);
//                  }
//                }
//              }
//            }
//          }
//        }
//      }
      
      
      std::cerr << "SpecGraph Size: " << seqan::numVertices(base_graph.graph) << " " << seqan::numEdges(base_graph.graph) << std::endl;
      std::cerr << "MatrixGraph Size: " << seqan::numVertices(graph) << " " << seqan::numEdges(graph) << " " << conflict_clusters.size() <<  std::endl;
      //exit(1);
      return 0;
    }
    
    UInt createMatrixGraph2(const SpectrumGraphSeqan &base_graph)
    {
      //identify border node
      Size base_last_id = base_graph.vertices_.size() - 1;
      DoubleReal pm = base_graph.vertices_.back().real_mass;
      DoubleReal pm_2 = pm / 2;
      matrix_vertices.assign(base_graph.size(), std::vector<VertexDescriptor>(base_graph.size(), -1u));
      
      //border node is first node >= than 0.5 parent mass
      for (border_node = 0; base_graph.vertices_[border_node].real_mass < pm_2; ++border_node);
      
      TGraph& graph_tmp = graph;
      
//      std::vector<std::pair<Size,Size> > matrix_vertices_inv_tmp;
      matrix_vertices_inv.clear();
      matrix_vertices_inv.reserve(base_graph.size() * base_graph.size());
      
      //create vertex set
      for (Size i = 0; i < border_node; ++i)
      {
        for (Size j = base_last_id; j >= border_node; --j)
        {
          if ( !base_graph.pairwise_conflicting[i][j] &&
              fabs(base_graph.vertices_[i].real_mass + base_graph.vertices_[j].real_mass - pm) <= 400)
          {
            matrix_vertices[i][j] = seqan::addVertex(graph_tmp);
            matrix_vertices_inv.push_back(std::make_pair(i,j));
          }
        }
      }
      
      VertexDescriptor matrix_source = matrix_vertices[0][base_last_id];
      VertexDescriptor matrix_sink = seqan::addVertex(graph_tmp);
      
      
      DoubleReal max_dist_e = 0;
      
      std::vector<std::pair<EdgeDescriptor, DoubleReal> > edge_weights_tmp;
      edge_weights_tmp.reserve(seqan::numEdges(base_graph.graph) * base_graph.size());
      
      EdgeIterator base_it(base_graph.graph);
      while (!seqan::atEnd(base_it))
      {
        VertexDescriptor s = seqan::sourceVertex(base_it);
        VertexDescriptor t = seqan::targetVertex(base_it);
        max_dist_e = std::max(max_dist_e, fabs(base_graph.vertices_[s].real_mass - base_graph.vertices_[t].real_mass));
        
        //std::cerr << "s: " << s << " t: " << t << std::endl;
        
        if (t < border_node) // within lower half of spectrum graph
        {
          for (Size i = border_node; i < base_graph.size(); ++i)
          {
            if (matrix_vertices[s][i] == -1 || matrix_vertices[t][i] == -1)
              continue;
            
            //check whether target vertex is below diagonal
            if (base_graph.vertices_[t].real_mass + base_graph.vertices_[i].real_mass > pm)
            {
              //std::cerr << "(" << s << ',' << i << ") -> (" << t << ',' << i << ")" << std::endl;
              EdgeDescriptor edge = seqan::addEdge(graph_tmp, matrix_vertices[s][i], matrix_vertices[t][i]);
              edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(base_graph.edge_weights_, *base_it)));
            }
          }
        }
        else if (s >= border_node) // within upper half of spectrum graph
        {
          for (Size i = 0; i < border_node; ++i)
          {
            if (matrix_vertices[i][s] == -1 || matrix_vertices[i][t] == -1)
              continue;
            
            //check whether target vertex is above diagonal
            if (base_graph.vertices_[s].real_mass + base_graph.vertices_[i].real_mass <= pm)
            {
              //std::cerr << "(" << i << ',' << t << ") -> (" << i << ',' << s << ")" << std::endl;
              EdgeDescriptor edge = seqan::addEdge(graph_tmp, matrix_vertices[i][t], matrix_vertices[i][s]);
              edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(base_graph.edge_weights_, *base_it)));
            }
          }
        }
        else //possible accepting vertex in matrix spectrum graph
        {
          if (matrix_vertices[s][t] != -1)
          {
            EdgeDescriptor edge = seqan::addEdge(graph_tmp, matrix_vertices[s][t], matrix_sink);
            edge_weights_tmp.push_back(std::make_pair(edge, seqan::property(base_graph.edge_weights_, *base_it)));
          }
        }
        seqan::goNext(base_it);
      }
      
      // edge weights
      seqan::resizeEdgeMap(graph_tmp, edge_weights_);
      for (Size i = 0; i < edge_weights_tmp.size(); ++i)
      {
        seqan::property(edge_weights_, edge_weights_tmp[i].first) = edge_weights_tmp[i].second;
      }
      
      std::cerr << "SpecGraph Size: " << seqan::numVertices(base_graph.graph) << " " << seqan::numEdges(base_graph.graph) << std::endl;
      std::cerr << "MatrixGraph Size: " << seqan::numVertices(graph_tmp) << " " << seqan::numEdges(graph_tmp) << std::endl;
      
      seqan::resize(ordering_, seqan::numVertices(graph));
      for (Size i = 0; i < seqan::numVertices(graph); ++i)
        ordering_[i] = i;
      
      
      conflict_clusters.clear();
      conflict_clusters_inv.clear();
      
      conflict_clusters.resize(base_graph.conflict_clusters.size());
      conflict_clusters_inv.resize(seqan::numVertices(graph));
      
      conflicting_nodes.resize(seqan::numVertices(graph));
      
      
      return 0;
    }
    
    
    void getMassesForPaths(const SpectrumGraphSeqan &base_graph,
                         const std::vector<std::vector<VertexDescriptor> > &paths,
                         std::vector<std::vector<UInt> >&prms) const
    {
      prms.resize(paths.size());

      for (Size i = 0; i < paths.size(); ++i)
      {
        //get path without fake vertices
        std::vector<VertexDescriptor> path;
        path.reserve(paths[i].size());
        
        for (Size j = 0; j < paths[i].size(); ++j)
        {
          if (matrix_vertices_inv[paths[i][j]].first != -1u)
            path.push_back(paths[i][j]);
        }
        
        prms[i].clear();
        prms[i].resize(path.size()-1);
        Size lpos = 0;
        Size rpos = path.size()-2;

        Size j = 0;
        for (; j+2 < path.size(); ++j)
        {
          const std::pair<VertexDescriptor, VertexDescriptor> &ij = matrix_vertices_inv[path[j]];
          const std::pair<VertexDescriptor, VertexDescriptor> &kl = matrix_vertices_inv[path[j+1]];
          
          if (ij.first != kl.first) //vertical edge
          {
            prms[i][lpos] = base_graph.getIntEdgeMass(ij.first, kl.first);
            ++lpos;
          }
          else //horizontal edge
          {
            prms[i][rpos] = base_graph.getIntEdgeMass(kl.second, ij.second);
            --rpos;
          }
        }
        //the diagonal spanning edge
        const std::pair<VertexDescriptor, VertexDescriptor> &ij = matrix_vertices_inv[path[j]];
        prms[i][lpos] = base_graph.getIntEdgeMass(ij.first, ij.second);
        
        for (Size j = 0; j < prms[i].size(); ++j)
          std::cout << prms[i][j] << std::endl;
      }
    }

    
    
    
  };
  

  template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
  void
  SpectrumGraphSeqan::dagShortestPathST(seqan::Graph<TSpec> const& g,
                                        TVertexDescriptor const source,
                                        TVertexDescriptor const target,
                                        TWeightMap const& weight,
                                        TPredecessorMap& predecessor,
                                        TDistanceMap& distance,
                                        seqan::String<TVertexDescriptor> const& order)
  {
    using namespace seqan;

    typedef typename seqan::EdgeDescriptor<seqan::Graph<TSpec> >::Type TEdgeDescriptor;
    typedef typename seqan::Iterator<seqan::Graph<TSpec>, seqan::EdgeIterator>::Type TEdgeIterator;
    typedef typename seqan::Iterator<seqan::Graph<TSpec>, seqan::OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename seqan::Iterator<const seqan::String<TVertexDescriptor>, Rooted>::Type TStringIterator;

    // Initialization
    resizeVertexMap(g,predecessor);
    resizeVertexMap(g,distance);

    _initializeSingleSource(g, source, weight, predecessor, distance);

    //DAG Shortest Paths
    TStringIterator it = begin(order);
    //jump over all vertices before source
    while(getValue(it) != source && !atEnd(it))
    {
      goNext(it);
    }
    //stop at target node
    while(getValue(it) != target && !atEnd(it)) {
      TOutEdgeIterator itout(g, getValue(it));
      for(;!atEnd(itout);++itout) {
        _relax(g,weight,predecessor, distance, getValue(it), getValue(itout));
      }
      goNext(it);
    }
  }
//
//  template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
//  void
//  SpectrumGraphSeqan::dagShortestPathST(seqan::Graph<TSpec> const& g,
//                                        TVertexDescriptor const source,
//                                        TVertexDescriptor const target,
//                                        TWeightMap const& weight,
//                                        TPredecessorMap& predecessor,
//                                        TDistanceMap& distance)
//  {
//    // Topological sort
//    seqan::String<TVertexDescriptor> order;
//    seqan::topologicalSort(g, order);
//
//    //call dagShortestPath with ordering
//    dagShortestPathST(g, source, target, weight, predecessor, distance, order);
//  }

  
  //Compute shortest source - target path in dag
  //vertex descriptors are implicitly assumed to have topological order - true for spectrum graphs by construction
  template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
  void SpectrumGraphSeqan:: dagShortestPathST(seqan::Graph<TSpec> const& g,
                                TVertexDescriptor const source,
                                TVertexDescriptor const target,
                                TWeightMap const& weight,
                                TPredecessorMap& predecessor,
                                TDistanceMap& distance,
                                seqan::String<bool> const & is_vertex_inactive,
                                seqan::String<bool> const & is_edge_inactive,
                                std::set<TVertexDescriptor> const &forced_vertices,
                                seqan::String<TVertexDescriptor> const& order
                                )
  {
    using namespace seqan;

    typedef typename seqan::EdgeDescriptor<seqan::Graph<TSpec> >::Type TEdgeDescriptor;
    typedef typename seqan::Iterator<seqan::Graph<TSpec>, seqan::VertexIterator>::Type TVertexIterator;
    typedef typename seqan::Iterator<seqan::Graph<TSpec>, seqan::EdgeIterator>::Type TEdgeIterator;
    typedef typename seqan::Iterator<seqan::Graph<TSpec>, seqan::OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename seqan::Iterator<const seqan::String<TVertexDescriptor>, Rooted>::Type TStringIterator;
   

    typedef typename std::set<TVertexDescriptor>::const_iterator TSetIter;

    // Initialization
    resizeVertexMap(g,predecessor);
    resizeVertexMap(g,distance);

    _initializeSingleSource(g, source, weight, predecessor, distance);

    TStringIterator it = begin(order);
    //TVertexIterator it(g);
    
    //jump over all vertices before source
    while(getValue(it) != source && !atEnd(it))
    {
      goNext(it);
    }
    
    
    TSetIter forced_it = forced_vertices.begin();
    const TSetIter forced_end = forced_vertices.end();
    TVertexDescriptor next_forced = -1u;

    if (forced_it != forced_end)  //if nonempty forced vertices
    {
      next_forced = *forced_it;
    }

    while(getValue(it) != target && !atEnd(it))
    {
      if (getProperty(is_vertex_inactive, getValue(it)))
      {
        goNext(it);
        continue;
      }
      
      if (getValue(it) == next_forced)
      {
        ++forced_it;
        
        if (forced_it != forced_end)
          next_forced = *forced_it;
        else
          next_forced = -1u;
      }
    

      TOutEdgeIterator itout(g, getValue(it));
      for( ; !atEnd(itout); ++itout)
      {
        if (targetVertex(itout) > next_forced)  //edge must not cross next forced vertex
          continue;
          
        if(!getProperty(is_edge_inactive, getValue(itout)))
        {
          _relax(g, weight, predecessor, distance, getValue(it), getValue(itout));
        }
      }
      goNext(it);
    }
  }


}//namespace



#endif // SPECTRUMGRAPH_SEQAN_H
