/***************************************************************************
 *   Copyright (C) 2005 by Gunnar W. Klau, Markus Bauer                    *
 *   gunnar@mi.fu-berlin.de, mbauer@inf.fu-berlin.de                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef LISAKLONGESTLAGRANGEPROBLEM_BOOST_H
#define LISAKLONGESTLAGRANGEPROBLEM_BOOST_H


#include <OpenMS/ANALYSIS/DENOVO/AntilopeLagrangeProblem.h>
#include <OpenMS/ANALYSIS/DENOVO/LagrangeProblemBase.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeSpectrumGraph.h>
//#include <Common/CommonTypes.h>


namespace OpenMS {


class DeNovoLagrangeProblemBoost : public LagrangeProblem {

  //typedefs
  typedef std::vector<DoubleReal> DVector;
  typedef SpectrumGraphSeqan::VertexDescriptor VertexDescriptor;
  typedef SpectrumGraphSeqan::EdgeDescriptor EdgeDescriptor;
  typedef std::set<VertexDescriptor> FilterNodeMap;
  typedef std::set<EdgeDescriptor> FilterEdgeMap;


  public:

    struct path_score_pair
    {
        std::vector<VertexDescriptor> path;
        float score;

        bool operator<(const path_score_pair &b) const
        {
          return (score > b.score);
        }
    };


    ///constructor
    DeNovoLagrangeProblemBoost(const SpectrumGraphSeqan *G_in):
                          LagrangeProblem(),
                          G(G_in),
                          forbidden_nodes_(),                          
                          edge_weights_(),
                          edge_weights_bk_(),
                          spanning_edges_(),
                          prefix_(),
                          lower_bound_(MINUS_INF)
    {      
      best_feasible_solution_.score = MINUS_INF;
      best_infeasible_solution_.score = PLUS_INF;

      seqan::resizeEdgeMap(G->graph, forbidden_edges_, std::vector<bool>(seqan::numEdges(G->graph), false));
      seqan::resizeVertexMap(G->graph, forbidden_nodes_, std::vector<bool>(seqan::numVertices(G->graph), false));

      resetWeights(true);
    }


    /// copy constructor
    DeNovoLagrangeProblemBoost(const DeNovoLagrangeProblemBoost &dnlp):
                          LagrangeProblem(dnlp),
                          G(dnlp.G),
                          edge_weights_(dnlp.edge_weights_),
                          edge_weights_bk_(dnlp.edge_weights_bk_),
                          forbidden_nodes_(dnlp.forbidden_nodes_),
                          forbidden_edges_(dnlp.forbidden_edges_),
                          spanning_edges_(dnlp.spanning_edges_),
                          prefix_(dnlp.prefix_),
                          lower_bound_(MINUS_INF)
    {      
      best_feasible_solution_.score = MINUS_INF;
      best_infeasible_solution_.score = PLUS_INF;
      resetWeights(false);
    }

    /// assignment operator
    DeNovoLagrangeProblemBoost operator= (const DeNovoLagrangeProblemBoost &dnlp)
    {
      LagrangeProblem::operator= (dnlp);
      G = dnlp.G;
      edge_weights_ = dnlp.edge_weights_;
      edge_weights_bk_ = dnlp.edge_weights_bk_;
      forbidden_nodes_ = dnlp.forbidden_nodes_;
      forbidden_edges_ = dnlp.forbidden_edges_;
      spanning_edges_ = dnlp.spanning_edges_;
      prefix_ = dnlp.prefix_;
      lower_bound_ = dnlp.lower_bound_;

      best_feasible_solution_.score = MINUS_INF;
      best_infeasible_solution_.score = PLUS_INF;
      resetWeights(false);
    }

    /// Destructor
    ~DeNovoLagrangeProblemBoost() {;}



    int EvaluateProblem( const DVector& Dual,
         list<int>& DualIndices,
         double& DualValue,
         double& PrimalValue,
         DVector& Subgradient,
         list<int>& SubgradientIndices,
         DVector& PrimalSolution,
         DVector& PrimalFeasibleSolution
         );


    int ComputeFeasibleSolution(
               DVector Dual,
               DVector Primal
               ){return 0;}


    void set_lower_bound(double LB_in){lower_bound_ = LB_in;}


    DoubleReal get_lower_bound(){ return lower_bound_;}

  private:

    //private default constructor
    DeNovoLagrangeProblemBoost(){;}

    //remove this. seems to be conly for testing reasons
    double lower_bound_;

    ///SpectrumGraph
    const SpectrumGraphSeqan *G;

    /// the best longest path of the primal problem
    path_score_pair best_feasible_solution_;

    /// best (lowest scoring) dual path
    path_score_pair best_infeasible_solution_;

    ///subgradient vector
    std::vector<double> subgradient_;

    // nodes that can not be selected for the path either because of Yen or because of branching
    seqan::String<bool> forbidden_nodes_;

    // edges that can not be selected for the path because of Yen
    seqan::String<bool> forbidden_edges_;

    // nodes that must be selected for the path because of branching
    std::set<VertexDescriptor> forced_nodes_;

    ///the last(k-1)th path
    path_score_pair prefix_;

    seqan::String<DoubleReal> edge_weights_;

    seqan::String<DoubleReal> edge_weights_bk_;

    std::vector<std::vector<EdgeDescriptor> > spanning_edges_;


  public:
    /// method computes the longest path
    void computeLongestPath(DoubleReal &dual_score, DoubleReal &primal_score, const DVector& Lambda);

    /// check wether found solution is feasible
    bool check_feasibility();

    /// method computes the score of the primal problem, 0 if infeasible
    double compute_primal_score();

    ///return the list of subgradient indices
    void computeSubgradientIndices(std::list<int> &subgradient_ind);

    //Getter Methods

    path_score_pair get_longest_path();

    double get_best_primal_score()
    {
      return best_feasible_solution_.score;
    }

    //TODO check whether here is really the right place. returning the path in terms of nodes should be sufficient
    //vector<float> get_best_solution();

    int num_of_duals()
    {
      return G->getNumberOfClusters();
    }

    void set_forbidden_edges(seqan::String<bool> &forb_edges_in)
    {
      forbidden_edges_ = forb_edges_in;
    }

    void setPrefix(const path_score_pair& prefix_in)
    {
      prefix_ = prefix_in;
      std::cout<<"Last Prefix Node:" << prefix_.path.back() <<std::endl;
    }

    void forbidConflictingNodes(VertexDescriptor v);

    void forceNode(VertexDescriptor node);

    void forbidNode(VertexDescriptor node);

    void getBestIllegalPath(std::vector<VertexDescriptor> &path)
    {
      path = best_infeasible_solution_.path;
    }

    void reset();

    void getViolatedClusters(std::vector<Size> &clusters, const std::vector<VertexDescriptor> &path);

    void updateWeights(const DVector& Lambda);

    void resetWeights(bool init);

    void forbidEdge(EdgeDescriptor edge);

    void forbidEdge(VertexDescriptor source, VertexDescriptor target);

    const SpectrumGraphSeqan& getSpectrumGraph()
    {
      return *G;
    }


}; //k_longest_lagrangeProblem
}//namespace

#endif
