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

//#include <iostream>
//#include <ctime>


#include <OpenMS/ANALYSIS/DENOVO/AntilopeILP.h>
// ILOG stuff
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloalg.h>


#define debug 1
//#define PAIRWISE
#undef PAIREWISE
using namespace std;

namespace OpenMS {

//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------Constructors------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------

  //default constructor
  de_novo_ILP::de_novo_ILP():DefaultParamHandler("de_novo_ILP")
  {
    defaults_.setValue("suboptimals",(UInt)20,"The number of suboptimal solutions to be computed");
    defaults_.setValue("parenttolerance",2.0,"The maximum mass difference between predicted peptide and parent mass in query spectrum");
    defaults_.setMinFloat("parenttolerance",0.5);

    defaultsToParam_();
  }




  //copy constructor
  de_novo_ILP::de_novo_ILP(const de_novo_ILP & ilp):DefaultParamHandler(ilp)
  {
  }




  //assignment Operator
  de_novo_ILP & de_novo_ILP::operator=(const de_novo_ILP& in)
  {
    if(this !=&in)
    {
      DefaultParamHandler::operator=(in);
    }
    return *this;
  }


//--------------------------------------------------------------------------------------------------------------
//---------------------------------------Method  computeCandidates----------------------------------------------
//--------------------------------------------------------------------------------------------------------------
void de_novo_ILP::computeCandidates(const SpectrumGraphSeqan &graph, CandSetInt &cands, DRealVec &scores)
{

  //parameter extraction
  UInt num_suboptimals=param_.getValue("suboptimals");


    std::cout << "Algorithm: de_novo_ILP_solver" << std::endl;

    //CandSetReal result_edge_masses_;
    //CandSetInt result_edge_masses_int_;
    //DRealVec result_scores_;

    //begin 1.zeitmessung!!
    double time1 = 0.0, time2 = 0.0, tstart1, tstart2; // time measurment variables
    tstart1 = clock(); // start
    //---------------------------------------



    //---------------------------------------------------------------------------------------------------------
    //------------------------------------------Building ILP Model---------------------------------------------
    //---------------------------------------------------------------------------------------------------------

    IloEnv env; // get an ILOG environment (has something to do with the license)

    //ILP variables for directed edges
    IloModel M(env); // get an empty ILOG model
  
    UInt num_of_nodes = graph.size();
    std::vector <IloBoolVar> x;
    x.reserve(num_of_nodes);
    typedef std::list<std::pair<Size, Size> > TOutEdgeList;
    typedef TOutEdgeList::const_iterator TOutEdgeListIter;
    std::vector<TOutEdgeList> adj_list(num_of_nodes);
    std::vector<TOutEdgeList> adj_list_rev(num_of_nodes);
  
    //std::map<std::pair<UInt, UInt>, IloBoolVar> x;
  


    //------------------------------------objective function-----------------------------------------------

    UInt num_variables = 0;
    IloExpr obj(env);
    IloExpr delta_const(env);

    std::cout << "NUMOFNODES " << num_of_nodes << std::endl;
    for (UInt i = 0; i < num_of_nodes - 1; i++)
    {
      for (UInt j = i + 1; j < num_of_nodes; j++)
      {
        if (graph.findDirectedEdge(i, j))
        {
          x.push_back(IloBoolVar(env));
          adj_list[i].push_back(make_pair(j, num_variables));
          adj_list_rev[j].push_back(make_pair(i, num_variables));
          
          num_variables++;
          std::cout<<"edge_mass:: "<<graph.getRealEdgeMass(i, j)<<std::endl;
          obj += x.back() * -(graph.getEdgeWeight(i, j));
          //cout<<"create x-variable"<<i<<":"<<j<<": "<<graph.getEdgeWeight(i,j)<<endl;
        }
      }
    }
    //adding obejctive function in order to maximize
    M.add(IloMaximize(env, obj));
    obj.end();

    if (debug) cout << "objective function added" << endl;

    //------------------------------------adding constraints--------------------------------------------------
    //

    //constraint that garantuees, that the total mal error is under a certain threshold, here 2 Dalton
    //   cout<<"PARENTMASS in ILP "<< graph->get_parent_mass()-19<<endl;

    //   M.add(delta_const- (G->get_parent_mass()-19)<=2);
    //   M.add(delta_const- (G->get_parent_mass()-19)>=-2);
    //   delta_const.end();
    //   if (debug) cout<< "delta const added"<<endl;

    // constraint forbidding selection of two complementary nodes except source and sink (node 0 and last node)
    // the selection of a node is here modelled by the selection of incomiung edges
#ifdef PAIRWISE
    for (UInt i = 1; i < num_of_nodes - 2; i++)
    {
      for (UInt j = i + 1; j < num_of_nodes - 1; j++)
      {
        if (graph.hasUndirectedEdge(i, j))
        {
          IloExpr compl_const(env);
          
          for (TOutEdgeListIter in_edge_it = adj_list_rev[i].begin();
               in_edge_it != adj_list_rev[i].end();
               ++in_edge_it)
          {
            compl_const += x[in_edge_it->second]; //summing up all incoming edges for i
          }
          
          for (TOutEdgeListIter in_edge_it = adj_list_rev[j].begin();
               in_edge_it != adj_list_rev[j].end();
               ++in_edge_it)
          {
            compl_const += x[in_edge_it->second]; //summing up all incoming edges for j
          }          

          std::cout << "added u_edge: " << i << " --> " << j << std::endl;
          M.add(compl_const <= 1); // at most one of all complementary nodes can be in the solution
          compl_const.end();
        }
      }
    }
#else

  typedef SpectrumGraphSeqan::VertexDescriptor VertexDescriptor;
  for (UInt c = 0; c < graph.getNumberOfClusters(); ++c)
  {
    const std::vector<VertexDescriptor>& cluster = graph.getCluster(c);
    IloExpr compl_const(env);
    
    for (UInt i = 0; i < cluster.size(); ++i)
    {
        VertexDescriptor v = cluster[i];

        for (TOutEdgeListIter in_edge_it = adj_list_rev[v].begin(); in_edge_it != adj_list_rev[v].end(); ++in_edge_it)
        {
          compl_const += x[in_edge_it->second]; //summing up all incoming edges for v
        }
    }
    std::cout << "clique: " << c << std::endl;
    M.add(compl_const <= 1); // at most one of all complementary nodes can be in the solution
    compl_const.end();
  }
#endif

    if (debug) cout << "constraint 1 done" << endl;

    // Constraint guarantees for input = output for each node (FLOW CONSERVATION)
    for (UInt i = 1; i < num_of_nodes - 1; i++)
    {
      // for all nodes but source and sink
      IloExpr input_minus_output(env);

      for (TOutEdgeListIter in_edge_it = adj_list_rev[i].begin();
           in_edge_it != adj_list_rev[i].end();
           ++in_edge_it)
      {
        input_minus_output += x[in_edge_it->second]; // adding inputs
      }

      for (TOutEdgeListIter out_edge_it = adj_list[i].begin();
           out_edge_it != adj_list[i].end();
           ++out_edge_it)
      {
        input_minus_output -= x[out_edge_it->second]; // subtracting outputs
      }

      M.add(input_minus_output == 0);
      input_minus_output.end();
    }

    if (debug) cout << "constraint 2 done" << endl;

    // constraints for source and sink nodes
    IloExpr source_out(env); //all edges out of the source

    IloExpr sink_in(env); //all edges ending in the sink

    for (TOutEdgeListIter in_edge_it = adj_list_rev[num_of_nodes - 1].begin();
         in_edge_it != adj_list_rev[num_of_nodes - 1].end();
         ++in_edge_it)
    {
      sink_in += x[in_edge_it->second];
    }
  
    for (TOutEdgeListIter out_edge_it = adj_list[0].begin();
         out_edge_it != adj_list[0].end();
         ++out_edge_it)
    {
      source_out += x[out_edge_it->second];
    }

    M.add(source_out == 1); //exactly ONE edge leaving the source node
    M.add(sink_in == 1); //exactly ONE edge entering the sink node
    source_out.end();
    sink_in.end();

    if (debug) cout << "constraint 3 done" << endl;

    //Ende1 zeitmessung!!
    time1 += clock() - tstart1;
    if (debug) std::cout << "  Zeit vor cplex = " << time1 / CLOCKS_PER_SEC << " sec." << std::endl;
    if (debug) std::cout << "CLOCKS_PER_SEC" << CLOCKS_PER_SEC << std::endl;


    //---------------------------------------------------------------------------------------------------------
    //--------------------------------------Solving and querying result----------------------------------------
    //---------------------------------------------------------------------------------------------------------

    //put model into a cplex object
    IloCplex cplex(M);

    std::cout << "model created" << std::endl;
    cplex.exportModel("de_novo_ILP.lp");
    //cplex.setParam(IloCplex::EpInt,0.0);

    DoubleReal time_compl_solv = 0;
    DoubleReal tstart5 = clock();
    // start

    std::vector<UInt> result_Vec;
    std::vector<UInt> result_vec_int;

    //TODO number of suboptimal solutions as param
    //loop computes the r optimal + suboptimal solutions by iterativlely solving and modifying the ILP
    for (UInt r = 0; r < num_suboptimals; r++)
    {
      cout << "ITERATION " << r << endl;

      //Begin zweite Zeitmessung
      tstart2 = clock(); // start

      //cplex.clear()
      //solve model
      if (!cplex.solve())
      {
        env.error() << "Fehler beim Loesen" << endl;
        //throw(-1);
        break;
      }
      else
      {
        //Ende2 zeitmessung!!
        time2 = clock() - tstart2;
        time_compl_solv = (clock() - tstart5) / CLOCKS_PER_SEC;
        if (debug) cout << "  Zeit fuer solving = " << time2 / CLOCKS_PER_SEC << " sec." << "::::" << time2 << endl;
        //if (debug) cout << "  Zeit fuer solving TTT = " << time_compl_solv/CLOCKS_PER_SEC << " sec." <<"::::"<<time_compl_solv<<endl;
        if (debug) cout << "  Zeit fuer solving TTT = " << time_compl_solv << " sec." << "::::" << time_compl_solv << endl;

        cout << "after timing" << endl;
        cplex.out() << "solution status= " << cplex.getStatus() << endl;

        //the update for the next iteration (forbid the actual solution)
        IloExpr update(env); // used to exclude current path from future iterations

        IloInt count_true_edges = 0;

        result_Vec.clear();
        result_vec_int.clear();

        DoubleReal total_score = 0.0;
        Size source_v = 0;

        while (source_v != num_of_nodes-1)
        {
          for (TOutEdgeListIter out_edge_it = adj_list[source_v].begin();
               out_edge_it != adj_list[source_v].end();
               ++out_edge_it)
          {
            if (cplex.getValue(x[out_edge_it->second]) >= 1.0 - cplex.getParam(IloCplex::EpInt))
            {
              Size target_v = out_edge_it->first;
              update += x[out_edge_it->second];
              ++count_true_edges;
              result_vec_int.push_back(graph.getIntEdgeMass(source_v, target_v));
              total_score += graph.getEdgeWeight(source_v, target_v);
              source_v = target_v;
              break;
            }
          }
        }

        cands.push_back(result_vec_int);

        if(r > 0 and scores.back()<cplex.getObjValue())
        {
            std::cout<<"Warning, better score in later iteration!"<<std::endl;
        }

        scores.push_back(cplex.getObjValue());

        cout << "score " << cplex.getObjValue() << endl;

        //add update constraint to model (forbid actual solution)
        M.add(update <= count_true_edges - 1);
        update.end();
      }
    }

    double time5 = clock() - tstart5;
    if (debug) cout << "  TOTAL Zeit fuer solving = " << time5 / CLOCKS_PER_SEC << " sec." << endl;

    cout << "after getObjValue" << endl;

    //clear model
    x.clear();
    cplex.end();

    //compiler problems
    //env.end();
  }

}// namespace




