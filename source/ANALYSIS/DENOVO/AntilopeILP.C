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
    std::map<std::pair<UInt, UInt>, IloBoolVar> x;
    UInt num_of_nodes = graph.size();


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
          num_variables++;
          x[std::make_pair(i, j)] = IloBoolVar(env);
          std::cout<<"edge_mass:: "<<graph.getRealEdgeMass(i, j)<<std::endl;
          obj += x[std::make_pair(i, j)] * (graph.getEdgeWeight(i, j));
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
    for (UInt i = 1; i < num_of_nodes - 2; i++)
    {
      for (UInt j = i + 1; j < num_of_nodes - 1; j++)
      {
        if (graph.hasUndirectedEdge(i, j))
        {
          IloExpr compl_const(env);

          for (UInt k = 0; k < i; k++)
          {
            if (graph.findDirectedEdge(k, i)) compl_const += x[std::make_pair(k, i)]; //summing up all incoming edges for i
          }
          for (UInt k = 0; k < j; k++)
          {
            if (graph.findDirectedEdge(k, j)) compl_const += x[std::make_pair(k, j)]; //summing up all incoming edges for j
          }

          std::cout << "added u_edge: " << i << " --> " << j << std::endl;
          M.add(compl_const <= 1); // at most one of all complementary nodes can be in the solution
          compl_const.end();
        }
      }
    }

    if (debug) cout << "constraint 1 done" << endl;

    // Constraint guarantees for input = output for each node (FLOW CONSERVATION)
    for (UInt i = 1; i < num_of_nodes - 1; i++)
    {
      // for all nodes but source and sink
      IloExpr input_minus_output(env);

      for (UInt k = 0; k < i; k++)
      {
        if (graph.findDirectedEdge(k, i))
          input_minus_output += x[std::make_pair(k, i)]; // adding inputs
      }

      for (UInt l = i + 1; l < num_of_nodes; l++)
      {
        if (graph.findDirectedEdge(i, l))
          input_minus_output -= x[std::make_pair(i, l)]; // subtracting outputs
      }

      M.add(input_minus_output == 0);
      input_minus_output.end();
    }

    if (debug) cout << "constraint 2 done" << endl;

    // constraints for source and sink nodes
    IloExpr source_out(env); //all edges out of the source

    IloExpr sink_in(env); //all edges ending in the sink

    for (UInt i = 0; i < num_of_nodes; i++)
    {
      if (graph.findDirectedEdge(0, i))
      {
        source_out += x[std::make_pair(0, i)];
      }

      if (graph.findDirectedEdge(i, num_of_nodes - 1))
      {
        sink_in += x[std::make_pair(i, num_of_nodes - 1)];
      }
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
        IloExpr update(env);

        IloInt count_true_edges = 0;

        result_Vec.clear();
        result_vec_int.clear();

        DoubleReal total_score = 0.0;
        UInt from = 0;
        UInt to = 1;

        while (to != num_of_nodes)
        {
          if (graph.findDirectedEdge(from, to))
          {
            if (cplex.getValue(x[std::make_pair(from, to)]) >= 1.0 - cplex.getParam(IloCplex::EpInt))
            {
              cout << from << "!:!" << to << endl;

              update += x[std::make_pair(from, to)];
              ++count_true_edges;

              //result_Vec.push_back(graph.getRealEdgeMass(from, to));
              result_vec_int.push_back(graph.getIntEdgeMass(from, to));

              total_score += graph.getEdgeWeight(from, to);
              from = to;
            }
          }
          ++to;
        }

        //cands.push_back(result_Vec);
        cands.push_back(result_vec_int);
        //result_edge_masses_int_.push_back(result_Vec_int);

        if(r>0 and scores.back()<cplex.getObjValue())
        {
            std::cout<<"Error, better score in later iteration!"<<std::endl;
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




