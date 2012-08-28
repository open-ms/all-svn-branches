/***************************************************************************
 *   Copyright (C) 2006 by Gunnar W. Klau, Markus Bauer, Patrick May       *
 *   gunnar@mi.fu-berlin.de, mbauer@inf.fu-berlin.de, patrick.may@zib.de   *
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

#ifndef LISALAGRANGESUBGRADIENTSOLVER_H
#define LISALAGRANGESUBGRADIENTSOLVER_H

#include <OpenMS/ANALYSIS/DENOVO/AntilopeLagrangeProblem.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/StopWatch.h>


#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

// including stl stuff
#include <string>
#include <cstdlib>

using namespace std;

namespace OpenMS
{

  /**
   * SubgradientSolver implements the subgradient method for solving the
   * Lagrangian (dual) problem
   */
  class SubgradientSolver: public DefaultParamHandler
  {

    typedef std::vector<DoubleReal> DVector;

  public:

    enum SolverProgress
    {
      SOLVED,               /** solved! */
      UNSOLVED,             /** unsolved for some reason */
      UNSOLVED_INTERRUPT,   /** unsolved (caught interrupt) */
      UNSOLVED_ITERATIONS,  /** unsolved (no. of iterations exceeded) */
      UNSOLVED_MU,          /** unsolved (mu below precision) */
      WEIRD,                /** weird (this should not happen) */
      INSANE,               /** initialization does not look good */
      UTIME_LIMIT           /** unsolved (aborted due to exceeded user-time limit)*/
    };

    protected:
      DVector _dual;
      DVector _primal;
      DVector _multiplierLowerBound;
      DVector _multiplierUpperBound;
      double _currentLowerBound;
      double _currentUpperBound;
      double _bestLowerBound;
      double _bestUpperBound;

      double _my;
      double _stepsize;

      int _noOfIterations;
      int _noOfNondecreasingIterations;
      int _noOfLastIteration;

      int _verbose; // verbose output?

      int _utimeLimit; // limit for the user time
      StopWatch _curIter;
      double _totalIters_time;

      /// Calculates the current stepsize given the vector of subgradients
      /// and _dual
      double _calculateStepsize(double upper, double lower, int noOfSubgradients);

      /// Calculates the stepsize according to the harmonic series, i.e. we have
      /// _provable_ convergence
      double _calculateStepsize();

      /// Updates the values of the dual variables given the stepsize as
      /// computed by _calculateStepsize
      int _updateDualVariables(void);

      /// Before starting the optimization process check whether all required
      /// variables have been initialized
      bool _sanityCheck(void);

      ///sets the default parameters
      void setDefaultParams_();

      void updateMembers_()
      {
        _noOfIterations = param_.getValue("noofiterations");
        _noOfNondecreasingIterations = param_.getValue("noofnondecreasingiterations");
        _my = param_.getValue("my");
        _verbose = param_.getValue("verbosesolver");
        _utimeLimit = param_.getValue("utimelimit");
      }

    public:
      /// default constructor --- problem initialization is done in InitProblem
      SubgradientSolver(void) :
        DefaultParamHandler("SubgradientSolver")
      {
        setDefaultParams_();
        defaultsToParam_();
      }

      ~SubgradientSolver(void)
      {
      }
      ;

      /**
       * methods declared virtual in LagrangeSolverInterface.h
       */
      void InitProblem(int PrimalDimension, int DualDimension, DVector LBounds, DVector UBounds);
      void InitProblem(int PrimalDimension, int DualDimension);
      int SetMultiplierLowerBound(DVector LBounds);
      int SetMultiplierUpperBound(DVector UBounds);
      SolverProgress Solve(LagrangeProblem& l);

      int SolveSingleIteration(LagrangeProblem& l);
      DVector GetFinalSolution(void) const;
      DVector GetFinalHeuristicSolution(void) const;
      DVector GetDualVariables(void) const;
      DVector GetPrimalVariables(void) const;
      double GetBestUpperBound() const
      {
        return _bestUpperBound;
      }
      double GetBestLowerBound() const
      {
        return _bestLowerBound;
      }
      double GetCurrentUpperBound(void) const;
      double GetCurrentLowerBound(void) const;
      int GetPrimalDimension(void) const;
      int GetDualDimension(void) const;

      /**
       * SubgradientSolver specific methods
       */
      void InitProblem(int PrimalDimension, int DualDimension, DVector LBounds, DVector UBounds, int NoOfIterations, int NoOfNondecreasingIterations, double My);

      void InitProblem(int PrimalDimension, int DualDimension, int NoOfIterations, int NoOfNondecreasingIterations, double My);

      double GetCurrentStepsize() const
      {
        return _stepsize;
      }
      ;

      double GetMy() const
      {
        return _my;
      }
      ;
      void SetMy(double m)
      {
        _my = m;
      }
      ;

      int GetNoOfIterations() const
      {
        return _noOfIterations;
      }
      ;
      void SetNoOfIterations(int n)
      {
        _noOfIterations = n;
      }
      ;

      int GetNoOfNondecreasingIterations() const
      {
        return _noOfNondecreasingIterations;
      }
      ;
      void SetNoOfNondecreasingIterations(int n)
      {
        _noOfNondecreasingIterations = n;
      }
      ;

      int GetVerbose() const
      {
        return _verbose;
      }
      void SetVerbose(int v)
      {
        _verbose = v;
      }

      /*
       *  If <code>EvaluateProblem</code> returnes a value different from
       *  0 (:=LagrangeSolverBase::UNSOLVED), the iteration index returned by this
       *  function might be different from the return value of method
       *  <code>GetNoOfIterations</code>, since optimisation is finished at that
       *  moment.
       *  @return index of last iteration
       */
      int GetLastIterationIdx() const;


      /**
       * @brief Interprets given solver progress.
       *
       * @param p  a solver progress enum
       */
      std::string interpret(SolverProgress p) {
        if (p == SOLVED) return "solved";
        if (p == UNSOLVED) return "unsolved (unknown reason)";
        if (p == UNSOLVED_INTERRUPT) return "unsolved (caught interrupt)";
        if (p == UNSOLVED_ITERATIONS) return "unsolved (no. of iterations exceeded)";
        if (p == UNSOLVED_MU) return "unsolved (mu below precision)";
        if (p == UTIME_LIMIT) return "unsolved (aborted due to exceeded user-time limit)";
        if (p == WEIRD) return "weird (this should not happen)";
        if (p == INSANE) return "initialization does not look good";
        return "unknown solver progress code";
      }
  };
}

#endif
