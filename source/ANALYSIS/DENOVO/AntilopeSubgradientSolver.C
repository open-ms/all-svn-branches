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

// STL stuff
#include <iostream>
#include <stdio.h>
#include <limits>

#include <OpenMS/ANALYSIS/DENOVO/AntilopeSubgradientSolver.h>
#include <OpenMS/ANALYSIS/DENOVO/LagrangeProblemBase.h>

//#include <OpenMS/CONCEPT/Exception.h>
//#include <cstdlib>


#define EPSILON 0.0001

using namespace std;


namespace OpenMS {


  void SubgradientSolver::InitProblem( int PrimalDimension, int DualDimension, DVector LBounds, DVector UBounds )
  {
    // initialize problem with the correct dimension, set upper and lower bounds
    // on the value of the multipliers
    _primal = DVector( PrimalDimension );
    _dual = DVector( DualDimension );
    _multiplierLowerBound = LBounds;
    _multiplierUpperBound = UBounds;

    // initialize the rest
    _currentUpperBound = LagrangeProblem::PLUS_INF;
    _currentLowerBound = LagrangeProblem::MINUS_INF;
    _bestUpperBound = LagrangeProblem::PLUS_INF;
    //cout << "setting u to " << _bestUpperBound << endl;
    _bestLowerBound = LagrangeProblem::MINUS_INF;

    // switch on signal stuff
    //if (signal(SIGINT, catchSignals) == SIG_ERR) {  cerr << "error catching signals" << endl; exit(-1); }

    // initialize with user-specific values
    //      _noOfIterations = defaults_.getValue("noofiterations");
    //      _noOfNondecreasingIterations = defaults_.getValue("noofnondecreasingiterations");
    //      _my = defaults_.getValue("my");
    //      _verbose = defaults_.getValue("verbosesolver");
    //      _utimeLimit = defaults_.getValue("utimelimit");
  }


  void SubgradientSolver::InitProblem( int PrimalDimension, int DualDimension )
  {
    // initialize problem with the correct dimension
    _primal = DVector( PrimalDimension );
    _dual = DVector( DualDimension );

    _multiplierLowerBound = DVector( DualDimension );
    _multiplierUpperBound = DVector( DualDimension );

    for( int j=0;j<DualDimension;j++)
    {
      // initialize lower and upper bounds appropriately
      _multiplierLowerBound[j] = LagrangeProblem::MINUS_INF;
      _multiplierUpperBound[j] = LagrangeProblem::PLUS_INF;
    }


    // initialize the rest
    _currentUpperBound = LagrangeProblem::PLUS_INF;
    _currentLowerBound = LagrangeProblem::MINUS_INF;
    _bestUpperBound = LagrangeProblem::PLUS_INF;
    _bestLowerBound = LagrangeProblem::MINUS_INF;

    // switch on signal stuff
    //if (signal(SIGINT, catchSignals) == SIG_ERR) {  cerr << "error catching signals" << endl; exit(-1); }


    //    // initialize with user-specific values
    //    _noOfIterations = defaults_.getValue("noofiterations");
    //    _noOfNondecreasingIterations = defaults_.getValue("noofnondecreasingiterations");
    //    _my = defaults_.getValue("my");
    //    _verbose = defaults_.getValue("verbosesolver");
    //    _utimeLimit = defaults_.getValue("utimelimit");
  }


  /**
   * SubgradientSolver specific problem initialization
   */
  void SubgradientSolver::InitProblem(int PrimalDimension,int DualDimension,DVector LBounds,DVector UBounds,
                                      int NoOfIterations,int NoOfNondecreasingIterations, double My)
  {
    _primal = DVector( PrimalDimension );
    _dual = DVector( DualDimension );
    _multiplierLowerBound = LBounds;
    _multiplierUpperBound = UBounds;

    // initialize the rest
    _currentUpperBound = LagrangeProblem::PLUS_INF;
    _currentLowerBound = LagrangeProblem::MINUS_INF;
    _bestUpperBound = LagrangeProblem::PLUS_INF;
    _bestLowerBound = LagrangeProblem::MINUS_INF;
    _noOfIterations = NoOfIterations;
    _noOfNondecreasingIterations = NoOfNondecreasingIterations;
    _my = My;


    //    // initialize with user-specific values
    //   _verbose = defaults_.getValue("verbosesolver");
    //   _utimeLimit = defaults_.getValue("utimelimit");


  }

  /**
   * SubgradientSolver specific problem initialization
   */
  void SubgradientSolver::InitProblem(int PrimalDimension,int DualDimension,int NoOfIterations,
                                      int NoOfNondecreasingIterations, double My)
  {
    // initialize problem with the correct dimension
    _primal = DVector( PrimalDimension );
    _dual = DVector( DualDimension );

    _multiplierLowerBound = DVector( DualDimension );
    _multiplierUpperBound = DVector( DualDimension );

    for( int j=0;j<DualDimension;j++)
    {
      // initialize lower and upper bounds appropriately
      _multiplierLowerBound[j] = LagrangeProblem::MINUS_INF;
      _multiplierUpperBound[j] = LagrangeProblem::PLUS_INF;
    }

    // initialize the rest
    _currentUpperBound = LagrangeProblem::PLUS_INF;
    _currentLowerBound = LagrangeProblem::MINUS_INF;
    _bestUpperBound = LagrangeProblem::PLUS_INF;
    _bestLowerBound = LagrangeProblem::MINUS_INF;

    _noOfIterations = NoOfIterations;
    _noOfNondecreasingIterations = NoOfNondecreasingIterations;
    _my = My;

    // switch on signal stuff
    //if (signal(SIGINT, catchSignals) == SIG_ERR) {  cerr << "error catching signals" << endl; exit(-1); }

    //    // initialize with user-specific values
    //       _verbose = defaults_.getValue("verbosesolver");
    //       _utimeLimit = defaults_.getValue("utimelimit");
  }



  /**
   * LagrangeSolverInterface specific methods
   */
  typedef std::vector<DoubleReal> DVector;

  int SubgradientSolver::SetMultiplierLowerBound( DVector LBounds ) { _multiplierLowerBound = LBounds; }

  int SubgradientSolver::SetMultiplierUpperBound( DVector UBounds ) { _multiplierUpperBound = UBounds; }

  SubgradientSolver::DVector SubgradientSolver::GetFinalSolution( void ) const { return _dual; }

  DVector SubgradientSolver::GetFinalHeuristicSolution( void ) const { return _primal; }

  DVector SubgradientSolver::GetDualVariables( void ) const { return _dual; }

  DVector SubgradientSolver::GetPrimalVariables( void ) const { return _primal; }

  double SubgradientSolver::GetCurrentUpperBound( void ) const { return _currentUpperBound; }

  double SubgradientSolver::GetCurrentLowerBound( void ) const { return _currentLowerBound; }

  int SubgradientSolver::GetPrimalDimension( void ) const { return _primal.size(); }

  int SubgradientSolver::GetDualDimension( void ) const { return _dual.size(); }


  /**
   * SugradientSolver specific methods
   */

  /*
   *  If <code>EvaluateProblem</code> returns a value different from
   *  0 (:=LagrangeSolverBase::UNSOLVED), the iteration index returned by this
   *  function might be different from the return value of method
   *  <code>GetNoOfIterations</code>, since optimization is finished at that
   *  moment.
   *  @return index of last iteration
   */
  int SubgradientSolver::GetLastIterationIdx() const { return _noOfLastIteration; }

  /**
   * Solve iteratively the dual problem by subgradient optimization
   * @return 0 if successful, 1 otherwise
   */
  SubgradientSolver::SolverProgress SubgradientSolver::Solve( LagrangeProblem& l )
  {
    // all variables initialized?
    if( _sanityCheck() == false )
      return INSANE;

    int evalProb_returnVal = 0;

    int noNondecreasingRounds = -1;

    list<int> dualIndices;
    list<int> subgradientIndices;

    DVector subgradient( _dual.size() );
    DVector primalVector( _primal.size() );
    DVector primalfeasibleSolution( _primal.size() );

    std::numeric_limits<double> dlimits;

    for( int i=0;i<_noOfIterations; i++ )
    {
      _noOfLastIteration=i;

      //if (ctrlC_count >= 1) return UNSOLVED_INTERRUPT;
      //if (get_ctrlC_counter() >= 1)
      //return UNSOLVED_INTERRUPT;

      if( _verbose == 1 ) {
        cout << "(" << i << ")\tbest: " << _bestUpperBound << "\t/" << _bestLowerBound << "\tcurrent: ";
        cout.flush();
      }

      try {

        if( _utimeLimit+1 != 0 )
          _curIter.start(); // TIC

        // guwek: removed the dependency on the return code from
        // evaluateProblem. i did not understand.

        evalProb_returnVal = l.EvaluateProblem( _dual,
                                                dualIndices,
                                                _currentUpperBound,
                                                _currentLowerBound,
                                                subgradient,
                                                subgradientIndices,
                                                primalVector,
                                                primalfeasibleSolution
                                                );
        if( _utimeLimit+1 != 0 )
          _curIter.stop(); // TOC

        if( _verbose == 1 ) {
          cout << _currentUpperBound << "/\t" << _currentLowerBound//;
              << "\t(" << subgradientIndices.size() << ")";
        }
        // compare upper and lower bound
        if( _currentUpperBound < _bestUpperBound )
        {
          _bestUpperBound = _currentUpperBound;
          noNondecreasingRounds = -1;
        }

        if( _currentLowerBound > _bestLowerBound )
        {
          _bestLowerBound = _currentLowerBound;
          noNondecreasingRounds = -1;
          if( _verbose == 1 ) {
            cout << "(*)\n";
          }
        } else {
          if( _verbose == 1 )
            cout << "\n";
        }

        // now check whether we already reached an optimal solution
        if( (_bestUpperBound - _bestLowerBound)< EPSILON ){
          cout<<"case 1"<<endl;
          cout<<"SOLVED IN ITERATION: "<<i+1<<endl;
          return SOLVED;      }



        // increase the number of iterations in every single
        // iteration, since we're initializing noNonDecreasingRounds
        // with -1
        noNondecreasingRounds++;

        if( noNondecreasingRounds == _noOfNondecreasingIterations )
        {
          // there was nothing going on in the last couple of iterations,
          // half _my therefore
          _my /= 2;

          if (_verbose == 1) { cout << "Setting my to " << _my << endl; }
          noNondecreasingRounds = 0;
        }


        // calculate number of subgradients
        int noOfSubgradients = subgradientIndices.size();

#ifdef VERBOSE_OUTPUT
        cout << "number of subgradients = " << noOfSubgradients << endl;
#endif

        double temp = 0;
        list<int>::iterator sg_iter;
        for(sg_iter = subgradientIndices.begin(); sg_iter != subgradientIndices.end(); sg_iter++){
          temp += (subgradient[*sg_iter]*subgradient[*sg_iter]);
        }
        double ss= _my* ((_bestUpperBound - _bestLowerBound)/sqrt(temp));
        //hack sandro
        ss = min(ss,2.0);
        cout << "LagrangeITERATION: = " << i << endl;

#ifdef VERBOSE_OUTPUT
        cout << "stepsize = " << ss << endl;
#endif


        // stop optimizing when the stepsize falls below the machine
        // accuracy ...

        // guwek [apr 07]: replaced the folowing by the code
        // including 'changed'. see below. I have ssen instances that got solved
        // with ss < epsilon, so I stop know if all the multipliers stay the same.

        if (ss < dlimits.epsilon()) return UNSOLVED_MU;

        // adapt dual variables according to stepsize
        std::list<int>::iterator it;

        bool changed = false;
        int count=0;

        for( it=subgradientIndices.begin();it!=subgradientIndices.end();it++) {
          count++;
          double oldValue = _dual[*it];
          double newValue = _dual[*it] - ss*subgradient[*it];
          //cout<<"old value: "<<oldValue<<" new value: "<< newValue<<endl;
          if( newValue < _multiplierLowerBound[*it] ) {
            _dual[*it] = _multiplierLowerBound[*it];
          } else if( newValue > _multiplierUpperBound[*it] ) {
            _dual[*it] = _multiplierUpperBound[*it];
          } else {
            _dual[*it] = newValue;
          }
          /*	    if (_dual[*it] != oldValue) {
            changed = true;
            //cout<<"changed set true here! "<<*it<<" to: "<<_dual[*it]<<endl;
        }
*/        
        } // for
        /*	  if (!changed){
          cout<<"hier gehts raus!"<<endl;
          return UNSOLVED_MU;
      }
*/      
			} catch (Exception::BaseException e) {
				cerr << "SubgradientSolver::EvaluateProblem " << e.what() << endl;
			}

#ifdef VERBOSE_OUTPUT
			cout << "subgradient.size  = " << subgradient.size() << endl;
			cout << "primalVector.size = " << primalVector.size() << endl;
			cout << "primalfeasibleSolution.size = " << primalfeasibleSolution.size() << endl;
			cout << "_multiplierLowerBound.size = " << _multiplierLowerBound.size() << endl;
			cout << "_multiplierUpperBound.size = " << _multiplierUpperBound.size() << endl;
#endif

			// assign the indices of the last iteration to the new one
			dualIndices = subgradientIndices;
			subgradientIndices.clear();
			// reset subgradient vector
			for( std::list<int>::iterator it = dualIndices.begin(); it != dualIndices.end(); ++it )
				subgradient[*it] = 0.0;

			//primalVector.clear();
			//primalVector.resize( _primal.size() );
			//primalfeasibleSolution.clear();
			//primalfeasibleSolution.resize( _primal.size() );  //TODO check whether this is correct to be commented

			if (i == _noOfIterations - 1)
				return UNSOLVED_ITERATIONS;

			if ( _utimeLimit+1 != 0 ) // _utimeLimit != -1
			{
        _totalIters_time += _curIter.getUserTime();
        if( _totalIters_time >= _utimeLimit ) // cast should be stable as long as we do not want to compute until "St. Nimmerleins Tag"
					return UTIME_LIMIT;
			}
		}
		cout<<"case2"<<endl;
    return SOLVED;
  }

  /**
   * Perform one single subgradient step
   * @return 0 if successful, 1 otherwise
   */
  int SubgradientSolver::SolveSingleIteration( LagrangeProblem& l )
  {
    return 0;
  }


  double SubgradientSolver::_calculateStepsize(double upper, double lower, int noOfSubgradients)
  {
    return _my*( (upper-lower)/noOfSubgradients );
  }

  double SubgradientSolver::_calculateStepsize()
  {
    // _noOfLastIterations+1 because we start from 0 for the number of iterations
    return (1/(double)(_noOfLastIteration+1));
  }

  bool SubgradientSolver::_sanityCheck( void )
  {
    if( (_my == -1.0) ||
        (_noOfIterations == -1) ||
        (_noOfNondecreasingIterations == -1) ||
        (_verbose == -1 )
      ) {

      throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Variables not initialized. Check your parameter file.");
    }

    return true;
  }

  void SubgradientSolver::setDefaultParams_()
  {
    defaults_.setValue("noofiterations", 200, "Number of iterations for the subgradient method");
    defaults_.setMinInt("noofiterations", 0);

    defaults_.setValue("noofnondecreasingiterations",20,"number of non-decreasing iteration before abort");
    defaults_.setMinInt("noofnondecreasingiterations", 1);

    defaults_.setValue("my", 1.0, "The stepsize in the subgradient method");
    defaults_.setMinFloat("my", 0.05);

    defaults_.setValue("verbosesolver",1,"set verbose output of the solver");
    defaults_.setMaxInt("verbosesolver",1);
    defaults_.setMinInt("verbosesolver",0);

    defaults_.setValue("utimelimit", -1,"The maximal time allowed before optimization is abortet. choose -1 to deactivate");
    defaults_.setMinInt("utimelimit", -1);

    defaultsToParam_();
  }

}
