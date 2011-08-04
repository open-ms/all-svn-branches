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

#ifndef LISALAGRANGEPROBLEM_H
#define LISALAGRANGEPROBLEM_H

// STL stuff
#include <list>
#include <vector>
#include <OpenMS/CONCEPT/Types.h>


using std::list;


namespace OpenMS {

  /**
   * LagrangeSolverInterface defines the interface for all solvers of
   * the Lagrangian dual
   */
  class LagrangeProblem {
  public:

    const static DoubleReal PLUS_INF;
    const static DoubleReal MINUS_INF;

    typedef std::vector<DoubleReal> DVector;

    virtual ~LagrangeProblem() {}

    /**
       @param Dual (IN) holds the multipliers
       @param DualIndices (IN) holds a list of the multipliers that got changed
       @param DualValue (OUT) holds the solution value for the relaxed problem
       @param PrimalValue (OUT) holds the solution value for the original problem, might be -inf if non-existent
       @param Subgradient (OUT) holds the subgradient at the current position
       @param SubgradientIndices (OUT) holds a list of subgradient indices for the _current_ iteration, where {abs(Subgradient[x]) = 1, x in SubgradientIndices}
       @param PrimalSolution (OUT) holds the primal solution for the current position
    */
    virtual int EvaluateProblem( const DVector& Dual,
				 list<int>& DualIndices,
				 double& DualValue,
				 double& PrimalValue,
				 DVector& Subgradient,
				 list<int>& SubgradientIndices,
				 DVector& PrimalSolution,
				 DVector& PrimalFeasibleSolution
				 ) =0;


    /**
       @param Dual (IN) holds the multipliers
       @param Primal (OUT) holds the feasible solution to the original problem
    */
    virtual int ComputeFeasibleSolution(
					DVector Dual,
					DVector Primal
					) = 0;

  };
}; // namespace lisa

#endif
