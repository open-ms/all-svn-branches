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

#ifndef LISALAGRANGESOLVERBASE_H
#define LISALAGRANGESOLVERBASE_H

#include <Solver/Lagrange/Base/LagrangeCommon.h>
#include <Solver/Lagrange/Base/LagrangeSolverInterface.h>


namespace OpenMS {

  /**
   * LagrangeSolverBase defines the base class for all solvers of the
   * Lagrangian dual
   */
  class LagrangeSolverBase : public LagrangeSolverInterface {
  protected:
    DVector _dual;
    DVector _primal;
    DVector _multiplierLowerBound;
    DVector _multiplierUpperBound;
    double _currentLowerBound;
    double _currentUpperBound;
    double _bestLowerBound;
    double _bestUpperBound;

    int _verbose; // verbose output?

  public:

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


};

#endif
