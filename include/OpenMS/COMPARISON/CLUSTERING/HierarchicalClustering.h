// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut, Steffen Sass $
// --------------------------------------------------------------------------


#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H

namespace OpenMS
{
  /**
   * @brief generic n-dimensional hierarchical clustering
   */
  template <typename PointInfo, int Dim = 2>
  class HierarchicalClustering
  {
    public:
      typedef boost::array<DoubleReal, Dim> Point;

      /**
       * @brief Subset of points.
       *
       * Describes a set of points in the grid.
       */
      typedef typename boost::unordered_multimap<Point, PointInfo> Subset;
      typedef HashGrid<Subset, Dim> Grid;

      Grid grid;

      HierarchicalClustering(const Point &max_delta)
        : grid(max_delta)
      { }

      /**
       * @brief Insert new Point into grid.
       * @param d Point to insert.
       * @param info Associated caller specified info.
       * @return iterator to inserted subset.
       */
      typename Grid::local_iterator insertPoint(const Point &d, const PointInfo &info)
      {
        typename Grid::local_iterator it = insertSubset(d);
        it->second.insert(std::make_pair(d, info));
        return it;
      }

    protected:
      /**
       * @brief Insert new Subset into grid.
       * @param d Point to insert.
       * @return iterator to inserted subset.
       */
      typename Grid::local_iterator insertSubset(const Point &d)
      {
        return grid.insert(std::make_pair(d, Subset()));
      }
  };
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H */
