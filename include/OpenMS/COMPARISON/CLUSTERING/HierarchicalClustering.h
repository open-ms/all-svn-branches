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


#ifndef OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>

namespace OpenMS
{
  template <typename PointInfo, int Dim = 2>
  class HierarchicalClustering
  {
    public:
      typedef typename boost::array<DoubleReal, Dim> Coord;

    protected:
      template <typename T>
      class CoordHash
      {
        public:
          std::size_t operator() (const T &p) const
          {
            boost::hash<Coord> hasher;
            return hasher(p.coord);
          }
      };

    public:
      class Point
      {
        public:
          const Coord coord;
          const PointInfo info;

          Point(const Coord &coord, const PointInfo &info)
            : coord(coord), info(info)
          { }

          // XXX
          bool operator==(const Point &) const { return false; }
      };

      class Subset
        : public boost::unordered_multiset<Point, CoordHash<Point> >
      {
        public:
          typedef typename boost::unordered_multiset<Point, CoordHash<Point> > Points;

          const Coord coord;

          Subset(const Coord &coord)
            : coord(coord)
          { }

          Subset(const Coord &coord, const Points &p)
            : Points(p), coord(coord)
          { }

          // XXX
          bool operator==(const Subset &d) { return false; }
      };

      typedef HashGrid<Subset, Dim, CoordHash<Subset> > Grid;

      const Coord max_delta;

      Grid grid;

      HierarchicalClustering(const Coord &max_delta)
        : max_delta(max_delta)
      { }

      /**
       * Insert new Point into correct grid cell.
       */
      void insertPoint(const Coord &d, const PointInfo &info)
      {
        typename Subset::Points p;
        p.insert(Point(d, info));
        newSubset(d, p);
      }

    protected:
      /**
       * Insert new Subset into correct grid cell.
       * Returns iterator to inserted item.
       */
      typename Grid::mapped_type::iterator newSubset(const Coord &d, const typename Subset::Points &p)
      {
        typename Grid::Point cellcoord;
        for (UInt i = 0; i < Dim; ++i) cellcoord[i] = d[i] / max_delta[i];
        typename Grid::mapped_type &cell = grid[cellcoord];
        return cell.insert(Subset(d, p));
      }
  };
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H */
