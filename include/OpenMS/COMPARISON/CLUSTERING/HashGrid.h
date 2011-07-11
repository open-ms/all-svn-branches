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
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <boost/array.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <OpenMS/COMPARISON/CLUSTERING/Hasher.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HASHGRID_H
#define OPENMS_COMPARISON_CLUSTERING_HASHGRID_H

namespace OpenMS
{
  /**
   * @brief Container for (n-dimensional point, value) pairs.
   *
   * It stores all points in cells with a discrete size.
   */
  template <typename Value, std::size_t Dim>
  class OPENMS_DLLAPI HashGrid
  {
    public:
      /**
       * @brief Point for stored pairs.
       */
      // XXX: Check is there is another type handy in OpenMS allready
      typedef boost::array<DoubleReal, Dim> Point;
      /**
       * @brief Index for cells.
       */
      typedef boost::array<UInt, Dim> CellIndex;

      typedef typename boost::unordered_multimap<Point, Value> Cell;
      typedef boost::unordered_map<CellIndex, Cell> CellMap;

      typedef typename Cell::key_type key_type;
      typedef typename Cell::mapped_type mapped_type;
      typedef typename Cell::value_type value_type;

      typedef typename CellMap::const_iterator const_cell_iterator;
      typedef typename CellMap::iterator cell_iterator;
      typedef typename Cell::const_iterator const_local_iterator;
      typedef typename Cell::iterator local_iterator;
      typedef typename Cell::size_type size_type;

    private:
      CellMap cells_;
      CellIndex max_key_;

    public:
      /**
       * @brief Size of each cell.
       */
      const Point max_delta;
      /**
       * @brief Upper-right corner of key space for cells.
       */
      const CellIndex &max_key;

    public:
      HashGrid(const Point &max_delta)
        : max_delta(max_delta), max_key(max_key_)
      {
        // XXX: constructor?
        for (typename CellIndex::iterator it = max_key_.begin(); it != max_key_.end(); ++it) *it = 0;
      }

      /**
       * @brief Inserts a std::pair.
       * @param v Pair to be inserted.
       * @return Iterator that points to the inserted pair.
       */
      local_iterator insert(const value_type &v)
      {
        const CellIndex cellkey = key_to_cellkey(v.first);
        Cell &cell = cells_[cellkey];
        update_max_key(cellkey);
        return cell.insert(v);
      }

      /**
       * @brief Erases elements matching the key.
       * @param x Key of element to be erased.
       * @return Number of elements erased.
       */
      size_type erase(const key_type &x)
      {   
        const CellIndex cellkey = key_to_cellkey(x);
        try
        {
          Cell &cell = cells_.at(cellkey);
          return cell.erase(x);
        }
        catch (std::out_of_range &) { }
        return 0;
      }

      void clear() { cells_.clear(); }

      const_cell_iterator cell_begin() const { return cells_.begin(); }
      const_cell_iterator cell_end() const { return cells_.end(); }

      // XXX: Currently needes non-const
      typename CellMap::mapped_type &cell_at(const CellIndex &x) { return cells_.at(x); }
      const typename CellMap::mapped_type &cell_at(const CellIndex &x) const { return cells_.at(x); }

      void cell_clear(const CellIndex &x)
      {
        try
        {
          cells_.at(x).clear();
        }
        catch (std::out_of_range &) { }
      }

    private:
      // XXX
      CellIndex key_to_cellkey(const Point &key)
      {
        CellIndex ret;
        typename CellIndex::iterator it = ret.begin();
        typename Point::const_iterator lit = key.begin(), rit = max_delta.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit / *rit;
        return ret;
      }

      void update_max_key(const CellIndex &d)
      {
        typename CellIndex::const_iterator it1 = d.begin();
        typename CellIndex::iterator it2 = max_key_.begin();
        for (; it1 != d.end(); ++it1, ++it2)
        {
          if (*it1 > *it2)
            *it2 = *it1;
        }
      }
  };
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HASHGRID_H */
