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


#include <boost/array.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <OpenMS/COMPARISON/CLUSTERING/Hasher.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HASHGRID_H
#define OPENMS_COMPARISON_CLUSTERING_HASHGRID_H

namespace OpenMS
{
  template <typename Value, std::size_t Dim>
  class OPENMS_DLLAPI HashGrid
  {
    public:
      typedef boost::array<DoubleReal, Dim> Point;
      typedef boost::array<UInt, Dim> CellPoint;

      typedef typename boost::unordered_multimap<Point, Value> Cell;
      typedef boost::unordered_map<CellPoint, Cell> CellMap;

      typedef typename Cell::key_type key_type;
      typedef typename Cell::mapped_type mapped_type;
      typedef typename Cell::value_type value_type;

      typedef typename Cell::const_iterator const_local_iterator;
      typedef typename Cell::iterator local_iterator;
      typedef typename Cell::size_type size_type;

    protected:
      CellMap cells_;
      CellPoint max_key_;

    public:
      const CellMap &cells;
      const Point max_delta;
      const CellPoint &max_key;

    public:
      HashGrid(const Point &max_delta)
        : cells(cells_), max_delta(max_delta), max_key(max_key_)
      {
        for (typename CellPoint::iterator it = max_key_.begin(); it != max_key_.end(); ++it) *it = 0;
      }

      local_iterator insert(const value_type& obj)
      {
        const CellPoint cellkey = key_to_cellkey(obj.first);
        Cell &cell = cells_[cellkey];
        update_max_key(cellkey);
        return cell.insert(obj);
      }

      size_type erase(const key_type& key)
      {   
        try
        {
          Cell &cell = cells_.at(key_to_cellkey(key));
          return cell.erase(key);
        }
        catch (std::out_of_range &)
        { }
        return 0;
      }

      void clear()
      {
          cells.clear();
      }

    protected:
      CellPoint key_to_cellkey(const Point &key)
      {
        CellPoint ret;
        typename CellPoint::iterator it = ret.begin();
        typename Point::const_iterator lit = key.begin(), rit = max_delta.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit / *rit;
        return ret;
      }

      void update_max_key(const CellPoint &d)
      {
        typename CellPoint::const_iterator it1 = d.begin();
        typename CellPoint::iterator it2 = max_key_.begin();
        for (; it1 != d.end(); ++it1, ++it2)
        {
          if (*it1 > *it2)
            *it2 = *it1;
        }
      }
  };
}

#endif /* HASHGRID_H_ */
