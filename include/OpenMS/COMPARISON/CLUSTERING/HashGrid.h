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

#include <OpenMS/COMPARISON/CLUSTERING/HashGridCell.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HASHGRID_H
#define OPENMS_COMPARISON_CLUSTERING_HASHGRID_H

// XXX: Move somewhere else
namespace boost
{
  template <class T, std::size_t N>
  std::size_t hash_value(const array<T, N>& b)
  {
    boost::hash<T> hasher;
    // XXX: Hash all items
    return hasher(b[0]);
  }
}

namespace OpenMS
{
  template <class Value, int Dim = 2>
  class OPENMS_DLLAPI HashGrid
    : public boost::unordered_map<boost::array<UInt, Dim>, OpenMS::HashGridCell<Value> >
  {
    public:
      typedef boost::array<UInt, Dim> point;
      typedef boost::unordered_map<boost::array<UInt, Dim>, OpenMS::HashGridCell<Value> > base_;

      typedef typename base_::key_type key_type;
      typedef typename base_::mapped_type mapped_type;
      typedef typename base_::value_type value_type;

      typedef typename base_::const_iterator const_iterator;
      typedef typename base_::iterator iterator;
      typedef typename base_::size_type size_type;

    protected:
      const point max_key_;

    public:
     HashGrid(const point &max_key)
       : max_key_(max_key)
     { }

     const point &max_key() const
     {
       return max_key_;
     }

     mapped_type& operator[](const key_type& key)
     {
       // XXX: Check key
       return base_::operator[](key);
     }

     iterator insert(const value_type& obj)
     {
       // XXX: Check key
       return base_::insert(obj);
     }
  };
}

#endif /* HASHGRID_H_ */
