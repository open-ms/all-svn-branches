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


#ifndef OPENMS_COMPARISON_CLUSTERING_HASHER_H
#define OPENMS_COMPARISON_CLUSTERING_HASHER_H

#include <boost/array.hpp>

namespace boost
{
  /** Hash Value for boost::array. */
  template <class T, std::size_t N>
  std::size_t hash_value(const array<T, N>& b)
  {
    boost::hash<T> hasher;
    std::size_t hash = 0;
    for (typename array<T, N>::const_iterator it = b.begin(); it != b.end(); ++it) hash ^= hasher(*it);
    return hash;
  }
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HASHER_H */
