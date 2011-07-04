// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Lars Nilse, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>

using namespace OpenMS;

class Value
{
};
std::size_t hash_value(const Value&)
{ return 0; }

typedef OpenMS::HashGrid<Value, 2> Test2D;

START_TEST(HashGrid, "$Id$")

START_SECTION(HashGrid::max_key)
{
  Test2D t;
  TEST_EQUAL(t.max_key()[0], 0);
  TEST_EQUAL(t.max_key()[1], 0);
}
END_SECTION

START_SECTION(mapped_type& operator[](const key_type& key))
{
  Test2D t;
  const Test2D::Point key = {{1, 2}};
  Test2D::mapped_type a = t[key];
  TEST_EQUAL(t.max_key()[0], key[0]);
  TEST_EQUAL(t.max_key()[1], key[1]);
}
END_SECTION

START_SECTION(iterator insert(const value_type& obj))
{
  Test2D t;
  const Test2D::Point key = {{1, 2}};
  t.insert(std::make_pair(key, Test2D::mapped_type()));
  TEST_EQUAL(t.max_key()[0], key[0]);
  TEST_EQUAL(t.max_key()[1], key[1]);
}
END_SECTION

END_TEST

