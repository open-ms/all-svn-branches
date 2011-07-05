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

typedef OpenMS::HashGrid<Value, 2> Test2D;
const Test2D::Point max_delta = {{1, 1}};

START_TEST(HashGrid, "$Id$")

START_SECTION(HashGrid::max_key)
{
  Test2D t(max_delta);
  TEST_EQUAL(t.max_key[0], 0);
  TEST_EQUAL(t.max_key[1], 0);
}
END_SECTION

START_SECTION(local_iterator insert(const value_type& obj))
{
  Test2D t(max_delta);
  const Test2D::Point key1 = {{1, 2}};
  Test2D::local_iterator it = t.insert(std::make_pair(key1, Test2D::mapped_type()));
  TEST_EQUAL(t.max_key[0], key1[0]);
  TEST_EQUAL(t.max_key[1], key1[1]);
  TEST_EQUAL(it->first[0], key1[0]);
  TEST_EQUAL(it->first[1], key1[1]);
  const Test2D::Point key2 = {{2, 3}};
  it = t.insert(std::make_pair(key2, Test2D::mapped_type()));
  TEST_EQUAL(t.max_key[0], key2[0]);
  TEST_EQUAL(t.max_key[1], key2[1]);
  TEST_EQUAL(it->first[0], key2[0]);
  TEST_EQUAL(it->first[1], key2[1]);
}
END_SECTION

START_SECTION(size_type erase(const key_type& key))
{
  Test2D t(max_delta);
  const Test2D::Point key = {{1, 2}};
  t.insert(std::make_pair(key, Test2D::mapped_type()));
  TEST_EQUAL(t.erase(key), 1);
}
END_SECTION

END_TEST

