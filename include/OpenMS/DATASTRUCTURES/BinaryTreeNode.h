// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef BINARYTREENODE_H_
#define BINARYTREENODE_H_

#include <OpenMS/DATASTRUCTURES/DataPoint.h>


namespace OpenMS
{

class BinaryTreeNode {
public:
	DataPoint* data1;
	DataPoint* data2;
	DoubleReal distance;
	BinaryTreeNode();
	BinaryTreeNode(DataPoint* data1_,DataPoint* data2_,DoubleReal distance_);
	bool operator==(const BinaryTreeNode &cp) const;
	bool operator<(const BinaryTreeNode &cp) const;

};
}

#endif /* BINARYTREENODE_H_ */
