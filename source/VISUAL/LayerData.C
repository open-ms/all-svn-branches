// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// $Authors: $
// -------------------------------------------------------------------------- 

#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/GridData.h>

using namespace std;

namespace OpenMS
{
  LayerData::LayerData()
	  : visible(true),
		  flipped(false),
		  type(DT_UNKNOWN),
		  name(),
		  filename(),
		  peaks(),
		  features(),
		  consensus(),
		  current_spectrum(0),
		  f1(false),
		  f2(false),
		  f3(false),
		  param(),
		  gradient(),
		  filters(),
		  annotations_1d(),
		  modifiable(false),
		  modified(false),
		  label(L_NONE)
  {
cout << "LayerData constructor" << endl;  
	  annotations_1d.resize(1);

    map = new MapData();  
cout << "B" << endl;
  }
		  
  LayerData::~LayerData()
	{
cout << "LayerData delete" << endl;  
	  delete map;
	}
		
	const std::string LayerData::NamesOfLabelType[] = {"None","Index","Label meta data","Peptide identification"};
	
	std::ostream& operator << (std::ostream& os, const LayerData& rhs)
	{
		os << "--LayerData BEGIN--"<<std::endl;
		os << "name: " << rhs.name << std::endl;
		os << "visible: " << rhs.visible << std::endl;
		os << "number of peaks: " << rhs.peaks.getSize() << std::endl;
		os << "--LayerData END--"<<std::endl;
		return os;
	}
}//Namespace

