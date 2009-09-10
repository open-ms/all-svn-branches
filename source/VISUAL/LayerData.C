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
		  label(L_NONE),
		  map_(NULL),
		  draw_mode_(DM_POINTS),
		  primitive_mode_(PM_POINTS)
  {
  }

  LayerData::LayerData(const LayerData& layer)
	  : visible(layer.visible),
		  flipped(layer.flipped),
		  type(layer.type),
		  name(layer.name),
		  filename(layer.filename),
		  peaks(layer.peaks),
		  features(layer.features),
		  consensus(layer.consensus),
		  current_spectrum(layer.current_spectrum),
		  f1(layer.f1),
		  f2(layer.f2),
		  f3(layer.f3),
		  param(layer.param),
		  gradient(layer.gradient),
		  filters(layer.filters),
		  annotations_1d(layer.annotations_1d),
		  modifiable(layer.modifiable),
		  modified(layer.modified),
		  label(layer.label),
		  map_(NULL),
		  draw_mode_(layer.draw_mode_),
		  primitive_mode_(layer.primitive_mode_)
  {
  }  
		  
  LayerData::~LayerData()
	{
	  if(NULL != map_)
	    delete map_;
	}
	
	LayerData& LayerData::operator= (const LayerData& layer)
	{
	  if (&layer == this) return *this;
	  
	  visible = layer.visible;
	  flipped = layer.flipped;
	  type = layer.type;
	  name = layer.name;
	  filename = layer.filename;
	  peaks = layer.peaks;
	  features = layer.features;
	  consensus = layer.consensus;
	  current_spectrum = layer.current_spectrum;
	  f1 = layer.f1;
	  f2 = layer.f2;
	  f3 = layer.f3;
	  param = layer.param;
	  gradient = layer.gradient;
	  filters = layer.filters;
	  annotations_1d = layer.annotations_1d;
	  modifiable = layer.modifiable;
	  modified = layer.modified;
	  label = layer.label;
	  map_ = NULL;
	  draw_mode_ = layer.draw_mode_;
	  primitive_mode_ = layer.primitive_mode_;
	  
	  return *this;	
	}
	
	MapData* LayerData::getMapData()
	{
	  if(NULL == map_)
	    map_ = new MapData();
	  map_->setDrawMode(draw_mode_);
	  return map_;
	}

  void LayerData::setDrawMode(const LayerData::DrawModes mode)
  {
    draw_mode_ = mode;
    if(NULL != map_)
      map_->setDrawMode(draw_mode_);
  }
  
  void LayerData::setPrimitiveMode(const LayerData::PrimitiveModes mode)
  {
    primitive_mode_ = mode;
  }	
	
	LayerData::DrawModes LayerData::getDrawMode() const
  {
    return draw_mode_;
  }
  
  LayerData::PrimitiveModes LayerData::getPrimitiveMode() const
  {
    return primitive_mode_;
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

