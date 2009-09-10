// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2009
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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/GridData.h>
#include <OpenMS/VISUAL/MultiGradient.h>

namespace OpenMS
{
  // publics membres
  
  MapData::MapData()
    : QThread(),
      valideData_(false), valideVertex_(false), valideNormals_(false), valideColors_(false),
      needVertex_(false), needNormals_(false), needColors_(false),
      mode_(LayerData::DM_POINTS),       
      cols_(0), rows_(0), 
      mz_min_(0.0), mz_max_(0.0), mz_width_(0.0),
      rt_min_(0.0), rt_max_(0.0), rt_width_(0.0),
      experiment_(NULL), gradient_(NULL)
  { 
  }
  
  MapData::~MapData()
  {
    clearAll_();
  }

  void MapData::run()
  { 
    // update data
    if(!valideData_)
      updateData_();
      
    if(!valideVertex_ && needVertex_)
      updateVertex_();
    
    if(!valideNormals_ && needNormals_)
      updateNormals_();
      
    if(!valideColors_ && needColors_)
      updateColors_();
      
    emit(finish());
  }
  
  // get and set membres
    
  void MapData::setDrawMode(const LayerData::DrawModes mode)
  {
    mode_ = mode;
    valideVertex_ = false;
    valideNormals_ = false;
    valideColors_ = false;
  }

  void MapData::setDataSize(const Size cols,const Size rows)
  {  
    cols_ = cols;
    rows_ = rows;
    invalidate();
  }
  
  void MapData::setRange(const double& mz_min, const double& mz_max, const double& rt_min, const double& rt_max)
  { 
    mz_min_ = mz_min;
    mz_max_ = mz_max;
    rt_min_ = rt_min;
    rt_max_ = rt_max;
    invalidate();
  }
  
  void MapData::setGradient(const MultiGradient* gradient)
  {    
    gradient_ = gradient;
  }

  void MapData::setData(const ExperimentType* experiment)
  {
    experiment_ = experiment;
    invalidate();
  }
      
  Size MapData::getRowSize()
  {
    if(valideVertex_ || valideNormals_ || valideColors_)
    {
      switch(mode_)
      {
	      case LayerData::DM_POINTS :
	      case LayerData::DM_PEAKS :	    
	        return rows_;
	      case LayerData::DM_LINES :
	        return rows_ - 1;
	      case LayerData::DM_MAP :
          return (rows_ * 2) - 1;
	      default : 
	        return 0;	        
      }
    }
    else
      return 0;
  }
  
  Size MapData::getColSize()
  {
    if(valideVertex_ || valideNormals_ || valideColors_)
    {  
      switch(mode_)
      {
	      case LayerData::DM_POINTS :
	        return rows_;
	      case LayerData::DM_PEAKS :
	      case LayerData::DM_LINES :
	      case LayerData::DM_MAP :
	        return rows_ * 2;
	      default : 
	        return 0;		      
      }
    }
    else
      return 0;
  }
      
  const Vector3d& MapData::getVertex()
  {
    return vertex_;
  }
  
  const Vector3d& MapData::getNormals()
  {
    return normals_; 
  }
  
  const VectorRGBA& MapData::getColors()
  {
    return colors_;
  }
  
  void MapData::invalidate(const bool data, const bool vertex, const bool normals, const bool colors)
  {
    valideData_ = !data;
    valideVertex_ = !vertex;
    valideNormals_ = !normals;
    valideColors_ = !colors;
  }
  
  bool MapData::isValideVertex()
  {
    return valideVertex_;
  }
  
  bool MapData::isValideNormals()
  {
    return valideNormals_;
  }
  
  bool MapData::isValideColors()
  {
    return valideColors_;
  }  

  void MapData::needVertex()
  { 
    needVertex_ = true;
  }
  
  void MapData::needNormals()
  {   
    needNormals_ = true;
  }
  
  void MapData::needColors()
  {  
    needColors_ = true;
  }
  
  // privates members

  void MapData::clearAll_()
  {
    if(!data_.empty())
      data_.clear();
      
    if(!vertex_.empty())
      vertex_.clear();
      
    if(!normals_.empty())
      normals_.clear();
      
    if(!colors_.empty())
      colors_.clear();    
      
    valideData_ = false;   
    valideVertex_ = false;
    valideNormals_ = false;
    valideColors_ = false;                   
  }

  void MapData::updateData_()
  {
    clearAll_();
       
    mz_width_ = (mz_max_ - mz_min_) / (cols_ - 1);
    rt_width_ = (rt_max_ - rt_min_) / (rows_ - 1);

    data_.resize(cols_*rows_, 0.0);

    for(AreaIt it=experiment_->areaBeginConst(rt_min_, rt_max_, mz_min_, mz_max_);
        it!=experiment_->areaEndConst();
        ++it)
    {
      if(it->getMZ()>=mz_min_ && it->getMZ()<=mz_max_ && it.getRT()>=rt_min_ && it.getRT()<=rt_max_)
      {		      
        if(it->getIntensity() > getData_(it->getMZ(), it.getRT()))
        {
          setData_(it->getMZ(), it.getRT(), it->getIntensity());
        }
      }
    }
    valideData_ = true;
  }
  
  void MapData::updateVertex_()
  {
    if(!data_.empty() && vertex_.empty())
    {
      // update vertex
	    switch(mode_)
	    {
		    case LayerData::DM_POINTS :    
		    {
		      for(Size iRow=0; iRow<rows_; ++iRow)
		      {
		        double rt = indexToRt_(iRow);
		        for(Size iCol=0; iCol<cols_; ++iCol)
		        {
		          double mz = indexToMz_(iCol);
		          double intensity = getData_(iCol, iRow);
		          vertex_.push_back(Struct3d(mz, rt, intensity));
		        }
		      }
		      break;
		    }
		    
		    case LayerData::DM_PEAKS :
		    {
		      for(Size iRow=0; iRow<rows_; ++iRow)
		      {
		        double rt = indexToRt_(iRow);
		        for(Size iCol=0; iCol<cols_; ++iCol)
		        {
		          double mz = indexToMz_(iCol);
		          double intensity = getData_(iCol, iRow);
		          vertex_.push_back(Struct3d(mz, rt, 0.0));
		          vertex_.push_back(Struct3d(mz, rt, intensity));
		        }
		      }
		      break;
		    }
		    
		    case LayerData::DM_LINES :
		    {
		      for(Size iRow=0; iRow<(rows_-1); ++iRow)
		      {
		        double rt1 = indexToRt_(iRow);
		        double rt2 = indexToRt_(iRow+1);
		        for(Size iCol=0; iCol<cols_; ++iCol)
		        {
		          double mz = indexToMz_(iCol);
		          double intensity1 = getData_(iCol, iRow);
		          double intensity2 = getData_(iCol, iRow+1);
		          vertex_.push_back(Struct3d(mz, rt1, intensity1));
		          vertex_.push_back(Struct3d(mz, rt2, intensity2));
		        }
		      }
		      break;
		    }
		    
		    case LayerData::DM_MAP :
		    {
		      for(Size iRow=0; iRow<rows_; ++iRow)
		      {
		        double rt = indexToRt_(iRow);

		        for(Size iCol=0; iCol<cols_; ++iCol)
		        {
		          double mz = indexToMz_(iCol);
		          double intensity = getData_(iCol, iRow);
		          vertex_.push_back(Struct3d(mz, rt - (rt_width_ / 2.0), intensity));
		          vertex_.push_back(Struct3d(mz, rt + (rt_width_ / 2.0), intensity));
		        }
		        
		        if(iRow < (rows_ - 1))
		        {
		          for(Size iCol=0; iCol<cols_; ++iCol)
		          {
		            double mz = indexToMz_(iCol);
		            double intensity1 = getData_(iCol, iRow);
		            double intensity2 = getData_(iCol, iRow+1);
		            vertex_.push_back(Struct3d(mz, rt + (rt_width_ / 2.0), intensity1));
		            vertex_.push_back(Struct3d(mz, rt + (rt_width_ / 2.0), intensity2));
		          }	     
		        }   
		      }
		      break;
		    }
		    
		    default : 
		      break;
		  }
    }
    valideVertex_ = true; 
  }
  
  void MapData::updateNormals_()
  {
    if(!data_.empty() && normals_.empty())
    {
      // update vertex
	    switch(mode_)
	    {
		    case LayerData::DM_POINTS :    
		    case LayerData::DM_PEAKS :
		      break;

		    case LayerData::DM_LINES :
		    {
		      for(Size iRow=0; iRow<(rows_-1); ++iRow)
		      {
		        for(Size iCol=0; iCol<cols_; ++iCol)
		        {
		          normals_.push_back(Struct3d(0.0, 0.0, 1.0));
		          normals_.push_back(Struct3d(0.0, 0.0, 1.0));
		        }
		      }
		      break;
		    }
		    
		    case LayerData::DM_MAP :
		    {
		      for(Size iRow=0; iRow<rows_; ++iRow)
		      {
		        for(Size iCol=0; iCol<cols_; ++iCol)
		        {
		          normals_.push_back(Struct3d(0.0, 0.0, 1.0));
		          normals_.push_back(Struct3d(0.0, 0.0, 1.0));
		        }
		        
		        if(iRow < (rows_ - 1))
		        {
		          for(Size iCol=0; iCol<cols_; ++iCol)
		          {
		          normals_.push_back(Struct3d(0.0, 0.0, 1.0));
		          normals_.push_back(Struct3d(0.0, 0.0, 1.0));
		          }	     
		        }   
		      }
		      break;
		    }

		    default : 
		      break;		    
		  }
    }
    valideNormals_ = true; 
  }
    
  void MapData::updateColors_()
  { 
    if(!data_.empty() && colors_.empty())
    {
      if(!valideVertex_)
        updateVertex_();
        
      for(Iterator3d it=vertex_.begin(); it!=vertex_.end(); ++it)
      {
        colors_.push_back(ColorRGBA(gradient_->interpolatedColorAt(it->intensity)));
      }
    }
    valideColors_ = true;   
  }  
                    
  Size MapData::getPosition_(const int col, const int row) const
  {
    int col_ = col < 0 ? 0 : col;
    col_ = col < (int) cols_ ? col : (int) cols_-1;
    
    int row_ = row < 0 ? 0 : row;
    row_ = row < (int) rows_ ? row : (int) rows_-1;
    
    return (Size) (col_ + cols_ * row_); 
  }

  Size MapData::getPosition_(const double mz, const double rt) const
  {
    return getPosition_(mzToIndex_(mz), rtToIndex_(rt));
  }

  Size MapData::getPosition_(const Size col, const Size row) const
  {
    return getPosition_((int) col, (int) row);
  }

  double MapData::indexToMz_(const Size col) const
  { 
    return (mz_min_ + (col<cols_ ? col : cols_-1) * mz_width_);
  }
  
  double MapData::indexToRt_(const Size row) const
  {
    return (rt_min_ + (row<rows_ ? row : rows_-1) * rt_width_);
  }
  
  Size MapData::mzToIndex_(const double mz) const
  {
    if(mz < mz_min_) return 0;
    if(mz > mz_max_) return (cols_-1);

    return (Size) floor(0.5 + (mz-mz_min_) / mz_width_);
  }
  
  Size MapData::rtToIndex_(const double rt) const
  {
    if(rt < rt_min_) return 0;
    if(rt > rt_max_) return (rows_-1);
       
    return (Size) floor(0.5 + (rt-rt_min_) / rt_width_);
  }

  void MapData::setData_(const Size col, const Size row, const double intensity)
  {   
    if(col<cols_ && row<rows_)
    {
      data_[getPosition_(col, row)] = intensity;
    }
  }
  
  double MapData::getData_(const Size col, const Size row) const
  { 
    return data_[getPosition_(col, row)];
  }

  void MapData::setData_(const double mz, const double rt, const double intensity)
  { 
    if(mz>=mz_min_ && mz<=mz_max_ && rt>=rt_min_ && rt<=rt_max_)
    {
      data_[getPosition_(mz, rt)] = intensity;
    }
  }
  
  double MapData::getData_(const double mz, const double rt) const
  { 
    return data_[getPosition_(mz, rt)];
  }
  
} // namespace OpenMS
