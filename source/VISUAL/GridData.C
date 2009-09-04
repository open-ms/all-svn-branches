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
#include <OpenMS/MATH/MISC/Struct3d.h>
#include <OpenMS/VISUAL/ColorRGBA.h>
#include <OpenMS/VISUAL/Arrow3d.h>

namespace OpenMS
{
  // publics membres
  
  MapData::MapData()
    : QThread(),
      valideData_(false), valideVertex_(false), valideNormals_(false), valideColors_(false),
      needVertex_(false), needNormals_(false), needColors_(false),
      mode_(MAP),       
      cols_(0), rows_(0), 
      mz_min_(0.0), mz_max_(0.0), mz_width_(0.0),
      rt_min_(0.0), rt_max_(0.0), rt_width_(0.0),
      begin_(), end_(), gradient_(NULL)
  {
cout << "MapData constructor" << endl;  
  }
  
  MapData::~MapData()
  {
cout << "MapData delete" << endl;  
    clearAll_();
  }

  void MapData::run()
  {
cout << "begin run..." << endl;  
    
    // update data
    if(!valideData_)
      updateData_();
      
    if(!valideVertex_ && needVertex_)
    {
      updateVertex_();
cout << "data: " << data_.size() << " - vertex:" << vertex_.size() << endl;
      // draw    
      Arrow3d arrow;
      for(vector<Struct3d>::iterator itVertex=vertex_.begin(); itVertex!=vertex_.end(); ++itVertex)
      {
        //arrow.draw(*itVertex, Struct3d(0.0, 1.0, 0.0), 50.0);                                    		        
      }
    }
    
    if(!valideNormals_ && needNormals_)
      updateNormals_();
      
    if(!valideColors_ && needColors_)
      updateColors_();

cout << "wait" << endl;
sleep(5);
cout << "end run..." << endl;       
  }
  
  // get and set membres
    
  void MapData::setInterpolationMode(const MapData::interpolationMode mode)
  {
    mode_ = mode;
    invalidate(false);
  }

  void MapData::setDataSize(const Size cols,const Size rows)
  {
cout << "set data size" << endl;    
    cols_ = cols;
    rows_ = rows;
    invalidate(false);
  }
  
  void MapData::setRange(const double& mz_min, const double& mz_max, const double& rt_min, const double& rt_max)
  {
cout << "set range" << endl;    
    mz_min_ = mz_min;
    mz_max_ = mz_max;
    rt_min_ = rt_min;
    rt_max_ = rt_max;
    invalidate(false);
  }
  
  void MapData::setGradient(const MultiGradient* gradient)
  {
cout << "set gradient" << endl;    
    gradient_ = gradient;
  }
  
  void MapData::setData(const AreaIt begin, const AreaIt end)
  {
cout << "set data" << endl;  
    begin_ = begin;
    end_ = end;
    invalidate();
  }
      
  Size MapData::getRowSize()
  {
    if(isValidVertex() || isValidNormals() || isValidColors())
    {
      switch(mode_)
      {
	      case LINE :
	      case BARPLOT :	    
	        return rows_;
	      case MAP :
	        return rows_ - 1;
	      case PSEUDOGEL :
          return (rows_ * 2) - 1;
	      case NONE :
	      default : 
	        return 0;	        
      }
    }
    else
      return 0;
  }
  
  Size MapData::getColSize()
  {
    if(isValidVertex() || isValidNormals() || isValidColors())
    {  
      switch(mode_)
      {
	      case LINE :
	        return rows_;
	      case BARPLOT :
	      case MAP :
	      case PSEUDOGEL :
	        return rows_ * 2;
	      case NONE :
	      default : 
	        return 0;		      
      }
    }
    else
      return 0;
  }
      
  const vector<Struct3d>& MapData::getVertex()
  {
    return vertex_;
  }
  
  const vector<Struct3d>& MapData::getNormals()
  {
    return normals_; 
  }
  
  const vector<ColorRGBA>& MapData::getColors()
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
  
  bool MapData::isValidVertex()
  {
    return valideVertex_;
  }
  
  bool MapData::isValidNormals()
  {
    return valideNormals_;
  }
  
  bool MapData::isValidColors()
  {
    return valideColors_;
  }  

  void MapData::needVertex()
  {
cout << "need vertex..." << endl;    
    needVertex_ = true;
  }
  
  void MapData::needNormals()
  {
cout << "need normals..." << endl;    
    needNormals_ = true;
  }
  
  void MapData::needColors()
  {
cout << "need colors..." << endl;  
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
  }

  void MapData::updateData_()
  {  
cout << "update data...";    
    valideVertex_ = false;
    valideNormals_ = false;
    valideColors_ = false;
    clearAll_();
       
    mz_width_ = (mz_max_ - mz_min_) / (cols_ - 1);
    rt_width_ = (rt_max_ - rt_min_) / (rows_ - 1);

    data_.resize(cols_*rows_, 0.0);
              
    for(AreaIt it=begin_; it!=end_; ++it)
    {
      if(it->getMZ()>=mz_min_ && it->getMZ()<=mz_max_ && it.getRT()>=rt_min_ && it.getRT()<=rt_max_)
      {		      
        if(it->getIntensity() > getData_(it->getMZ(), it.getRT()))
        {
          setData_(it->getIntensity(), it->getMZ(), it.getRT());
        }
      }
    }
    valideData_ = true;
cout << " done" << endl;     
  }
  
  void MapData::updateVertex_()
  {
cout << "update vertex...";  
    if(!data_.empty() && vertex_.empty())
    {
      // update vertex
	    switch(mode_)
	    {
		    case LINE :    
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
		    
		    case BARPLOT :
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
		    
		    case MAP :
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
		    
		    case PSEUDOGEL :
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
		    
		    case NONE :
		      break;
		    default : 
		      break;
		  }
    }
    valideVertex_ = true;
cout << " done" << endl;     
  }
  
  void MapData::updateNormals_()
  {
cout << "update normals...";  
    if(!data_.empty() && normals_.empty())
    {
      // update vertex
	    switch(mode_)
	    {
		    case LINE :    
		    case BARPLOT :
		      break;

		    case MAP :
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
		    
		    case PSEUDOGEL :
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

		    case NONE :
		      break;
		    default : 
		      break;		    
		  }
    }
    valideNormals_ = true;
cout << " done" << endl;     
  }
    
  void MapData::updateColors_()
  {
cout << "update colors...";  
    if(!data_.empty() && colors_.empty())
    {
      if(!valideVertex_)
        updateVertex_();
        
      for(vector<Struct3d>::iterator it=vertex_.begin(); it!=vertex_.end(); ++it)
      {
        colors_.push_back(ColorRGBA(gradient_->interpolatedColorAt(it->z)));
      }
    }
    valideColors_ = true;
cout << " done" << endl;    
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

  void MapData::setData_(const double value, const Size col, const Size row)
  {   
    if(col<cols_ && row<rows_)
    {
      data_[getPosition_(col, row)] = value;
    }
  }
  
  double MapData::getData_(const Size col, const Size row) const
  { 
    return data_[getPosition_(col, row)];
  }

  void MapData::setData_(const double value, const double mz, const double rt)
  { 
    if(mz>=mz_min_ && mz<=mz_max_ && rt>=rt_min_ && rt<=rt_max_)
    {
      data_[getPosition_(mz, rt)] = value;
    }
  }
  
  double MapData::getData_(const double mz, const double rt) const
  { 
    return data_[getPosition_(mz, rt)];
  }
  
} // namespace OpenMS
