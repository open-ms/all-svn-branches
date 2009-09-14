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

#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/MappingThread.h>
static int iRefTh = 0;
namespace OpenMS
{
  // publics membres
  
  MappingThread::MappingThread(LayerData* parent)
    : QThread(),
      cols_(0), rows_(0), 
      mz_min_(0.0), mz_max_(0.0), mz_width_(0.0),
      rt_min_(0.0), rt_max_(0.0), rt_width_(0.0),
      parent_(parent)
  {
cout << "mapping thread constructor" << endl;
iRef_ = ++iRefTh;
  }
  
  MappingThread::~MappingThread()
  {
cout << "mapping thread desconstructor " << iRef_ << endl;
  }

  void MappingThread::run()
  {
cout << "20 " << iRef_ << endl;
    // update data
    if(data_.empty())
      updateData_();
cout << "21" << endl;

    if(vertex_.empty())
      updateVertex_();
cout << "22" << endl;

    if(normals_.empty())
      updateNormals_();

cout << "23" << endl;
    emit(finish());
cout << "24" << data_.size() << " " << vertex_.size() << " " << normals_.size() << endl;
  }
  
  // get and set membres
    
  void MappingThread::setDataSize(const Size cols,const Size rows)
  {  
    cols_ = cols;
    rows_ = rows;
    data_.clear();
  }
  
  void MappingThread::setRange(const double& mz_min, const double& mz_max, const double& rt_min, const double& rt_max)
  { 
cout << "set range" << endl;
    mz_min_ = mz_min;
    mz_max_ = mz_max;
    rt_min_ = rt_min;
    rt_max_ = rt_max;
    data_.clear();
  }
      
  Size MappingThread::getRowSize()
  {
    switch(parent_->getMappingMode())
    {
      case MappingThread::MM_NONE :
        return 0;
      case MappingThread::MM_POINTS :
      case MappingThread::MM_PEAKS :	    
        return rows_;
      case MappingThread::MM_MAP :
        return rows_ - 1;
      case MappingThread::MM_PSEUDOGEL :
        return (rows_ * 2) - 1;
      default : 
        return 0;	        
    }
  }
  
  Size MappingThread::getColSize()
  {
    switch(parent_->getMappingMode())
    {
      case MappingThread::MM_NONE :
        return 0;
      case MappingThread::MM_POINTS :
        return rows_;
      case MappingThread::MM_PEAKS :
      case MappingThread::MM_MAP :
      case MappingThread::MM_PSEUDOGEL :
        return rows_ * 2;
      default : 
        return 0;		      
    }
  }
      
  const Vector3d& MappingThread::getVertex()
  {
    return vertex_;
  }
  
  const Vector3d& MappingThread::getNormals()
  {
    return normals_; 
  }
  
  void MappingThread::clearData()
  {
    data_.clear();
  }

  bool MappingThread::isValide()
  {
    bool valide = !data_.empty() && !vertex_.empty();
    switch(parent_->getMappingMode())
    {
      case MappingThread::MM_NONE :
	    case MappingThread::MM_POINTS :    
	    case MappingThread::MM_PEAKS :
	      break;
	    case MappingThread::MM_MAP :
	    case MappingThread::MM_PSEUDOGEL :
	      valide &= !normals_.empty();
	  }
	      
cout << "isValide: " << (data_.empty() ? 0 : data_.size()) << " " << (normals_.empty() ? 0 : normals_.size()) << " " << (vertex_.empty() ? 0 : vertex_.size()) << endl;
	  return valide;
  }
  
  // privates members

  void MappingThread::updateData_()
  {
cout << "update data " << iRef_ << endl;
    data_.clear();

    if(MappingThread::MM_NONE == parent_->getMappingMode())
      return;

    mz_width_ = (mz_max_ - mz_min_) / (cols_ - 1);
    rt_width_ = (rt_max_ - rt_min_) / (rows_ - 1);

    data_.resize(cols_*rows_, 0.0);
cout << "30 c:" << cols_ << " r:" << rows_ << endl;
cout << rt_min_ << " " << rt_max_ << " " << mz_min_ << " " << mz_max_ << endl;
int ii = 0;
    for(Spectrum3DCanvas::ExperimentType::ConstAreaIterator it = parent_->peaks.areaBeginConst(rt_min_, rt_max_, mz_min_, mz_max_);
        it != parent_->peaks.areaEndConst();
        ++it)
    {
++ii;
      if(it->getMZ()>=mz_min_ && it->getMZ()<=mz_max_ && it.getRT()>=rt_min_ && it.getRT()<=rt_max_)
      {		      
        if( it->getIntensity() > getData_(it->getMZ(), it.getRT()) )
        {
          setData_(it->getMZ(), it.getRT(), it->getIntensity());
        }
      }
    }
cout << "ii " << ii << endl;

		vertex_.clear();
		normals_.clear();
  }
  
  void MappingThread::updateVertex_()
  {
    switch(parent_->getMappingMode())
    {
      case MappingThread::MM_NONE :
        return;
        
	    case MappingThread::MM_POINTS :    
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
	    
	    case MappingThread::MM_PEAKS :
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
	    
	    case MappingThread::MM_MAP :
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
	    
	    case MappingThread::MM_PSEUDOGEL :
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
  
  void MappingThread::updateNormals_()
  {
    switch(parent_->getMappingMode())
    {
      case MappingThread::MM_NONE :
        return;
        	    
	    case MappingThread::MM_POINTS :    
	    case MappingThread::MM_PEAKS :
	      break;

	    case MappingThread::MM_MAP :
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
	    
	    case MappingThread::MM_PSEUDOGEL :
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
                
  Size MappingThread::getPosition_(const int col, const int row) const
  {
    int col_ = col < 0 ? 0 : col;
    col_ = col < (int) cols_ ? col : (int) cols_-1;
    
    int row_ = row < 0 ? 0 : row;
    row_ = row < (int) rows_ ? row : (int) rows_-1;
    
    return (Size) (col_ + cols_ * row_); 
  }

  Size MappingThread::getPosition_(const double mz, const double rt) const
  {
    return getPosition_(mzToIndex_(mz), rtToIndex_(rt));
  }

  Size MappingThread::getPosition_(const Size col, const Size row) const
  {
    return getPosition_((int) col, (int) row);
  }

  double MappingThread::indexToMz_(const Size col) const
  { 
    return (mz_min_ + (col<cols_ ? col : cols_-1) * mz_width_);
  }
  
  double MappingThread::indexToRt_(const Size row) const
  {
    return (rt_min_ + (row<rows_ ? row : rows_-1) * rt_width_);
  }
  
  Size MappingThread::mzToIndex_(const double mz) const
  {
    if(mz < mz_min_) return 0;
    if(mz > mz_max_) return (cols_-1);

    return (Size) floor(0.5 + (mz-mz_min_) / mz_width_);
  }
  
  Size MappingThread::rtToIndex_(const double rt) const
  {
    if(rt < rt_min_) return 0;
    if(rt > rt_max_) return (rows_-1);
       
    return (Size) floor(0.5 + (rt-rt_min_) / rt_width_);
  }

  void MappingThread::setData_(const Size col, const Size row, const double intensity)
  {
    if(col<cols_ && row<rows_)
    {
      data_[getPosition_(col, row)] = intensity;
    }
  }
  
  double MappingThread::getData_(const Size col, const Size row) const
  {
    return data_[getPosition_(col, row)];
  }

  void MappingThread::setData_(const double mz, const double rt, const double intensity)
  {
    if(mz>=mz_min_ && mz<=mz_max_ && rt>=rt_min_ && rt<=rt_max_)
    {
      data_[getPosition_(mz, rt)] = intensity;
    }
  }
  
  double MappingThread::getData_(const double mz, const double rt) const
  {
    return data_[getPosition_(mz, rt)];
  }
  
} // namespace OpenMS
