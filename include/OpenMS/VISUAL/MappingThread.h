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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MAPPINGTHREAD_H
#define OPENMS_VISUAL_MAPPINGTHREAD_H

#include <OpenMS/MATH/MISC/Struct3d.h>
#include <QThread>

namespace OpenMS
{ 
	class LayerData;

  class OPENMS_DLLAPI MappingThread
    : public QThread  
  {
    Q_OBJECT
      
    public:

			///Enumerate all avaiable mapping modes
			enum MappingModes 
			{
			  MM_NONE,
				MM_POINTS,
				MM_PEAKS,
				MM_MAP,
				MM_PSEUDOGEL
			};

      MappingThread(LayerData* parent);
      void run();
      
      void setDataSize(const Size cols, const Size rows);
      void setRange(const double& mz_min, const double& mz_max, const double& rt_min, const double& rt_max);

      Size getRowSize();
      Size getColSize();
      double getRtMin();
      double getRtMax();            
      
      const Vector3d& getVertex();
      const Vector3d& getNormals();
      
      void clearData();
      bool isValide();
        
    signals:
      void finish();
              
    private:              
      void updateData_();      
      void updateVertex_();
      void updateNormals_();
      
      Size getPosition_(const int col, const int row) const;
      Size getPosition_(const Size col, const Size row) const;      
      Size getPosition_(const double mz, const double rt) const;

      double indexToMz_(const Size col) const;
      double indexToRt_(const Size row) const;
      Size mzToIndex_(const double mz) const;
      Size rtToIndex_(const double rt) const;

      void setData_(const Size col, const Size row, const double intensity);
      double getData_(const Size col, const Size row) const;
      void setData_(const double mz, const double rt, const double intensity);
      double getData_(const double mz, const double rt) const;
      
    private:                    
      vector<double> data_;
      Vector3d vertex_;
      Vector3d normals_;
   
      Size cols_;
      Size rows_;
      double mz_min_;
      double mz_max_;
      double mz_width_;
      double rt_min_;
      double rt_max_;
      double rt_width_;
      LayerData* parent_;
  };

} // namespace OpenMS

#endif // OPENMS_VISUAL_MAPPINGTHREAD_H
