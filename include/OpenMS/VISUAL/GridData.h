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

#ifndef OPENMS_VISUAL_GRIDDATA_H
#define OPENMS_VISUAL_GRIDDATA_H

#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/MATH/MISC/Struct3d.h>
#include <OpenMS/VISUAL/ColorRGBA.h>
#include <QThread>

using namespace std;

namespace OpenMS
{
  class MultiGradient;
  
  class MapData
    : public QThread  
  {
    Q_OBJECT
      
    public:
      typedef SpectrumCanvas::ExperimentType ExperimentType;
      typedef ExperimentType::ConstAreaIterator AreaIt;
        
    public:
      MapData();
      ~MapData();
      void run();
        
      void setDrawMode(const LayerData::DrawModes mode);
      void setDataSize(const Size cols, const Size rows);
      void setRange(const double& mz_min, const double& mz_max, const double& rt_min, const double& rt_max);
      void setData(const ExperimentType* experiment);
      void setGradient(const MultiGradient* gradient);

      void needVertex();
      void needNormals();
      void needColors();
                 
      Size getRowSize();
      Size getColSize();
            
      const Vector3d& getVertex();
      const Vector3d& getNormals();
      const vector<ColorRGBA>& getColors();
      
      void invalidate(
        const bool data = true, 
        const bool vertex = true, 
        const bool normals = true,
        const bool colors = true);

      bool isValideVertex();
      bool isValideNormals();
      bool isValideColors();        

    signals:
      void finish();
              
    private:              
      void clearAll_();

      void updateData_();      
      void updateVertex_();
      void updateNormals_();
      void updateColors_();
      
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
      vector<double> data_; // set with AreaIT ; no get
      Vector3d vertex_; // no set ; get
      Vector3d normals_; // no set ; get
      vector<ColorRGBA> colors_; // no set ; get
      
      bool valideData_;
      bool valideVertex_; // no set ; get
      bool valideNormals_; // no set ; get
      bool valideColors_; // no set ; get

      bool needVertex_; // no set ; get
      bool needNormals_; // no set ; get
      bool needColors_; // no set ; get

      LayerData::DrawModes mode_; // set and get            
      Size cols_;
      Size rows_;
      double mz_min_;
      double mz_max_;
      double mz_width_;
      double rt_min_;
      double rt_max_;
      double rt_width_;
      const ExperimentType* experiment_;
      const MultiGradient* gradient_;
  };

} // namespace OpenMS

#endif // OPENMS_VISUAL_GRIDDATA_H
