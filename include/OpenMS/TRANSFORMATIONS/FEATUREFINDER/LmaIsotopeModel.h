// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{
  /** 
    @brief Isotope distribution approximated using linear interpolation.

    @htmlinclude OpenMS_LmaIsotopeModel.parameters
   */
    class OPENMS_DLLAPI LmaIsotopeModel
    : public InterpolationModel
    {

        public:
             
            typedef InterpolationModel::CoordinateType CoordinateType;
            typedef InterpolationModel::CoordinateType IntensityType;

            enum Averagines{C=0,H,N,O,S,AVERAGINE_NUM};

            /// Default constructor
            LmaIsotopeModel();

            ///  copy constructor
            LmaIsotopeModel(const LmaIsotopeModel& source);

            /// destructor
            virtual ~LmaIsotopeModel();

            /// assignment operator
            virtual LmaIsotopeModel& operator = (const LmaIsotopeModel& source);

            /// returns the charge
            UInt getCharge();

            /// create new LmaIsotopeModel object (needed by Factory)
            static BaseModel<1>* create()
            {
              return new LmaIsotopeModel();
            }

            /// name of the model (needed by Factory)
            static const String getProductName()
            {
              return "LmaIsotopeModel";
            }

            /** @brief set the offset of the model

                The whole model will be shifted to the new offset without being computing all over.
                This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
                standard deviations) but can get significant otherwise. In that case use setParameters()
                which enforces a recomputation of the model.
            */
            void setOffset(CoordinateType offset);

            /// returns the offset
            CoordinateType getOffset();

            /// set sample/supporting points of interpolation
            void setSamples();
            
            /** @brief get the monoisotopic mass of the Isotope model
                This is a m/z-value not necessarily the monoisotopic mass.
            */
            CoordinateType getCenter() const;

        protected:
          
            /// total intensity (area under curve)
            CoordinateType total_intensity_;	
            CoordinateType min_;
            CoordinateType max_;	

            CoordinateType isotope_stdev_;
            UInt charge_;
            CoordinateType mean_;
            CoordinateType monoisotopic_mz_;
            DoubleReal averagine_[AVERAGINE_NUM];
            Int max_isotope_;
            DoubleReal trim_right_cutoff_;
            DoubleReal isotope_distance_;

            void updateMembers_();
         };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEMODEL_H
