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

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACFILTER_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACFILTER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>



namespace OpenMS
{
/**
 * @brief Filter to use for SILACFiltering
 *
 * A SILACFilter searches for SILAC features, which correspond to the defined mass shifts and charge. A SILACFilter either can search
 * for doublets or triplets. Only peaks are taken into account, which were not blacklisted by other filters before e.g. are not part of
 * a feature yet.
 *
 * @see SILACFiltering
 */
class SILACFiltering;

class OPENMS_DLLAPI SILACFilter {
	friend class SILACFiltering;
private:
	/**
	 * @brief SILAC type of the filter. Either DOUBLE (2) or TRIPLE (3)
	 */
    Int silac_type;
    /**
     * @brief charge of the ions to search for
     */
    Int charge;
    /**
     * @brief distance between light and medium peaks. Used for triplet filtering
     */
    DoubleReal envelope_distance_light_medium;
    /**
     * @brief distance between light and heavy peaks
     */
    DoubleReal envelope_distance_light_heavy;
    /**
     * @brief distance between isotope peaks
     */
    DoubleReal isotope_distance;

    DataPoint next_element;

    /**
     * @brief holds the recognized features
     */
    std::vector<DataPoint> elements;

    /**
     * @brief maximal value of which a predicted SILAC feature may deviate from the averagine model
     */
    DoubleReal model_deviation;

    /**
     * @brief helper structure to compare double values
     */
    struct doubleCmp
    {
        bool operator ()(DoubleReal a, DoubleReal b) const;
    };
    /**
     * @brief holds the m/z values, which are ignored due to blacklisting of other filters
     */
    std::set<DoubleReal,doubleCmp> blacklist;
    /**
     * @brief holds the m/z values of the last identified feature. These values will be blacklisted in other filters
     */
    std::vector<DoubleReal> peak_values;

    /**
     * @brief returns the intensities of all peaks used in the identification of one SILAC feature
     * @param act_mz m/z value of the first peak
     * @param intensity interpolation
     * @param spline function of the interpolation
     */
    std::vector<DoubleReal> getIntensities(DoubleReal act_mz,gsl_interp_accel *acc, gsl_spline *spline);


    /**
     * @brief returns the predicted peak width at position mz
     * @param mz position of the peak
     */
    static DoubleReal getPeakWidth(DoubleReal mz);

    /**
     * @brief returns the quality of a potential feature at position act_mz by fitting the data to the isotope model
     * @param act_mz position of the potential feature
     */
    DoubleReal getPeakCorrelation(DoubleReal act_mz);



public:
	/**
	 * @brief double identifier (2)
	 */
    static const Int DOUBLE = 2;
    /**
     * @brief triple identifier (3)
     */
    static const Int TRIPLE = 3;
    /**
     * @brief returns the SILAC type of the filter. Either DOUBLE (2) or TRIPLE (3)
     */
    Int getSILACType();
    /**
     * @brief detailed constructor for double filtering
     * @param mass_separation_light_heavy distance between light and heavy peaks
     * @param charge_ charge of the ions to search for
     */
    SILACFilter(DoubleReal mass_separation_light_heavy,Int charge_,DoubleReal model_deviation_);
    /**
     * @brief detailed constructor for triple filtering
     * @param mass_separation_light_medium distance between light and medium peaks
     * @param mass_separation_light_heavy distance between light and heavy peaks
     * @param charge_ charge of the ions to search for
     */
	SILACFilter(DoubleReal mass_separation_light_medium,DoubleReal mass_separation_light_heavy,Int charge_,DoubleReal model_deviation_);
	/**
	 * @brief destructor
	 */
	virtual ~SILACFilter();
	/**
	 * @brief adds the given value to the blacklist
	 * @param the m/z value to add to the blacklist
	 */
	void blockValue(DoubleReal value);

	/**
	     * @brief returns if the given values is on the blacklist
	     * @param value m/z value to look up in the blacklist
	     */
	 bool blacklisted(DoubleReal value);

	/**
	 * @brief clears the blacklist
	 */
	void reset();
	/**
	 * @brief returns if there exists a SILAC feature at the given position, which corresponds to the filter's properties
	 * @param act_rt RT value of the position
	 * @param act_mz m/z value of the position
	 */
	bool isFeature(DoubleReal act_rt,DoubleReal act_mz);
	/**
	 * @brief gets the m/z values of all peaks , which belong the last identiefied feature
	 */
	std::vector<DoubleReal> getPeakValues();
	/**
	 * @brief returns the distance between light and medium peaks
	 */
	DoubleReal getEnvelopeDistanceLightMedium();
	/**
	 * @brief returns the distance between light and heavy peaks
	 */
	DoubleReal getEnvelopeDistanceLightHeavy();
	/**
	 * @brief returns the distance between two isotope peaks
	 */
	DoubleReal getIsotopeDistance();
	/**
	 * @brief returns all identiefied elements
	 */
	std::vector<DataPoint> getElements();
};
}

#endif /* SILACFILTER_H_ */
