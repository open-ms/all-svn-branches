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
#include <queue>
#include <list>

namespace OpenMS
{
/**
 * @brief Filter to use for SILACFiltering
 *
 * A SILACFilter searches for SILAC patterns, which correspond to the defined mass shifts and charge.
 * Only peaks are taken into account, which were not blacklisted by other filters before e.g. are not part of
 * a SILAC pair yet.
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
		 * @brief number of isotopes per peptide to search for
		 */
		Int isotopes_per_peptide;
    /**
     * @brief envelope distances within one feature
     */
    std::set<DoubleReal> envelope_distances;

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

    typedef std::multiset<DoubleReal,doubleCmp> Blacklist;

    /**
     * @brief holds the m/z values, which are ignored due to blacklisting of other filters
     */
    Blacklist blacklist;

    /**
     * @brief defines the lifetime of each position in the blacklist.
     * The positions with longest lifetime, i.e. whose iterators are located at the front of the list, will be deleted at the next reset
     */
    std::list<std::list<Blacklist::iterator> > blacklist_lifetime;
    /**
     * @brief holds the m/z values of the last identified feature. These values will be blacklisted in other filters
     */
    std::vector<DoubleReal> peak_positions;

    /**
     * @brief returns the predicted peak width at position mz
     * @param mz position of the peak
     */
    static DoubleReal getPeakWidth(DoubleReal mz);

    /**
     * @brief Computes the cross correlation of the area around current_mz to the spectrum in a given area
     * @param current_mz position of the potential feature
     * @param offset expected distance between monoisotpic peak and current peak
     * @param tolerance maximal deviation from the expected distance
     * @param data is filled during the computation
     */
    void computeCorrelation(DoubleReal current_mz,DoubleReal offset,DoubleReal tolerance,std::vector<DoubleReal>& data);

    /**
     * @brief Computes the exact position of a peak to the monoisotopic peak
     * The computation is based on the expected distance to the monoisotopic peak and their autocorrelation
     * @param current_mz m/z position of the monoisotopic peak
     * @param expected_distance expected distance of the current peak to the monoisotopic peak
     * @param tolerance maximal deviation from the expected distance
     * @param data autocorrelation data vector as calculated in computeCorrelation()
     */
    DoubleReal computeExactDistance(DoubleReal current_mz,DoubleReal expected_distance,DoubleReal tolerance,std::vector<DoubleReal> data);

/**
 * @brief Determines the quality of an isotope pattern by computing the Pearson correlation and the averagine model deviation
 * @param current_mz current m/z position
 * @param exact_positions the distances of each peak to the monoisotopic peak
 * @param intensities vector to be filled with the intensities of each peak
 * @param missing_peak is true if already a peak is missing in the SILAC pattern
 */
    bool checkPattern(DoubleReal current_mz, const std::vector<DoubleReal>& exact_positions_heathrow, std::vector<DoubleReal>& intensities, bool missing_peak);

    /*
    bool checkRatios(DoubleReal current_mz,const std::vector<DoubleReal>& light_positions, const std::vector<DoubleReal>& envelope_positions);
	*/


	/**
	 * @brief returns if there exists a SILAC feature at the given position, which corresponds to the filter's properties
	 * @param current_rt RT value of the position
	 * @param current_mz m/z value of the position
	 */
	bool isPair(DoubleReal current_rt,DoubleReal current_mz);
	/**
	 * @brief gets the m/z values of all peaks , which belong the last identiefied feature
	 */
	std::vector<DoubleReal> getPeakPositions();
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
     * @brief detailed constructor for SILAC pair filtering
     * @param mass_separations all mass shifts of the filter
     * @param charge_ charge of the ions to search for
     * @param model_deviation_ maximum deviation from the averagine model
		 * @param isotopes_per_peptide_ number of isotopes per petide to search for
     */
    SILACFilter(std::set<DoubleReal> mass_separations,Int charge_,DoubleReal model_deviation_, Int isotopes_per_peptide_);

    /**
         * @brief detailed constructor for singlet filtering
         * No mass shifts are given, so only singlets are searched.
         * @param charge_ charge of the ions to search for
         * @param model_deviation_ maximum deviation from the averagine model
         */
    SILACFilter(Int charge_,DoubleReal model_deviation_);


	/**
	 * @brief destructor
	 */
	virtual ~SILACFilter();
	/**
	 * @brief adds the given value to the blacklist
	 * @param the m/z value to add to the blacklist
	 */
	void  blockValue(DoubleReal value);

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
	 * @brief returns the distance between light and heavy peaks
	 */
	std::set<DoubleReal> getEnvelopeDistances();
	/**
	 * @brief returns the distance between two isotope peaks
	 */
	DoubleReal getIsotopeDistance();
	/**
	 * @brief returns all identiefied elements
	 */
	std::vector<DataPoint> getElements();

	/**
	 * @brief returns the charge of the filter
	 */
	Int getCharge();

  /**
   * @brief returns the number of isotopes per peptide of the filter
   */
  Int getIsotopesPerPeptide();
};
}

#endif /* SILACFILTER_H_ */
