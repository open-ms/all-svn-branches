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


#ifndef OPENMS_COMPARISON_CLUSTERING_QTCLUSTERING_H
#define QTCLUSTERING_H_

#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/QTSILACCluster.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS {
/**
 * @brief QT clustering based on geometric hashing.
 * Performs a QT clustering similar to the description of Heyer, Kruglyak and Yooseph (1999).
 * It uses a hash grid for the arrangement of the data points and computes a two-dimensional diameter for cluster (m/z-diameter, rt-diameter).
 * @see HashGrid
 * @ingroup SpectraClustering
 */
class OPENMS_DLLAPI QTClustering  : public ProgressLogger{
private:
	/**
	 * @brief the hash grid used for data arrangement
	 */
	HashGrid grid;
	/**
	 * @brief maximal rt diameter
	 * corresponds to the maximal gap in RT direction of cluster
	 */
	DoubleReal rt_diameter;
	/**
	 * @brief maximal m/z diameter
	 * corresponds to the maximal cluster extent in m/z direction
	 */
	DoubleReal mz_diameter;
	/**
	 * @brief list of identified clusters
	 */
	std::list<QTSILACCluster> clusters;
	/**
	 * @brief default constructor
	 */
	QTClustering();

	/**
	 * @brief recursive QT clustering method
	 * @param act_grid the data points to be clustered in the current step
	 */
	QTSILACCluster QTClust(HashGrid& act_grid);
public:
	/**
	 * @brief detailed constructor
	 * @param data the data to be clustered
	 * @param rt_diameter_ maximal rt diameter
	 * @param mz_diameter_ maximal m/z diameter
	 */
	QTClustering(std::vector<DataPoint>& data,DoubleReal rt_diameter_, DoubleReal mz_diameter_);
	/**
	 * @brief destructor
	 */
	virtual ~QTClustering();
	/**
	 * @brief performs the clustering on the given data and diameters and returns a vector of clusters
	 */
	std::vector<std::vector<DataPoint*> > performClustering();

	/**
			@brief Exception thrown if not enough data (<2) is used

			If the set of data to be clustered contains only one data point,
			clustering algorithms would fail for obvious reasons.
	 */
	class OPENMS_DLLAPI InsufficientInput : public Exception::BaseException
	{
	public:
		InsufficientInput(const char* file, int line, const char* function, const char* message= "not enough data points to cluster anything") throw();
		virtual ~InsufficientInput() throw();
	};
};
}
#endif /* QTCLUSTERING_H_ */
