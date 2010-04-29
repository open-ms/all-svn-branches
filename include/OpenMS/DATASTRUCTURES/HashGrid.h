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


#include<map>
#include<list>
#include<set>
#include<vector>
#include<iostream>

#include <OpenMS/DATASTRUCTURES/GridElement.h>



#ifndef HASHGRID_H_
#define HASHGRID_H_

namespace OpenMS
{

/**
		@brief A datastructure, which allows the arrangement of data points with an RT and m/z value in a two-dimensional grid.
		The size of each grid cell is determined by two values, namely <i>rt_threshold</i> and <i>mz_threshold</i>.
		<i>rt_threshold</i> defines the height of a grid cell and <i>mz_threshold</i> the width.
		The data points are stored in specific grid cells and are accessible via geometric hashing.
		So the corresponding cell of each data point can be calculated by dividing the rt and m/z values by its corresponding threshold.
		@ingroup Datastructures
	*/

class HashGrid {

private:

	DoubleReal rt_threshold;
	DoubleReal mz_threshold;
	DoubleReal hyp_threshold;
	DoubleReal rt_scaling;
	Int grid_size_x;
	Int grid_size_y;
	Int number_of_elements;

public:

	std::map<std::pair<Int,Int>, std::list<GridElement*> > elements;

/** @brief detailed constructor

	@param rt_threshold_ defines the height of each grid cell
	@param mz_threshold_ defines the width of each grid cell
*/
	HashGrid(DoubleReal rt_threshold_,DoubleReal mz_threshold_);
//Destructor
	~HashGrid();

/** @brief removes an element from the hash grid. The cell, in which the element may be contained, is specified:

	@param element_ the element to remove
	@param x x-value of the cell
	@param y y-value of the cell
*/
	void removeElement(GridElement* element_,Int x,Int y);
/** @brief removes an element from the hash grid. The cell, in which the element may be contained, is calculated out of the RT and m/z values of the element:

	@param element_ the element to remove
*/
	void removeElement(GridElement* element_);
/** @brief inserts a new element in the grid:

	@param element_ the element to remove
*/
	void insert(GridElement* element_);
/** @brief writes the content of the grid to the console:

*/
	void consoleOut();
/** @brief gets the number of cells

*/
	int size();
/** @brief gets the height of the cells

*/
	DoubleReal getRTThreshold() const;
/** @brief gets the width of the cells

*/
	DoubleReal getMZThreshold() const;
/** @brief gets the number of grids in m/z-direction

		*/
	Int getGridSizeX();

/** @brief gets the number of grids in RT-direction

		*/
	Int getGridSizeY();
	/**
	 * @brief gets the number of elements in the grid
	 */
	Int getNumberOfElements();
};
}




#endif /* HASHGRID_H_ */
