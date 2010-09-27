// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow $ 
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_COMPOMER_H
#define OPENMS_DATASTRUCTURES_COMPOMER_H

#include <vector>
#include <map>
#include <set>
#include <cstdlib>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

namespace OpenMS
{

/**
	@brief Holds information on an egde connecting two features from a (putative) charge ladder
	
	A compomer is storing information on the adducts used on LEFT and RIGHT nodes (Features) that are connected by the egde (i.e. ChargePair)
	holding the compomer. Additonally meta information on the edge (net_charge, edge score, id)
	which is kept up-to-date when adducts are added to either side is stored.
	
*/
class OPENMS_DLLAPI Compomer
{
public:
	/// side of compomer (LEFT ^ substract; RIGHT ^ add)
	enum SIDE {LEFT, RIGHT, BOTH};
	
	typedef Map<String,Adduct> CompomerSide; /// adducts and their abundance etc
	typedef std::vector< CompomerSide > CompomerComponents; /// container for the two sides [0]=left, [1]=right
	
	/// Default Constructor
	Compomer();

	/// Constructor with net-charge and mass
	Compomer(Int net_charge, DoubleReal mass, DoubleReal log_p);
	
	/// Copy C'tor
	Compomer(const Compomer& p);
	
	/// Assignment Operator
	Compomer& operator=(const Compomer& source);
	
	/// Add a.amount of Adduct @param a to Compomer's @param side and update its properties
	void add(const Adduct& a, UInt side);
	
	/**
	 *  indicates if these two compomers can coexist for one feature
	 * @param cmp The other Compomer we compare to
	 * @param side_this Indicates which "side"(negative or positive adducts) we are looking at. Negative adducts belong to the left side of the ChargePair.
	 * @param side_other See above.
	 */
	bool isConflicting(const Compomer& cmp, UInt side_this, UInt side_other) const;

	/// set an Id which allows unique identification of a compomer
	void setID(const Size& id);
	/// return Id which allows unique identification of this compomer
	const Size& getID() const;
	/// left and right adducts of this compomer
	const CompomerComponents& getComponent() const;

	/// net charge of compomer (i.e. difference between left and right side of compomer)
	const Int& getNetCharge() const;

	/// mass of all contained adducts
	const DoubleReal& getMass() const;

	/// summed positive charges of contained adducts
	const Int& getPositiveCharges() const;
	
	/// summed negative charges of contained adducts
	const Int& getNegativeCharges() const;
	
	/// return log probability
	const DoubleReal& getLogP() const;

	/// return log probability
	const DoubleReal& getRTShift() const;

	/// get adducts with their abundance as compact string for both sides
	String getAdductsAsString() const;

	/// get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
	/// @param side Use LEFT for left, RIGHT for right
	String getAdductsAsString(UInt side) const;
	
	/// check if Compomer only contains a single adduct on side @p side
	bool isSingleAdduct(Adduct& a, const UInt side) const;
	
	/**
		@brief Remove all adducts of type @p a
	
		Remove ALL instances of the given adduct,
		BUT use the given adducts parameters (charge, logp, mass etc) to update the compomers members
	**/
	Compomer removeAdduct(const Adduct& a) const;
	
	/**
		@brief Remove all adducts of type @p a from @p side (LEFT or RIGHT)
		
		remove ALL instances of the given adduct from the given side (LEFT or RIGHT),
		BUT use the given adducts parameters (charge, logp, mass etc) to update the compomers members
	*/
	Compomer removeAdduct(const Adduct& a, const UInt side) const;
	
	/**
		@brief Returns the adduct labels from @p side (LEFT or RIGHT)
		
		Get a list of labels for the @p side
		(useful for assigning channels (light, heavy etc) to features).
	*/
	StringList getLabels(const UInt side) const;
	
	
	/// Adds @p add_side to this compomer.
	void add(const CompomerSide& add_side, UInt side);
	
	/// Sort compomer by (in order of importance): net-charge, mass, probability
	friend OPENMS_DLLAPI bool operator< (const Compomer &c1, const Compomer &c2); 

	/// Print the contents of a Compomer to a stream.
	friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Compomer & cmp);

	/// Comparator
	friend OPENMS_DLLAPI bool operator==(const Compomer& a, const  Compomer& b);
	
private:	
	
	CompomerComponents cmp_; //< adducts of left and right side
	Int net_charge_;		 //< net charge (right - left)
	DoubleReal mass_;		 //< net mass (right - left)
	Int pos_charges_;		 //< net charges on the right
	Int neg_charges_;		 //< net charges on the left
	DoubleReal log_p_;   //< log probability of compomer
	DoubleReal rt_shift_; //< expected net RT shift of compomer (-shift_leftside + shift_rightside)
	Size id_;

}; // \Compomer

} // namespace OpenMS

#endif //OPENMS_DATASTRUCTURES_COMPOMER_H
