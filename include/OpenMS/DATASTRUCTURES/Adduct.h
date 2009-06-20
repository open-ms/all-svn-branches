
// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_ANALYSIS_DATASTRUCTURES_ADDUCT_H
#define OPENMS_ANALYSIS_DATASTRUCTURES_ADDUCT_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS {
    
class OPENMS_DLLAPI Adduct
{
public:

		typedef std::vector< Adduct > AdductsType;

    Adduct()
    : charge_(0),
    amount_(0),
    singleMass_(0),
    log_prob_(0),
    formula_()
    {
    }

    Adduct(Int charge)
    : charge_(charge),
    amount_(0),
    singleMass_(0),
    log_prob_(0),
    formula_()
    {
    }

    Adduct(Int charge, Int amount, DoubleReal singleMass, String formula, DoubleReal log_prob)
    : charge_(charge),
    amount_(amount),
    singleMass_(singleMass),
    log_prob_(log_prob),
    formula_(formula)
    {
    }

    Adduct operator *(Int m)
    {
        Adduct a = *this;
        a.amount_ *= m;
        return a;
    }
    Adduct operator +(const Adduct& rhs)
    {
        if (this->formula_ != rhs.formula_)
        {
          throw "Adduct::Operator +()  tried to add incompatible adduct!";
        }
        Adduct a = *this;
        a.amount_ += rhs.amount_;
        return a;
    }

    void operator +=(const Adduct& rhs)
    {
      if (this->formula_ != rhs.formula_)
      {
        throw "Adduct::Operator +=()  tried to add incompatible adduct!";
      }
      this->amount_ += rhs.amount_;
    }

		
		/// Print the contents of an Adduct to a stream.
    friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Adduct& a);

		/// Comparator
		friend OPENMS_DLLAPI bool operator==(const Adduct& a, const Adduct& b);
	
    //@{ Accessors
    const Int& getCharge() const
    {
      return charge_;
    }
    
    void setCharge(const Int& charge)
    {
      charge_ = charge;
    }

    const Int& getAmount() const
    {
      return amount_;
    }
    
    void setAmount(const Int& amount)
    {
      amount_ = amount;
    }    
 
    const DoubleReal& getSingleMass() const
    {
      return singleMass_;
    }
    
    void setSingleMass(const DoubleReal& singleMass)
    {
      singleMass_ = singleMass;
    }       

    const DoubleReal& getLogProb() const
    {
      return log_prob_;
    }
    
    void setLogProb(const DoubleReal& log_prob)
    {
      log_prob_ = log_prob;
    }  

    const String& getFormula() const
    {
      return formula_;
    }
    void setFormula(const String& formula)
    {
      formula_ = formula;
    }      
//}
    
private:
    Int charge_; //usually +1
    Int amount_; // number of entities (can be negative!)
    DoubleReal singleMass_; //mass of a single entity
    DoubleReal log_prob_;   // log probability of observing a single entity of this adduct
    String formula_;   // chemical formula (parsable by EmpiricalFormula)

};

} // namespace OpenMS


#endif

