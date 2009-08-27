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
// $Maintainer: Pavel Kudan $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <Isotope.H>
#include <OpenMS/CONCEPT/Types.h>

#ifndef SET_NUMBER_OF_NEUTRONS
#define SET_NUMBER_OF_NEUTRONS(pelement, mass) \
{\
        SignedSize number_of_neutrons =  ceil((mass - (Constants::PROTON_MASS_U + Constants::ELECTRON_MASS_U) * pelement->getAtomicNumber())/Constants::NEUTRON_MASS_U);\
        if (number_of_neutrons < 0)\
        {\
                throw Exception::SizeUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);\
        }\
        else\
        {\
                number_of_neutrons_=(Size) number_of_neutrons;\
        }\
}
#endif

#ifndef SET_DEFAULT_NAME_AND_SYMBOL
#define SET_DEFAULT_NAME_AND_SYMBOL \
{\
        Size integer_mass = number_of_neutrons_ + element_->getAtomicNumber();\
        String element_name = element_->getName();\
        name_= element_name;\
        name_ += "-";\
        name_ += integer_mass;\
        symbol_= integer_mass;\
        symbol_ += element_name;\
}
#endif



using namespace std;

namespace OpenMS
{
        Isotope::Isotope()
                :	element_(NULL),
                        mass_(0.0),
                        number_of_neutrons_(0),
                        name_(""),
                        symbol_("")
	{
	}

        Isotope::Isotope(const Isotope& isotope)
                :	element_(isotope.element_),
                        mass_(isotope.mass_),
                        number_of_neutrons_(isotope.number_of_neutrons_),
                        name_(isotope.name_),
                        symbol_(isotope.symbol_)
	{
//todo: check if number of neutrons are correct
	}

        Isotope::Isotope( Element* pelement,
                                DoubleReal mass)

                :	element_(pelement),
                        mass_(mass),
                        name_(""),
                        symbol_("")
        {
                if (pelement!=NULL && mass!=0.0)
                {
                        SET_NUMBER_OF_NEUTRONS(element_, mass_);
                        SET_DEFAULT_NAME_AND_SYMBOL;
                }
        }
        Isotope::Isotope( Element* pelement,
                                DoubleReal mass,
                                String name,
                                String symbol)

                :	element_(pelement),
                        mass_(mass),
                        name_(name),
                        symbol_(symbol)
	{
                if (pelement!=NULL && mass!=0.0)
                {
                        SET_NUMBER_OF_NEUTRONS (element_, mass_);
                }
	}

        Isotope::~Isotope()
	{
	}

//        void Isotope::setNumberOfNeutrons(UInt number_of_neutrons)
//	{
//                number_of_neutrons_ = number_of_neutrons;
//	}
	
        UInt Isotope::getNumberOfNeutrons() const
	{
                return number_of_neutrons_;
	}

//todo: clear()
        void Isotope::setMass(DoubleReal mass)
	{
                if (element_!=NULL && mass!=0.0)
                {
                        SET_NUMBER_OF_NEUTRONS (element_, mass);
                        SET_DEFAULT_NAME_AND_SYMBOL;
                }
                else
                {
                        number_of_neutrons_ = 0;
                        name_="";
                        symbol_="";
                }
                mass_ = mass;
	}
	
        DoubleReal Isotope::getMass() const
	{
                return mass_;
	}

        void Isotope::setElement(Element* pelement)
        {
                if (pelement!=NULL && mass_!=0.0)
                {
                        SET_NUMBER_OF_NEUTRONS (pelement, mass_);
                        SET_DEFAULT_NAME_AND_SYMBOL;
                }
                else
                {
                        number_of_neutrons_ = 0;
                        name_="";
                        symbol_="";
                }
                //todo: change name and symbol if changed
                element_ = pelement;
	}
	
        Element* Isotope::getElement() const
	{
                return element_;
	}

        void Isotope::setName(const String& name)
	{
		name_ = name;
	}
	
        String Isotope::getName() const
	{
                return name_;
	}

        void Isotope::setSymbol(const String& symbol)
	{
		symbol_ = symbol;
	}
	
       String Isotope::getSymbol() const
	{
                return symbol_;
	}

        Isotope& Isotope::operator = (const Isotope& isotope)
	{
                element_ = isotope.element_;
                mass_ = isotope.mass_;
                number_of_neutrons_ = isotope.number_of_neutrons_;
                symbol_ = isotope.symbol_;
                name_ = isotope.name_;
		return *this;
	}

        bool Isotope::operator == (const Isotope& isotope) const
	{
                //name or symbol may differs, but still, if element and isotope mass
                //are the same - it is the same isotope
                return 	mass_ == isotope.mass_ &&
                                        element_ == isotope.element_;
        }

        bool Isotope::operator != (const Isotope& isotope) const
	{
                return !(*this == isotope);
	}
	
        std::ostream& operator << (std::ostream& os, const Isotope& isotope)
	{
                os 	<< isotope.name_ << " "
                                << isotope.symbol_ << " "
                                << isotope.number_of_neutrons_ << " "
                                << isotope.mass_  << " "
                                << isotope.element_->getName();

		return os;
	}


} // namespace OpenMS

