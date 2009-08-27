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
// $Authors: Pavel Kudan $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_ISOTOPE_H
#define OPENMS_CHEMISTRY_ISOTOPE_H


#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CONCEPT/Constants.h>


namespace OpenMS
{
	/** @ingroup Chemistry
	
                        @brief Representation of an isotope
	*/
        class OPENMS_DLLAPI Isotope
	{
		public:

			/** @name Constructor and Destructors
			*/
			//@{
			/// default constructor
                        Isotope();

			/// copy constructor
                        Isotope(const Isotope& isotope);

                        /// detailed constructor with default name and symbol and calculated number_of_neutrons
                        Isotope(Element* pelement,
                                DoubleReal mass);

			/// detailed constructor
                        Isotope(Element* pelement,
                                DoubleReal mass,
                                UInt number_of_neutrons
                                String name,
                                String symbol);

			/// destructor
                        virtual ~Isotope();
			//@}

			/** @name Accessors
			*/
			//@{
                        /// sets the pointer to the element to which isotope belongs
                        void setElement(Element* pelement);

                        /// returns the pointer to the element to which isotope belongs
                        Element* getElement() const;

  //                      /// sets the unique number of neutrons
  //                      void setNumberOfNeutrons(UInt number_of_neutrons);

                        /// returns the unique number of neutrons
                        UInt getNumberOfNeutrons() const;
			
                        /// sets exact mass of the isotope
                        void setMass(DoubleReal mass);

                        /// returns the exact mass of the isotope
                        DoubleReal getMass() const;

                        /// set the name of the isotope
                        void setName(const String& name);

                        /// returns the name of the isotope
                        String getName() const;

                        /// sets symbol of the isotope
                        void setSymbol(const String& symbol);

                        /// returns symbol of the isotope
                        String getSymbol() const;




			//@}

			/** @name Assignment
			*/
			//@{
			/// assignment operator
                        Isotope& operator = (const Isotope& isotope);
			//@}

			/** @name Predicates
			*/
			//@{
			/// equality operator
                        bool operator == (const Isotope& isotope) const;

			/// inequality operator
                        bool operator != (const Isotope& isotope) const;
			//@}

                        /// writes the isotope to an output stream
                        friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Isotope& isotope);

		protected:
                        /// the pointer to the element to which the isotope belongs
                        Element* element_;

                        /// exact mass of isotope
                        DoubleReal mass_;

                        /// the unique number of neutrons of the isotope
                        UInt number_of_neutrons_;

                        /// name of the isotope
                        String name_;

                        /// symbol of the isotope
                        String symbol_;

        };

        OPENMS_DLLAPI std::ostream& operator << (std::ostream&, const Isotope&);

} // namespace OpenMS

#endif

