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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DEFAULTPARAMHANDLER_H
#define OPENMS_DATASTRUCTURES_DEFAULTPARAMHANDLER_H

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief A base class for all classes handling default parameters.

		This class facilitates the handling of parameters:
		- it manages default parameter (defaults_)
		- it checks for valid parameters:
			- unknown/misspelled parameter names
			- correct parameter type
		  - range restrictions of numeric parameters
		  - valid values for string parameters (enum) 
		- subsections that are passed to other classes can be excluded from the check (subsections_)
		- it keeps member variables in syncronicity with the parameters stored in param_
		- it helps to automatically create a doxygen documentation page for the parameters

		Extra member variables are needed if getting the value from param_ would be too slow
		e.g. when they are used in methods that are called very often.

		No matter if you have extra variables or not, do the following:
		- Set defaults_ and subsections_ in the derived classes' default constructor.
		- Make sure to set the 'advanced' flag of the parameters right in order to hide certain parameters from unexperienced users.
		- Set the range restrictions for numeric defaults and valid strings for string defaults (if possible)
		- Call defaultsToParam_() at the end of derived classes' default constructor.
			It copies the defaults to param_ (and calls updateMembers_()).

		If you have extra member variables you need to syncronize with param_, do the following:
		- Implement the updateMembers_() method. It is used after each change of param_
			in order to update the extra member variables. If the base class is a DefaultParamHandler as well
			make sure to call the updateMembers_() method of the base class in the updateMembers_() method.
		- Call updateMembers_() at the end of the derived classes' copy constructor and assignment operator.
		- If you need mutable access to the extra member variables, provide a set-method and make sure to set
		  the corresponding value in param_ as well!
		
		@b Base @b classes: @n
		If you create a class @a A that is derived from DefaultParamHandler and derive another class @a B
		for @a A, you should set use the setName(String) method to set the name used for error messages to @a B.
		
		@b Parameter @b documentation: @n
		Each default parameter has to be documented in a comprehensive way. This is done using the
		Param::setValue methods and the Param::setDescription method.
		
		@b Flags: @n
		Flags (boolean parameters) are not supported directly. It's best to implement them as a
		string parameter with valid strings 'true' and 'false'.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI DefaultParamHandler
	{
		public:
			/// Constructor with name that is displayed in error messages
			DefaultParamHandler(const String& name);

			/// Copy constructor
			DefaultParamHandler(const DefaultParamHandler& rhs);

			/// Destructor
			virtual ~DefaultParamHandler();

			/// Assignment operator.
			virtual DefaultParamHandler& operator= (const DefaultParamHandler& rhs);

			/// Equality operator
			virtual bool operator== (const DefaultParamHandler& rhs) const;

			/**
				@brief Sets the parameters.

				Before setting the parameters, missing parameters are filled up with default values.
				@n Then the parameters are checked for unknown parameters (warning) and violations of restrictions (exception)
				with the @ref Param::checkDefaults() method.
				
				@exception Exception::InvalidParameter is thrown if errors occur during the check.
			*/
			void setParameters(const Param& param);

			/// Non-mutable access to the parameters
			const Param& getParameters() const;

			/// Non-mutable access to the default parameters
			const Param& getDefaults() const;

			/// Non-mutable access to the name
    	const String& getName() const;

			/// Mutable access to the name
    	void setName(const String& name);

    	/// Non-mutable access to the registered subsections
    	const std::vector<String>& getSubsections() const;

		protected:
			/**
				@brief This method is used to update extra member variables at the end of the setParameters() method.

				Also call it at the end of the derived classes' copy constructor and assignment operator.

				The default implementation is empty.
			*/
			virtual void updateMembers_();

			///Updates the parameters after the defaults have been set in the constructor
			void defaultsToParam_();

			///Container for current parameters
			Param param_;

			/**
				@brief Container for default parameters. This member should be filled in the constructor of derived classes!

				@note call the defaultsToParam_() method at the end of the constructor in order to copy the defaults to param_.
			*/
			Param defaults_;

			/**
				@brief Container for registered subsections. This member should be filled in the constructor of derived classes!

				@note Do not add a ':' character at the end of subsections.
			*/
			std::vector<String> subsections_;

			/// Name that is displayed in error messages during the parameter checking
			String error_name_;

			/**
				@brief If this member is set to false no checking if parameters in done;

				The only reason to set this member to false is that the derived class has no parameters!
			*/
			bool check_defaults_;

		private:
			/// Hidden default C'tor (class name parameter is required!)
			DefaultParamHandler();

	}; //class

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DEFAULTPARAMHANDLER_H
