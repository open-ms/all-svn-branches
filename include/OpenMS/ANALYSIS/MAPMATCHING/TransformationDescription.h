// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <algorithm>

namespace OpenMS
{
	/**
	@brief Generic description of a coordinate transformation.
		
	This description stores the transformation name and parameters.
		
	The transformation can be applied to a double coordinate by using the @ref apply() method.
	@n Available transformations and parameters are:
	- none : \f$ f(x) = x \f$
	- linear : \f$ f(x) = \textit{intercept} + \textit{slope} * x \f$
	- interpolated_linear : Linear interpolation between pairs, extrapolation using first and last pair. At least two pairs must be given.
	- interpolated_cubic : Like interpolated_linear but cubic (?). <i>Not yet implemented!</i>
		
	You can also use transformation names not listed above, but you cannot apply
	the tranformation using this class then.
		
	Additionally corresponding coordinate pairs can be stored, e.g.
	to describe transformations that cannot be expressed as a simple function.
	When storing the pairs, but no function parameters, the name 'pairs' should be used.
		
	@ingroup MapAlignment
	*/
	class OPENMS_DLLAPI TransformationDescription
	{
		
	 public:
			
		///Coordnate pair vector type
		typedef std::vector< std::pair<Real,Real> > PairVector;
			
		/// Constructor
		TransformationDescription();
		/// Destructor
		~TransformationDescription();
				
		/// Copy constructor 
		TransformationDescription(const TransformationDescription& rhs);
		/// Assignment operator
		TransformationDescription& operator = (const TransformationDescription& rhs);
				
		/// Resets everything
		void clear();
			
		///Returns the name
		const String& getName() const
		{
			return name_;
		}
		///Sets the name
		void setName(const String& name)
		{
			delete trafo_;
			trafo_ = 0;
			name_ = name;
		}
			
		///Non-mutable access to the parameters
		const Param& getParameters() const
		{
			return param_;
		}
			
		///Sets the name
		void setParameters(const Param& param)
		{
			delete trafo_;
			trafo_ = 0;
			param_ = param;
		}
			
		///Returns the pairs
		const PairVector& getPairs() const
		{
			return pairs_;
		}
			
		///Sets the pairs
		void setPairs(const PairVector& pairs)
		{
			pairs_ = pairs;
		}
			
		/**
		@brief Convenience method to access double parameters
					
		@exception Exception::ElementNotFound is thrown if the parameter does not exist.
		*/
		DoubleReal getParam(const String& name) const
		{
			return param_.getValue(name);
		}
				
		/// Convenience method to set double parameters
		void setParam(const String& name, DoubleReal value)
		{
			delete trafo_;
			trafo_ = 0;
			param_.setValue(name,value);
		}
				
		/**
		@brief Apply the transformation to @p value .
					 
		@exception Exception::IllegalArgument is thrown if the transformation cannot be initialized according to the given name and parameters.
		*/
		void apply(DoubleReal& value)
		{
			//initialize transformation (if unset)
			if (!trafo_) init_();
			//apply transformation
			trafo_->operator()(value);
		}
				
	 protected:
			
		///Base class for all transformations
		struct Trafo_
		{
			virtual void operator ()(DoubleReal& value) const = 0;
		};
			
		///Tranformation name
		String name_;
		///Tranformation parameters
		Param param_;
		///Pairs of corrensponding values
		PairVector pairs_;
		///Poiter to actual transformation functor
		Trafo_ * trafo_;
				
		/// No transformation (i.e. identity)
		struct None_ : Trafo_
		{
			None_()
			{
			}
			virtual void operator ()(DoubleReal& ) const
			{
			}
		};
			
		///@brief Linear transformation that actually applies an affine transformation ;-)
		struct Linear_ : Trafo_
		{
			Linear_(DoubleReal slope, DoubleReal intercept)
				: slope_(slope),
					intercept_(intercept)
			{
			}
			virtual void operator ()(DoubleReal& value) const
			{
				value *= slope_;
				value += intercept_;
			}
			DoubleReal slope_;
			DoubleReal intercept_;
		};
		
		///@brief Piecewise linear transformation.
		/**
		In between the pairs, the interpolation uses the neighboring pairs.
		Outside the range spanned by the pairs, we extrapolate using a line
		through the first and the last pair.

		Each time this is applied, a binary search is performed.  We could
		precompute slopes for each segment, but it is not clear if this will pay
		off.
		*/
		struct InterpolatedLinear_ : Trafo_
		{
			InterpolatedLinear_(PairVector const& pairs)
				: pairs_(pairs)
			{
				if ( pairs_.size() < 2 )
				{
					throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"less than two pairs for 'interpolated_linear' transformation given");
				}
				std::sort(pairs_.begin(),pairs_.end());
			}
			virtual void operator ()(DoubleReal& value) const
			{
				if ( value <= pairs_.front().first )
				{
					DoubleReal slope = ( pairs_.back().second - pairs_.front().second ) / ( pairs_.back().first - pairs_.front().first );
					value = pairs_.front().second + ( value - pairs_.front().first ) * slope;
				}
				else if ( value >= pairs_.back().first )
				{
					DoubleReal slope = ( pairs_.back().second - pairs_.front().second ) / ( pairs_.back().first - pairs_.front().first );
					value = pairs_.back().second + ( value - pairs_.back().first ) * slope;
				}
				else
				{
					PairVector::value_type value_as_pair(value,0);
					PairVector::const_iterator right = std::lower_bound( pairs_.begin(), pairs_.end(), value_as_pair );
					PairVector::const_iterator left = right;
					--left;
					DoubleReal slope = ( right->second - left->second ) / ( right->first - left->first );
					value = left->second + ( value - left->first ) * slope;
				}
				return;
			}
			PairVector pairs_;
		};
		
		/**
		@brief Initialize the transformation according to the name and parameters.
		
		@exception Exception::IllegalArgument is thrown if the transformation cannot be initialized according to the name and parameters.
		*/
		void init_()
		{
			if ( trafo_ ) delete trafo_;
			trafo_ = 0;
			if (name_=="none")
			{
				trafo_ = new None_();
			}
			else if (name_=="linear")
			{
				if (!param_.exists("slope"))
				{
					throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"parameter 'slope' for 'linear' transformation not given");
				}
				if (!param_.exists("intercept"))
				{
					throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"parameter 'intercept' for 'linear' transformation not given");
				}
				trafo_ = new Linear_(param_.getValue("slope"),param_.getValue("intercept"));
			}
			else if (name_=="interpolated_linear")
			{
				trafo_ = new InterpolatedLinear_(pairs_);
			}
			else
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,(String("unknown transformation name '") + name_ + "'").c_str());
			}
		}
			
	};

	std::ostream& operator<<(std::ostream& os, TransformationDescription const & td);
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H

