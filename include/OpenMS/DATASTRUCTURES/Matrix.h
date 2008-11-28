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

#ifndef OPENMS_DATASTRUCTURES_MATRIX_H
#define OPENMS_DATASTRUCTURES_MATRIX_H

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cmath> // pow()
#include <iomanip>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

namespace OpenMS
{

  /**
     @brief A two-dimensional matrix.  Similar to std::vector, but uses a binary
     operator(,) for element access.

		 Think of it as a random access container.  You can also generate gray
     scale images.  This is not designed to be used for linear algebra.

		 The following member functions of the base class std::vector<ValueType>
		 can also be used:

		 - begin
		 - end
		 - rbegin
		 - rend
		 - front
		 - back
		 - assign
		 - empty
		 - size
		 - capacity
		 - max_size
		 .

		 (It seems that Doxygen does not parse pure access declarations, so we
		 list them here.)

		 @ingroup Datastructures
  */
  template <typename Value>
	class OPENMS_DLLAPI Matrix : protected std::vector < Value >
	{
	 protected:
		typedef std::vector < Value > Base;

	 public:

		///@name STL compliance type definitions
		//@{
		typedef Base container_type;

		typedef typename Base::difference_type difference_type;
		typedef typename Base::size_type size_type;

		typedef typename Base::const_iterator const_iterator;
		typedef typename Base::const_reverse_iterator const_reverse_iterator;
		typedef typename Base::iterator iterator;
		typedef typename Base::reverse_iterator reverse_iterator;

		typedef typename Base::const_reference const_reference;
		typedef typename Base::pointer pointer;
		typedef typename Base::reference reference;
		typedef typename Base::value_type value_type;

		typedef typename Base::allocator_type allocator_type;
		//@}

		///@name OpenMS compliance type definitions
		//@{
		typedef Base ContainerType;
		typedef difference_type DifferenceType;
		typedef size_type SizeType;

		typedef const_iterator ConstIterator;
		typedef const_reverse_iterator ConstReverseIterator;
		typedef iterator Iterator;
		typedef reverse_iterator ReverseIterator;

		typedef const_reference ConstReference;
		typedef pointer Pointer;
		typedef reference Reference;
		typedef value_type ValueType;

		typedef allocator_type AllocatorType;
		//@}

		///@name Constructors, assignment, and destructor
		//@{
		Matrix ()
			: Base(),
				rows_(0),
				cols_(0)
		{}

		Matrix (const SizeType rows, const SizeType cols, ValueType value = ValueType())
			: Base(rows*cols,value),
				rows_(rows),
				cols_(cols)
		{}

		Matrix (const Matrix & source)
			: Base(source),
				rows_(source.rows_),
				cols_(source.cols_)
		{}

		Matrix & operator = (const Matrix & rhs)
		{
			Base::operator=(rhs);
			rows_ = rhs.rows_;
			cols_ = rhs.cols_;
			return *this;
		}

		~Matrix() {}
		//@}

		///@name Accessors
		//@{
		const_reference operator() (size_type const i, size_type const j) const
		{
			return getValue(i,j);
		}

		reference operator() (size_type const i, size_type const j)
		{
			return getValue(i,j);
		}

    const_reference getValue(size_type const i, size_type const j) const
		{
			return Base::operator[](index(i,j));
		}

    reference getValue(size_type const i, size_type const j)
		{
			return Base::operator[](index(i,j));
		}

    void setValue(size_type const i, size_type const j, value_type value)
		{
			Base::operator[](index(i,j)) = value;
		}

		//@}


		/**
			 @name Pure access declarations

			 These make begin(), end() etc. from container_type accessible.
		*/
		//@{
	 public:

		Base::begin;
		Base::end;
		Base::rbegin;
		Base::rend;

		Base::front;
		Base::back;
		Base::assign;

		Base::empty;
		Base::size;

		Base::capacity;
		Base::max_size;

		//@}

		void clear()
		{
			Base::clear();
			rows_ = 0;
			cols_ = 0;
		}

		void resize(size_type i, size_type j, value_type value = value_type())
		{
			rows_ = i;
			cols_ = j;
			Base::resize(rows_*cols_, value);
		}

		void resize(std::pair<UInt,UInt> const & size_pair, value_type value = value_type())
		{
			rows_ = size_pair.first;
			cols_ = size_pair.second;
			Base::resize(rows_*cols_, value);
		}

		/// Number of rows
		SizeType rows() const 
		{
			return rows_;
		}

		/// Number of columns
		SizeType cols() const 
		{
			return cols_;
		}

		std::pair<UInt,UInt> sizePair() const
		{
			return std::pair<UInt,UInt>(rows_,cols_);
		}

		/**@brief Calculate the index into the underlying vector from row and
			 column.  Note that Matrix uses the (row,column) lexicographic ordering
			 for indexing.
		*/
		SizeType const index(SizeType row, SizeType col) const
		{
#ifdef OPENMS_DEBUG
			if ( row >= rows_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,row, rows_);
			if ( col >= cols_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,col, cols_);
#endif
			return row * cols_ + col;
		}

		/**@brief Calculate the row and column from an index into the underlying
			 vector.  Note that Matrix uses the (row,column) lexicographic ordering
			 for indexing.
		*/
		std::pair<UInt,UInt> const indexPair(UInt index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return std::pair<SizeType,SizeType>(index/cols_,index%cols_);
		}

		/**@brief Calculate the column from an index into the underlying vector.
			 Note that Matrix uses the (row,column) lexicographic ordering for
			 indexing.
		*/
		SizeType colIndex(SizeType index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return index%cols_;
		}

		/**@brief Calculate the row from an index into the underlying vector.
			 Note that Matrix uses the (row,column) lexicographic ordering for
			 indexing.
		*/
		SizeType rowIndex(SizeType index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif

			return index/cols_;
		}

		/**@brief Equality comparator.

		If matrices have different row or colmn numbers, throws a precondition exception.
		*/
		bool operator == ( Matrix const & rhs ) const
			throw (Exception::Precondition)
		{
			OPENMS_PRECONDITION(cols_ == rhs.cols_,
													"Matrices have different row sizes.");
			OPENMS_PRECONDITION(rows_ == rhs.rows_,
													"Matrices have different column sizes.");
 			return static_cast < typename Matrix<Value>::Base const &>(*this) == static_cast < typename Matrix<Value>::Base const &>(rhs);
		}

		/**@brief Less-than comparator.  Comparison is done lexicographically: first by row, then by column.

		If matrices have different row or column numbers, throws a precondition exception.
		*/
		bool operator < (Matrix const & rhs) const
			throw (Exception::Precondition)
		{
			OPENMS_PRECONDITION(cols_ == rhs.cols_,
													"Matrices have different row sizes.");
			OPENMS_PRECONDITION(rows_ == rhs.rows_,
													"Matrices have different column sizes.");
			return static_cast < typename Matrix<Value>::Base const &>(*this) < static_cast < typename Matrix<Value>::Base const &>(rhs);
		}

		/// set matrix to 2D arrays values
		template <int ROWS, int COLS>
		void setMatrix (const ValueType matrix[ROWS][COLS])
		{
			resize(ROWS, COLS);
			for (SizeType i=0; i<rows_; ++i)
			{
				for (SizeType j=0; j<cols_; ++j)
				{
					setValue(i,j,matrix[i][j]);
				}	
			}
		}

		/**
		 *	allocate and return an equivalent GSL matrix
		 *  @note Works only for Matrix<double>
		 *	@note Clean up the gsl_matrix using gsl_matrix_free (gsl_matrix * m)
		 */
		gsl_matrix* toGslMatrix();
		
	protected:

		///@name Data members
		//@{
		/// Number of rows (height of a column)
		SizeType rows_;
		/// Number of columns (width of a row)
		SizeType cols_;
		//@}

	}; // class Matrix

	template<>
	gsl_matrix* Matrix<double>::toGslMatrix();
	
	/**@brief Print the contents to a stream.

	@relatesalso Matrix
	*/
	template <typename Value>
	std::ostream& operator << (std::ostream& os, const Matrix<Value>& matrix)
	{
		typedef typename Matrix<Value>::size_type size_type;
	  for ( size_type i = 0; i < matrix.rows(); ++i ) {
	    for ( size_type j = 0; j < matrix.cols(); ++j ) {
				os << std::setprecision(6) << std::setw(6) << matrix(i,j) << ' ';
	    }
	    os << std::endl;
	  }
		return os;
	}


} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_MATRIX_H
