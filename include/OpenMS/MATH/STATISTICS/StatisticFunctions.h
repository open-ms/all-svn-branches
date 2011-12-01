// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <numeric>
#include <algorithm>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H
#define OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H

namespace OpenMS
{

	namespace Math
	{
		/**
			 @brief Calculates the sum of a range of values
			 
			 @ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType>
		static DoubleReal sum(IteratorType begin, IteratorType end)
		{
			return std::accumulate(begin, end, 0.0);
		}

		/**
			 @brief Calculates the mean of a range of values

			 @exception Exception::InvalidRange is thrown if the range is empty
			 
			 @ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType>
		static DoubleReal mean(IteratorType begin, IteratorType end)
		{
			SignedSize size = std::distance(begin, end);
			if (size <= 0)
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
			return sum(begin, end) / size;
		}
		
		/**
			 @brief Calculates the median of a range of values

			 @param sorted Is the range already sorted? (If not, it may be reordered.)

			 @exception Exception::InvalidRange is thrown if the range is empty
			 
			 @ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType>
		static DoubleReal median(IteratorType begin, IteratorType end, 
														 bool sorted=FALSE)
		{
			Size size = std::distance(begin, end);
			if (size <= 0)
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
			if (size % 2 == 0)  // even size => average two middle values
			{
				IteratorType it1 = begin, it2;
				std::advance(it1, size / 2 - 1);
				if (!sorted)
				{
					std::nth_element(begin, it1, end); // partition: smaller, nth, greater
					it2 = it1++;
					it1 = std::min_element(it1, end); // smallest in the greater partition
				}
				return (*it1 + *it2) / 2.0;
			}
			else 
			{
				IteratorType it = begin;
				std::advance(it, (size - 1) / 2);
				if (!sorted) std::nth_element(begin, it, end);
				return *it;
			}
		}

		/**
			@brief Calculates the mean square error for the values in [begin_a, end_a) and [begin_b, end_b)

			Calculates the mean square error for the data given by the two iterator ranges.

   		@exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

			@ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType1, typename IteratorType2>
	  static DoubleReal meanSquareError(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
		{
			//no data or different lengths
			SignedSize dist = std::distance(begin_a, end_a);
			if (dist == 0 || dist != std::distance(begin_b, end_b))
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}

			DoubleReal error = 0;
			while (begin_a != end_a)
			{
				error += pow(*begin_a - *begin_b, 2);
				++begin_a;
				++begin_b;
			}

			return error / dist;
		}

		/**
			@brief Calculates the classification rate for the values in [begin_a,	end_a) and [begin_b, end_b)

			Calculates the classification rate for the data given by the two iterator ranges.

	    @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

			@ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType1, typename IteratorType2>
	  static DoubleReal classificationRate(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
		{
			//no data or different lengths
		  SignedSize dist = std::distance(begin_a, end_a);
			if (dist == 0 || dist != std::distance(begin_b, end_b))
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}

			DoubleReal correct = (DoubleReal) dist;
			while (begin_a != end_a)
			{
				if ((*begin_a < 0 && *begin_b >= 0) || (*begin_a >= 0 && *begin_b < 0))
				{
					--correct;
				}
				++begin_a;
				++begin_b;
			}

			return correct / dist;
		}

		/**
			@brief Calculates the Matthews correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

			Calculates the Matthews correlation coefficient for the data given by the two iterator ranges. The values in [begin_a, end_a) have to be the predicted labels and the values in [begin_b, end_b) have to be the real labels.

	    @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

			@ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType1, typename IteratorType2>
	  static DoubleReal matthewsCorrelationCoefficient(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
		{
			//no data or different lengths
			Int dist = std::distance(begin_a, end_a);
			if (dist == 0 || dist != std::distance(begin_b, end_b))
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}

			DoubleReal tp = 0;
			DoubleReal fp = 0;
			DoubleReal tn = 0;
			DoubleReal fn = 0;

			while (begin_a != end_a)
			{
				if (*begin_a < 0 && *begin_b >= 0)
				{
					++fn;
				}
				else if (*begin_a < 0 && *begin_b < 0)
				{
					++tn;
				}
				else if (*begin_a >= 0 && *begin_b >= 0)
				{
					++tp;
				}
				else if (*begin_a >= 0 && *begin_b < 0)
				{
					++fp;
				}

				++begin_a;
				++begin_b;
			}

			return ((tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)));
		}

		/**
			@brief Calculates the Pearson correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

			Calculates the linear correlation coefficient for the data given by the two iterator ranges.

			If one of the ranges contains only the same values 'nan' is returned.

	    @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

			@ingroup MathFunctionsStatistics
		*/
		template <typename IteratorType1, typename IteratorType2>
	  static DoubleReal pearsonCorrelationCoefficient(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
		{
			//no data or different lengths
		  SignedSize dist = std::distance(begin_a, end_a);
			if (dist == 0 || dist != std::distance(begin_b, end_b))
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}

			//calculate average
			DoubleReal avg_a = std::accumulate(begin_a, end_a, 0.0) / dist;
			DoubleReal avg_b = std::accumulate(begin_b, end_b, 0.0) / dist;

			DoubleReal numerator = 0;
			DoubleReal denominator_a = 0;
			DoubleReal denominator_b = 0;
			while (begin_a != end_a)
			{
				DoubleReal temp_a = *begin_a - avg_a;
				DoubleReal temp_b = *begin_b - avg_b;
				numerator += (temp_a * temp_b);
				denominator_a += (temp_a * temp_a);
				denominator_b += (temp_b * temp_b);
				++begin_a;
				++begin_b;
			}

			return numerator / sqrt(denominator_a * denominator_b);
		}

    /// Computes the rank of the sorted vector @p w
    template <typename Value>
    static void computeRank(std::vector<Value>& w)
    {
      Size i = 0; // main index
      Size z  = 0;	// "secondary" index
      Value rank = 0;
      Size n = (w.size() - 1);

      while (i < n)
      {
				// test for equality with tolerance:
        if (fabs(w[i+1] - w[i]) > 0.0000001 * fabs(w[i+1])) // no tie
        {
          w[i] = Value(i);
          ++i;
        }
        else // tie, replace by mean rank
        {
					// count number of ties
          for (z = i + 1; z <= n && fabs(w[z] - w[i]) <= 0.0000001 * fabs(w[z]) ; ++z);
					// compute mean rank of tie
          rank = 0.5 * (i + z - 1);
					// replace intensities by rank
          for (Size v = i; v <= z - 1 ; ++v) w[v] = rank;

          i = z;
        }
      }
      if (i == n) w[n] = Value(n);
    }

    /**
	    @brief calculates the rank correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

	    Calculates the rank correlation coefficient for the data given by the two iterator ranges.

	    If one of the ranges contains only the same values 'nan' is returned.

	    @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

	    @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static DoubleReal rankCorrelationCoefficient(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
			//no data or different lengths
      SignedSize dist = std::distance(begin_a, end_a);
			if (dist == 0 || dist != std::distance(begin_b, end_b))
			{
				throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}

      // store and sort intensities of model and data
      std::vector<DoubleReal> ranks_data;
      ranks_data.reserve(dist);
      std::vector<DoubleReal> ranks_model;
      ranks_model.reserve(dist);

      while (begin_a != end_a)
      {
        ranks_model.push_back(*begin_a);
        ranks_data.push_back(*begin_b);
        ++begin_a;
        ++begin_b;
      }

      // sort ranks
      std::sort(ranks_data.begin(), ranks_data.end());
      std::sort(ranks_model.begin(), ranks_model.end());

      // compute rank
      computeRank(ranks_data);
      computeRank(ranks_model);

      DoubleReal mu = DoubleReal(ranks_data.size() + 1) / 2.; // mean of ranks
      // Was the following, but I think the above is more correct ... (Clemens)
      // DoubleReal mu = (ranks_data.size() + 1) / 2;

      DoubleReal sum_model_data = 0;
      DoubleReal sqsum_data = 0;
      DoubleReal sqsum_model = 0;

      for (Int i = 0; i < dist; ++i)
      {
        sum_model_data += (ranks_data[i] - mu) * (ranks_model[i] - mu);
        sqsum_data += (ranks_data[i] - mu) * (ranks_data[i] - mu);
        sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
      }

		  // check for division by zero
      if (!sqsum_data || !sqsum_model ) return 0;

      return sum_model_data / (sqrt(sqsum_data) * sqrt(sqsum_model));
    }

	} // namespace Math
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H
