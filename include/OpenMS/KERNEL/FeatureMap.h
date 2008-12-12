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
// $Maintainer: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_FEATUREMAP_H
#define OPENMS_KERNEL_FEATUREMAP_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <algorithm>
#include <vector>

namespace OpenMS
{

	/**	
		@brief A container for features.
		
		A map is a container holding 2-dimensional features,
		which in turn represent chemical entities (peptides, proteins, etc.) found
		in a 2-dimensional experiment.
		
		Maps are implemented as vectors of features and have basically the same interface
		as an STL vector has (model of Random Access Container and Back Insertion Sequence).
		
		Feature maps are typically created from peak data of 2D runs through the FeatureFinder.
		
		@improvement Add list of unassigned peptide features; allow loading and storing; change IDMapper; add to TextExport (Hiwi)
		
		@ingroup Kernel
	*/
	template <typename FeatureT = Feature >
	class OPENMS_DLLAPI FeatureMap
		: public std::vector<FeatureT>,
			public RangeManager<2>,
			public DocumentIdentifier
	{
		public:
			/**	
				 @name Type definitions
			*/
			//@{
			typedef FeatureT FeatureType;
			typedef RangeManager<2> RangeManagerType;
			typedef std::vector<FeatureType> Base;
			typedef typename Base::iterator Iterator;
			typedef typename Base::const_iterator ConstIterator;
			typedef typename Base::reverse_iterator ReverseIterator;
			typedef typename Base::const_reverse_iterator ConstReverseIterator;
			typedef FeatureType& Reference;
			typedef const FeatureType& ConstReference;
			//@}
			/**	
				 @name Constructors and Destructor
			*/
			//@{
			
			/// Default constructor
			FeatureMap()
				: Base(),
					RangeManagerType(),
					DocumentIdentifier(),
					protein_identifications_(),
					data_processing_()
			{
			}
			
			/// Copy constructor
			FeatureMap(const FeatureMap& source) 
				: Base(source),
					RangeManagerType(source),
					DocumentIdentifier(source),
					protein_identifications_(source.protein_identifications_),
					data_processing_(source.data_processing_)
			{
			}
			
			/// Destructor
			virtual ~FeatureMap()
			{
			}
			
			//@}
				
			/// Assignment operator
			FeatureMap& operator = (const FeatureMap& rhs)
			{
				if (&rhs==this) return *this;
					
				Base::operator=(rhs);
				RangeManagerType::operator=(rhs);
				DocumentIdentifier::operator=(rhs);
				protein_identifications_ = rhs.protein_identifications_;
				data_processing_ = rhs.data_processing_;

				return *this;
			}
	
			/// Equality operator
			bool operator == (const FeatureMap& rhs) const
			{
				return
					std::operator==(*this, rhs) &&
					RangeManagerType::operator==(rhs) &&
					DocumentIdentifier::operator==(rhs) &&
					protein_identifications_==rhs.protein_identifications_ &&
					data_processing_ == rhs.data_processing_
					;
			}
				
			/// Equality operator
			bool operator != (const FeatureMap& rhs) const
			{
				return !(operator==(rhs));
			}

			/**	
				@name Sorting.
				These simplified sorting methods are supported in addition to	
				the standard sorting methods of std::vector.
			*/
			//@{
			/// Sorts the peaks according to ascending intensity.
			void sortByIntensity(bool reverse=false)
			{ 
				if (reverse)
				{
					std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::IntensityLess()) );
				}
				else
				{
					std::sort(this->begin(), this->end(), typename FeatureType::IntensityLess() ); 
				}
			}
				
			///Sort features by position. Lexicographical comparison (first RT then m/z) is done.
			void sortByPosition() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::PositionLess() );
			}
			
			///Sort features by RT position.
			void sortByRT() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::RTLess() );
			}

			///Sort features by m/z position.
			void sortByMZ() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::MZLess() );
			}
			
			///Sort features by ascending overall quality.
			void sortByOverallQuality(bool reverse=false) 
			{
				if (reverse)
				{
					std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::OverallQualityLess()) );
				}
				else
				{
					std::sort(this->begin(), this->end(), typename FeatureType::OverallQualityLess() );
				}
			}
			//@}
			
			// Docu in base class
			void updateRanges()
			{
				this->clearRanges();
				updateRanges_(this->begin(),this->end());
				
				//enlarge the range by the convex hull points
				for (UInt i=0; i<this->size(); ++i)
				{
					DBoundingBox<2> box = this->operator[](i).getConvexHull().getBoundingBox();
					if (!box.isEmpty())
					{
						//update RT
						if (box.min()[Peak2D::RT] < this->pos_range_.min()[Peak2D::RT])
						{
							this->pos_range_.setMinX(box.min()[Peak2D::RT]);
						}
						if (box.max()[Peak2D::RT] > this->pos_range_.max()[Peak2D::RT])
						{
							this->pos_range_.setMaxX(box.max()[Peak2D::RT]);
						}
						//update m/z
						if (box.min()[Peak2D::MZ] < this->pos_range_.min()[Peak2D::MZ])
						{
							this->pos_range_.setMinY(box.min()[Peak2D::MZ]);
						}
						if (box.max()[Peak2D::MZ] > this->pos_range_.max()[Peak2D::MZ])
						{
							this->pos_range_.setMaxY(box.max()[Peak2D::MZ]);
						}
					}
				}
			}

			/// Swaps the content of this map with the content of @p from
			void swap(FeatureMap& from)
			{
				FeatureMap tmp;

				//range information
				tmp.RangeManagerType::operator=(*this);
				this->RangeManagerType::operator=(from);
				from.RangeManagerType::operator=(tmp);

				//swap actual features
				Base::swap(from);
				
				// swap DocumentIdentifier
				DocumentIdentifier::swap(from);
				
				// swap the remaining members
				protein_identifications_.swap(from.protein_identifications_);
				data_processing_.swap(from.data_processing_);
			}
			
			/// non-mutable access to the protein identifications
		 	const std::vector<ProteinIdentification>& getProteinIdentifications() const
		 	{
		  	return protein_identifications_;	   		
		 	}	
		 		    	
			/// mutable access to the protein identifications
		  std::vector<ProteinIdentification>& getProteinIdentifications()
		  {
		  	return protein_identifications_;	
		  }

			/// sets the protein identifications
		  void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications)
		  {
		  	protein_identifications_ = protein_identifications;
		  }
		  
			/// adds a protein identifications
		  void addProteinIdentification(ProteinIdentification& protein_identification)
		  {
		  	protein_identifications_.push_back(protein_identification);
		  }

			/// returns a const reference to the description of the applied data processing 
			const std::vector<DataProcessing>& getDataProcessing() const
			{
				return data_processing_; 
			}

			/// returns a mutable reference to the description of the applied data processing 
			std::vector<DataProcessing>& getDataProcessing()
			{
				return data_processing_; 
			}
			
			/// sets the description of the applied data processing 
			void setDataProcessing(const std::vector<DataProcessing>& processing_method)
			{
				data_processing_ = processing_method; 
			}

		protected:
			
			/// protein identifications
			std::vector<ProteinIdentification> protein_identifications_;
			
			/// applied data processing
			std::vector<DataProcessing> data_processing_;
			
	};
	
	/// Print content of a feature map to a stream.
	template <typename FeatureType >
	std::ostream& operator << (std::ostream& os, const FeatureMap<FeatureType>& map)
	{
		os << "# -- DFEATUREMAP BEGIN --"<< std::endl;
		os << "# POSITION \tINTENSITY\tOVERALLQUALITY\tCHARGE" << std::endl; 
		for (typename FeatureMap<FeatureType>::const_iterator iter = map.begin(); iter!=map.end(); iter++)
		{
			os << iter->getPosition() << '\t'
				 << iter->getIntensity() << '\t'
				 << iter->getOverallQuality() << '\t'
				 << iter->getCharge()
				 << std::endl;
		}
		os << "# -- DFEATUREMAP END --"<< std::endl;
		return os;
	}
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATUREMAP_H
