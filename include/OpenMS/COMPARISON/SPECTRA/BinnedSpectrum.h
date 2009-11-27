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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/SparseVector.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <boost/math/distributions/binomial.hpp> // for binomial distribution

#include <cmath>

namespace OpenMS
{

	class Peak1D;
	/**
	@brief This is a binned representation of a PeakSpectrum

		@param sz the size of the bins and
		@param sp number of neighboring bins to both sides affected by a peak contribution
		@param ps the peakspectrum, that shall be represented

		sz denotes the size of a bin in @p Th, thereby deciding the number of bins(all of size sz) the spectrum is discretized to.
		Each bin will represent a certain @p Th range and the peaks will be put in the respective bins and sum up inside.
		sp denotes the number of neighboring bins to the left and the number of neighboring bins to the right a peak is also added to.
		E.g. a BinnedSpectrum with binsize of 0.5 @p Th will have a peak at 100 @p Th in bin no. 200, a peak at 100.1 @p Th will be in bin no. 201.
		If the binspread is 1, the peak at 100 Th will be added to bin no. 199, 200 and 201.
		If the binspread is 2, the peak at 100 @p Th will also be added to bin no. 198 and 202, and so on.

		@ingroup SpectraComparison
	*/
	template <typename PeakT = Peak1D>
	class BinnedSpectrum
	: /* public DefaultParamHandler, */
		public MSSpectrum<PeakT>
	{

	private:

	UInt bin_spread_;
	Real bin_size_;
	SparseVector<Real> bins_;

	public:

	/**
		@brief 	Exception which is thrown if BinnedSpectrum bins are accessed and no PeakSpektrum has been
				integrated yet i.e. bins_ is empty
	*/
	class OPENMS_DLLAPI NoSpectrumIntegrated
		: public Exception::BaseException
	{
	public:
		NoSpectrumIntegrated(const char* file, int line, const char* function, const char* message ="BinnedSpectrum hasn't got a PeakSpectrum to base on yet") throw()
          : BaseException(file, line, function, "BinnedSpectrum::NoSpectrumIntegrated", message)
	{
	}

		~NoSpectrumIntegrated() throw()
	{
	}

	};

	typedef SparseVector<Real>::const_iterator const_bin_iterator;
	typedef SparseVector<Real>::iterator bin_iterator;

	/// default constructor
	BinnedSpectrum()
	   : MSSpectrum<PeakT>(), /* DefaultParamHandler("BinnedSpectrum"), */ bin_spread_(1), bin_size_(2.0), bins_()
	{
	}

	/// detailed constructor
	BinnedSpectrum(Real size, UInt spread, const PeakSpectrum& ps)
	   : MSSpectrum<PeakT>(ps), /* DefaultParamHandler("BinnedSpectrum"), */ bin_spread_(spread), bin_size_(size), bins_()
	{
		setBinning();
	}

	/// detailed constructor
	BinnedSpectrum(Real size, UInt spread, std::vector< MSSpectrum<PeakT> >& unmerged)
	   : MSSpectrum<PeakT>(), /* DefaultParamHandler("BinnedSpectrum"), */ bin_spread_(spread), bin_size_(size), bins_()
	{
		//~  find largest mz in given spectra
		DoubleReal max_size(max_element(unmerged.begin(),unmerged.end(),typename MSSpectrum<PeakT>::PMLess())->getPrecursors().front().getMZ());
		//~ set sparsevector size accordingly
		bins_ = SparseVector<Real>((UInt)ceil(max_size/bin_size_) + bin_spread_ ,0,0); // aka overall intensity of peaks in corresponding bins
		std::vector< SparseVector<Real> > binned(unmerged.size(),bins_); // aka each spectrums binned version
		SparseVector<int> synthetic_peaks((UInt)ceil(max_size/bin_size_) + bin_spread_ ,0,0); // additional for number of sythetic peaks in each bin
		SparseVector<int> overall_number((UInt)ceil(max_size/bin_size_) + bin_spread_ ,0,0); // overall number of peaks in corresponding bins
		SparseVector<Real> overall_center((UInt)ceil(max_size/bin_size_) + bin_spread_ ,0,0); // overall m/z center of peaks in corresponding bins

		if(unmerged.empty())
		{
			throw Exception::IllegalArgument (__FILE__, __LINE__, __PRETTY_FUNCTION__, "no spectra to merge given");
		}

		//~ will point to the integerdataarray named "synthetic peaks"
		Size ida_spa = 0;
		for(; ida_spa < unmerged.front().getIntegerDataArrays().size(); ++ida_spa)
		{
			if(unmerged.front().getIntegerDataArrays()[ida_spa].getName()=="synthetic peaks")
			{
				break;
			}
		}

		for(Size i = 0; i < unmerged.size(); ++i)
		{
			//put all peaks into bins
			UInt bin_number;
			for (Size j = 0; j < unmerged[i].size(); ++j)
			{
				if(unmerged[i][j].getMZ() > max_size)
				{
					break;
				}

				//bin_number counted form 0 -> floor
				bin_number = (UInt)floor(unmerged[i][j].getMZ()/bin_size_);
				//e.g. bin_size_ = 1.5: first bin covers range [0,1.5] so peak at 1.5 falls in first bin (index 0)
				if(unmerged[i][j].getMZ()/bin_size_ == (DoubleReal)bin_number and unmerged[i][j].getMZ() != 0)
				{
					--bin_number;
				}

				// add peak to corresponding bin
				/*debug std::cout << "bin_number " <<  bin_number << std::endl;*/
				binned[i][bin_number] = binned[i].at(bin_number) + unmerged[i][j].getIntensity();

				overall_number[bin_number] = overall_number.at(bin_number) + 1;
				overall_center[bin_number] = overall_center.at(bin_number) + unmerged[i][j].getMZ();

				if(unmerged[i].getIntegerDataArrays().size()>ida_spa and unmerged[i].getIntegerDataArrays()[ida_spa].getName()=="synthetic peaks")
				{
					synthetic_peaks[bin_number] = synthetic_peaks.at(bin_number) + unmerged[i].getIntegerDataArrays()[ida_spa][j];
				}

				// add peak to neighboring binspread many
				for (Size k = 0; k < bin_spread_; ++k)
				{
					binned[i][bin_number+k+1] = binned[i].at(bin_number+k+1) + unmerged[i][j].getIntensity();
					overall_number[bin_number] = overall_number.at(bin_number) + 1;
					overall_center[bin_number] = overall_center.at(bin_number) + unmerged[i][j].getMZ();
					if(unmerged[i].getIntegerDataArrays().size()>ida_spa and unmerged[i].getIntegerDataArrays()[ida_spa].getName()=="synthetic peaks")
					{
						synthetic_peaks[bin_number] = synthetic_peaks.at(bin_number) + unmerged[i].getIntegerDataArrays()[ida_spa][j];
					}
					// we are not in one of the first bins (0 to bin_spread)
					// not working:  if (bin_number-k-1 >= 0)
					if (bin_number >= k+1)
					{
						binned[i][bin_number-k-1] = binned[i].at(bin_number-k-1) + unmerged[i][j].getIntensity();
						overall_number[bin_number] = overall_number.at(bin_number) + 1;
						overall_center[bin_number] = overall_center.at(bin_number) + unmerged[i][j].getMZ();
						if(unmerged[i].getIntegerDataArrays().size()>ida_spa and unmerged[i].getIntegerDataArrays()[ida_spa].getName()=="synthetic peaks")
						{
							synthetic_peaks[bin_number] = synthetic_peaks.at(bin_number) + unmerged[i].getIntegerDataArrays()[ida_spa][j];
						}
					}
				}
			}
		}

		//~ merge down all bins according to given statistical hypothesis
		Real noise_propability(0.001); // noise prob. for bernoulli trial nois peak generation
		//~ trial number n, to succeed k times in a bernoulli trial follows a "negativen Binomialverteilung" but number of r successful trials out of n follows a normal "Binomialverteilung"
		Real propability(0.99); // prob. seeing k in n trials by chance is 0.01
		int trials(unmerged.size());
		//~afaik not right like that for boost: Real success(1-(pow(1-noise_propability,(1+2*spread))));
		Real success(pow(1-noise_propability,(1+2*spread)));
			//~ from boost docu: The smallest number of successes that may be observed from n (or 'trials') trials with success fraction p (or ('success'), at probability P (or 'propability'). Note that the value returned is a real-number, and not an integer. Depending on the use case you may want to take either the floor or ceiling of the result. For example:
		int min_peak_concurrence = floor(boost::math::quantile(boost::math::complement(boost::math::binomial(trials, success), propability))); /// @improvement find out whether to ceil or to floor;
		/// @improvement run over the SparseVector with Iterators to speed up
		for(Size n = 0; n < overall_number.size(); ++n)
		{
			if(overall_number.at(n) >= (Real)min_peak_concurrence)
			{
				for(Size y = 0; y < binned.size(); ++y)
				{
					bins_[n] = bins_.at(n) + binned[y].at(n);
				}
			}
		}
		for(Size n = 0; n < overall_number.size(); ++n)
		{
			if((int)overall_number.at(n) >= min_peak_concurrence)
			{
				bins_[n] = bins_.at(n) / (Real)overall_number.at(n);
				overall_center[n] = overall_center.at(n) / (Real)overall_number.at(n);
			}
		}

		std::vector<Real> synthetics;
		//~ from merged bins fill peaks in integrated MSSpectrum
		for(Size n = 0; n < bins_.size(); ++n)
		{
			if((int)overall_number.at(n) >= min_peak_concurrence and (int)overall_number.at(n) > 0)
			{
				PeakT tmp;
				tmp.setIntensity(bins_[n]);
				tmp.setMZ(overall_center[n]);
				this->push_back(tmp);
				synthetics.push_back((Real)synthetic_peaks.at(n)/(Real)overall_number.at(n));
			}
		}
		this->setRT(unmerged.front().getRT());
		this->setPrecursors(unmerged.front().getPrecursors());
		this->setMSLevel(unmerged.front().getMSLevel());
		this->getFloatDataArrays().resize(1); // imperative that we have one
		std::swap(this->getFloatDataArrays().front(), synthetics); // imperative that inserted befor sorting
		this->getFloatDataArrays().front().setName("synthetic consensus peaks");
		this->sortByPosition();
	}

	/// copy constructor
	BinnedSpectrum(const BinnedSpectrum& source)
		: MSSpectrum<PeakT>(source), /* DefaultParamHandler(source), */ bin_spread_(source.getBinSpread()), bin_size_(source.getBinSize()), bins_(source.getBins())
	{
	}

	/// destructor
	~BinnedSpectrum()
	{
	}

	/// assignment operator
	BinnedSpectrum& operator= (const BinnedSpectrum& source)
	{
		if (&source != this)
		{
			setBinSize(source.getBinSize());
			setBinSpread(source.getBinSpread());
			bins_ = source.getBins();
			MSSpectrum<PeakT>::operator=(source);
		}
		return *this;
	}

	/// assignment operator for PeakSpectra
	BinnedSpectrum& operator= (const PeakSpectrum& source)
	{
		if (!MSSpectrum<PeakT>::operator==(source))
		{
			MSSpectrum<PeakT>::operator=(source);
			setBinning();
		}
		return *this;
	}

	/// equality operator
	bool operator== (const BinnedSpectrum& rhs) const
	{
		return
			(MSSpectrum<PeakT>::operator==(rhs) &&
			rhs.getBinSize()==this->bin_size_ &&
			rhs.getBinSpread()==this->bin_spread_)
			;
	}

	/// inequality operator
	bool operator!= (const BinnedSpectrum& rhs) const
	{
		return !(operator==(rhs));
	}

	/// equality operator for PeakSpectra
	bool operator== (const PeakSpectrum& rhs) const
	{
		return MSSpectrum<PeakT>::operator==(rhs);
	}

	/// inequality operator for PeakSpectra
	bool operator!= (const PeakSpectrum& rhs) const
	{
		return !(operator==(rhs));
	}

	/// get the BinSize
	inline double getBinSize() const
	{
		return this->bin_size_;
	}

	/// get the BinSpread
	inline UInt getBinSpread() const
	{
		return this->bin_spread_;
	}

	/// get the BinNumber, number of Bins
	inline UInt getBinNumber() const
	{
		return (UInt)this->bins_.size();
	}

	/// get the FilledBinNumber, number of filled Bins
	inline UInt getFilledBinNumber() const
	{
		return (UInt)this->bins_.nonzero_size();
	}

	/** unmutable access to the Bincontainer

			@throw NoSpectrumIntegrated is thrown if no spectrum was integrated
	*/
	inline const SparseVector<Real>& getBins() const
	{
		if(bins_.size() == 0)
		{
			throw NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return bins_;
	}

	/** mutable access to the Bincontainer

			@throw NoSpectrumIntegrated is thrown if no spectrum was integrated
	*/
	inline SparseVector<Real>& getBins()
	{
		if(bins_.size() == 0)
		{
			try
			{
				this->setBinning();
			}
			catch(...)
			{
				throw NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}
		return bins_;
	}

	/// returns the const begin iterator of the container
	inline const_bin_iterator begin() const
	{
		return bins_.begin();
	}
	/// returns the const end iterator of the container
	inline const_bin_iterator end() const
	{
		return bins_.end();
	}

	/// returns the begin iterator of the container
	inline bin_iterator begin()
	{
		return bins_.begin();
	}
	/// returns the end iterator of the container
	inline bin_iterator end()
	{
		return bins_.end();
	}

	/** sets the BinSize_ (and rebinnes)

			@param s defines the size of the bins
			@throw NoSpectrumIntegrated is thrown if no spectrum is integrated
	*/
	inline void setBinSize(double s)
	{
		if(this->bin_size_ != s)
		{
			this->bin_size_ = s;
			try
			{
				this->setBinning();
			}
			catch(...)
			{
				throw NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}
	}


	/** sets the BinSpread_ (and rebinnes)

			@param s defines the binning spread, given as positive integer
			@throw NoSpectrumIntegrated is thrown if no spec was integrated into the instance
	*/
	inline void setBinSpread(UInt s)
	{
		if(this->bin_spread_ != s)
		{
			this->bin_spread_ = s;
			try
			{
				this->setBinning();
			}
			catch(...)
			{
				throw NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}
	}


	/** makes the binning: all Peaks of the containing PeakSpectrum are summed up in the bins corresponding to @p m/z ranges

			@throw NoSpectrumIntegrated is thrown if no spectrum was integrated before
	*/
	void setBinning()
	{
	 	if (this->MSSpectrum<PeakT>::empty())
	 	{
		 	throw NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		bins_.clear();

	 	//make all necessary bins accessible
		this->sortByPosition();
	 	bins_ = SparseVector<Real>((UInt)ceil(this->back().getMZ()/bin_size_) + bin_spread_ ,0,0);
	 	//~ bins_ = SparseVector<Real>((UInt)ceil(this->getPrecursors().front().getMZ()/bin_size_) + bin_spread_ ,0,0);

		//put all peaks into bins
		UInt bin_number;
		for (Size i = 0; i < this->size(); ++i)
		{
			//bin_number counted form 0 -> floor
			bin_number = (UInt)floor(this->operator[](i).getMZ()/bin_size_);
			//e.g. bin_size_ = 1.5: first bin covers range [0,1.5] so peak at 1.5 falls in first bin (index 0)

			if(this->operator[](i).getMZ()/bin_size_ == (DoubleReal)bin_number and this->operator[](i).getMZ() != 0)
			{
				--bin_number;
			}

			//add peak to corresponding bin
			bins_[bin_number] = bins_.at(bin_number) + this->operator[](i).getIntensity();

			//add peak to neighboring binspread many
			for (Size j = 0; j < bin_spread_; ++j)
			{
				bins_[bin_number+j+1] = bins_.at(bin_number+j+1) + this->operator[](i).getIntensity();
				// we are not in one of the first bins (0 to bin_spread)
				//not working:  if (bin_number-j-1 >= 0)
				if (bin_number >= j+1)
				{
					bins_[bin_number-j-1] = bins_.at(bin_number-j-1) + this->operator[](i).getIntensity();
				}
			}
		}

	}

	/// function to check comparability of two BinnedSpectrum objects, i.e. if they have equal bin size and spread
		//yields false if given BinnedSpectrum size or spread differs from this one (comparing those might crash)
	bool checkCompliance(const BinnedSpectrum<PeakT>& bs) const
	{
		return
			(this->bin_size_ == bs.getBinSize()) &&
			(this->bin_spread_ == bs.getBinSpread());
	}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
