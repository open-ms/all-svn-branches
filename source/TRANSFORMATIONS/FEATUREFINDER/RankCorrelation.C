// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/RankCorrelation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	RankCorrelation::RankCorrelation()
		:	BaseQuality()
	{
		setName(getProductName());
		check_defaults_ = false;
	}

	RankCorrelation::~RankCorrelation()
	{
	}

	double RankCorrelation::evaluate(const IndexSet& set, const BaseModel<2>& model)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		
		unsigned int n = set.size(); 						
		gsl_vector* dv = gsl_vector_alloc(n);
		gsl_vector* mv = gsl_vector_alloc(n);

		unsigned int i=0;
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			gsl_vector_set(dv,i,traits_->getPeakIntensity(*it) );
			gsl_vector_set(mv,i,model.getIntensity(  traits_->getPeakPos(*it) ));
			++i;
		}
		// compute ranks 
    gsl_permutation * perm1 = gsl_permutation_alloc(n);
    gsl_permutation * perm2 = gsl_permutation_alloc(n);
    
		gsl_permutation * rank1 = gsl_permutation_alloc(n);
   	gsl_permutation * rank2 = gsl_permutation_alloc(n);

    gsl_sort_vector_index(perm1,dv);
    gsl_permutation_inverse(rank1,perm1);
    
		gsl_sort_vector_index(perm2,mv);
    gsl_permutation_inverse(rank2,perm2);

		double sum_model_data, sqsum_data, sqsum_model = 0;
		double mu = (n + 1) / 2;

		for (unsigned int i=0; i<n;++i)
		{
			//cout << "ranks data " << rank1->data[i] << endl;
			//cout << "ranks model " << rank2->data[i] << endl;
			//cout << "------------------------" << endl;
			sum_model_data += (rank1->data[i] - mu) * (rank2->data[i] - mu);
			sqsum_data += (rank1->data[i] - mu) * (rank1->data[i] - mu);
			sqsum_model += (rank2->data[i] - mu) * (rank2->data[i] - mu);
		}	
																																															
		// check for division by zero
		if ( ! sqsum_data || ! sqsum_model ) return 0;		
		
		//cout << "sum_model_data " << sum_model_data << endl;
		//cout <<  "sqsum_data " << sqsum_data << endl;
		//cout << "sqsum_model " << sqsum_model << endl;

		double rs = sum_model_data /  sqrt(sqsum_data * sqsum_model); 
		//cout << "Spearman correlation: " << fabs(rs) << endl;
			
		// free gsl data structures
		gsl_vector_free(dv);
		gsl_vector_free(mv);
		gsl_permutation_free(perm1);
		gsl_permutation_free(rank1);
		gsl_permutation_free(perm2);
		gsl_permutation_free(rank2);
							
		UInt df = set.size()-1;
		double t_stat = sqrt(df) * rs; 
		
		// t_stat follows Normal Gaussian distribution 
		pval_ = (1 - gsl_cdf_ugaussian_P(t_stat));	
				
		return ( fabs(rs));
	}
	
	double RankCorrelation::evaluate(const IndexSet& set, const BaseModel<1>& model, UInt dim)
	{
			if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	
			unsigned int n = set.size(); 						
			gsl_vector* dv = gsl_vector_alloc(n);
			gsl_vector* mv = gsl_vector_alloc(n);

			unsigned int i=0;
			for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
			{
				const CoordinateType coord = traits_->getPeakPos(*it)[dim];
				gsl_vector_set(dv,i,traits_->getPeakIntensity(*it) );
				gsl_vector_set(mv,i,model.getIntensity( coord ) );
				++i;
			}
			// compute ranks 
    	gsl_permutation * perm1 = gsl_permutation_alloc(n);
    	gsl_permutation * perm2 = gsl_permutation_alloc(n);
    
			gsl_permutation * rank1 = gsl_permutation_alloc(n);
   		gsl_permutation * rank2 = gsl_permutation_alloc(n);

    	gsl_sort_vector_index(perm1,dv);
    	gsl_permutation_inverse(rank1,perm1);
    
			gsl_sort_vector_index(perm2,mv);
    	gsl_permutation_inverse(rank2,perm2);

			double sum_model_data, sqsum_data, sqsum_model = 0;
			double mu = (n + 1) / 2;

			for (unsigned int i=0; i<n;++i)
			{
				//cout << "ranks data " << rank1->data[i] << endl;
				//cout << "ranks model " << rank2->data[i] << endl;
				//cout << "------------------------" << endl;
				sum_model_data += (rank1->data[i] - mu) * (rank2->data[i] - mu);
				sqsum_data += (rank1->data[i] - mu) * (rank1->data[i] - mu);
				sqsum_model += (rank2->data[i] - mu) * (rank2->data[i] - mu);
			}	
																																															
			// check for division by zero
			if ( ! sqsum_data || ! sqsum_model ) return 0;		
		
			//cout << "sum_model_data " << sum_model_data << endl;
			//cout <<  "sqsum_data " << sqsum_data << endl;
			//cout << "sqsum_model " << sqsum_model << endl;

			double rs = sum_model_data /  sqrt(sqsum_data * sqsum_model); 
			//cout << "Spearman correlation: " << fabs(rs) << endl;
			
			// free gsl data structures
			gsl_vector_free(dv);
			gsl_vector_free(mv);
			gsl_permutation_free(perm1);
			gsl_permutation_free(rank1);
			gsl_permutation_free(perm2);
			gsl_permutation_free(rank2);
							
			UInt df = set.size()-1;
			double t_stat = sqrt(df) * rs; 
		
			// t_stat follows Normal Gaussian distribution 
			pval_ = (1 - gsl_cdf_ugaussian_P(t_stat));	
				
			return (fabs(rs));
		}
		
		double RankCorrelation::evaluate(const Math::LinearInterpolation<double,double>& lint, const BaseModel<1>& isomodel)
		{
			if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
						
			unsigned int n = lint.getData().size(); 						
			gsl_vector* dv = gsl_vector_alloc(n);
			gsl_vector* mv = gsl_vector_alloc(n);
			
			double data_sum, model_sum = 0.0;
			for (UInt i=0;i<lint.getData().size();++i)
			{
				data_sum += lint.getData()[i];
				model_sum += isomodel.getIntensity(lint.index2key(i));
			}
			
 			for (UInt i=0;i<lint.getData().size();++i)
			{
				double x = lint.getData()[i] / data_sum;
				double y = isomodel.getIntensity(lint.index2key(i)) / model_sum;

				gsl_vector_set(dv,i,x);
				gsl_vector_set(mv,i,y);

				//cout << "r : "  << lin_int.getData()[i] / data_sum << endl;
				//cout << "r : " << isomodel.getIntensity(lin_int.index2key(i)) /  model_sum << endl;
				//cout << "-----------------------------------------------------------" << endl;
			}
			// compute ranks 
    	gsl_permutation * perm1 = gsl_permutation_alloc(n);
    	gsl_permutation * perm2 = gsl_permutation_alloc(n);
    
			gsl_permutation * rank1 = gsl_permutation_alloc(n);
   		gsl_permutation * rank2 = gsl_permutation_alloc(n);

    	gsl_sort_vector_index(perm1,dv);
    	gsl_permutation_inverse(rank1,perm1);
    
			gsl_sort_vector_index(perm2,mv);
    	gsl_permutation_inverse(rank2,perm2);

			double sum_model_data, sqsum_data, sqsum_model = 0;
			double mu = (n + 1) / 2;

			for (unsigned int i=0; i<n;++i)
			{
				//cout << "ranks data " << rank1->data[i] << endl;
				//cout << "ranks model " << rank2->data[i] << endl;
				//cout << "------------------------" << endl;
				sum_model_data += (rank1->data[i] - mu) * (rank2->data[i] - mu);
				sqsum_data += (rank1->data[i] - mu) * (rank1->data[i] - mu);
				sqsum_model += (rank2->data[i] - mu) * (rank2->data[i] - mu);
			}	
																																															
			// check for division by zero
			if ( ! sqsum_data || ! sqsum_model ) return 0;		
		
			//cout << "sum_model_data " << sum_model_data << endl;
			//cout <<  "sqsum_data " << sqsum_data << endl;
			//cout << "sqsum_model " << sqsum_model << endl;

			double rs = sum_model_data /  sqrt(sqsum_data * sqsum_model); 
			//cout << "Spearman correlation: " << fabs(rs) << endl;
			
			// free gsl data structures
			gsl_vector_free(dv);
			gsl_vector_free(mv);
			gsl_permutation_free(perm1);
			gsl_permutation_free(rank1);
			gsl_permutation_free(perm2);
			gsl_permutation_free(rank2);
						
			return fabs(rs);
		}
		
}
