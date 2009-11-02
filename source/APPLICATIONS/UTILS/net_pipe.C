#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/XCorrelation.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <boost/math/distributions/normal.hpp> // for binomial distribution

#include <iostream>
#include <fstream>

typedef int32_t LONG;
typedef float sup;
typedef LONG TT;

using namespace OpenMS;
using namespace std;

template<class T> int loadBinArray(char *filename, vector<vector<T> >& data)
{
  ifstream is(filename, ios::in | ios::binary);  //open file stream
  if(!is)
  {
    // An error occurred!
    is.clear();
    cout << "bgwaaak!" << endl; // todo throw real error
  }

  UInt c; // read header information of array size (m x n)
  is.read(reinterpret_cast<char*>(&c), sizeof(UInt) );
  const UInt m = c;
  is.read(reinterpret_cast<char*>(&c), sizeof(UInt) );
  const UInt n = c;

  //~ int length;
  //~ streampos pos_type;
  //~ pos_type = is.tellg();
  //~ // get length of file:
  //~ is.seekg(0, ios::end);
  //~ length = is.tellg();
  //~ is.seekg(pos_type);
  //~ cout << length/(sizeof(T)) <<endl;

  T *buff = new T[n]; //read the file according to the header into a buffer
  vector< vector<T> > inread(m, vector<T>(n));
  for(Size rows = 0; rows < m; ++rows)
  {
    is.read(reinterpret_cast<char*>(buff), n*sizeof(T) );
    if(is.good())
    {
      for(Size cols = 0; cols < n; ++cols)
      {
        inread[rows][cols] = buff[cols];
      }
    }
    else
    {
      cerr << "Data is incorrect\n"; return EXIT_FAILURE;
    }
  }

  //~ c= 666;
  //~ is.read(reinterpret_cast<char*>(&c), sizeof(UInt) ); cerr << c << "§" << is.good() << endl;
  //~ is.read(reinterpret_cast<char*>(&c), sizeof(UInt) ); cerr << c << "§" << is.good() << endl << true << endl ;

  delete[] buff; //return all resources to where they belong
  is.close();
  data = inread;
  return EXIT_SUCCESS;
}

int loadPklBin(char *filename, PeakMap& data)
{
  ifstream is(filename, ios::in | ios::binary);  //open file stream
  if(!is)
  {
    // An error occurred!
    is.clear();
    cout << "bgwaaak!" << endl; // todo throw real error
  }

  UInt c; // read header information of map size
  is.read(reinterpret_cast<char*>(&c), sizeof(int32_t) );
  const UInt num_specs = c;
  int16_t *num_buff = new int16_t[num_specs]; //read the spectra header according to the file header into a buffer
  is.read(reinterpret_cast<char*>(num_buff), num_specs*sizeof(int16_t) );
  if(!is.good())
  {
    cerr << "Data is incorrect\n"; return EXIT_FAILURE;
  }
  for(Size i = 0; i < num_specs; ++i)
  {
    float *peak_buff = new float[ (2*num_buff[i]+2) ]; //read the spectra according to the file and spectrum header into a buffer
    is.read(reinterpret_cast<char*>(peak_buff), (2*num_buff[i]+2)*sizeof(float) );
    if(is.good())
    {
      PeakSpectrum tmp_spec;
      Precursor tmp_prec;
      tmp_prec.setIsolationWindowLowerOffset(peak_buff[0]);
      tmp_prec.setIsolationWindowUpperOffset(peak_buff[0]);
      if(peak_buff[1]<4)
      {
        tmp_prec.setCharge(peak_buff[1]);
      }
      //~ else
      //~ {
        //~ tmp_prec.setCharge(2);
      //~ }
      tmp_spec.getPrecursors().push_back(tmp_prec);

      for(Size j = 2; j < (2*num_buff[i]+2); ++j)
      {
        Peak1D tmp_peak;
        tmp_peak.setMZ(peak_buff[j]);
        ++j;
        tmp_peak.setIntensity(peak_buff[j]);
        tmp_spec.push_back(tmp_peak);
      }
      data.push_back(tmp_spec);
    }
    else
    {
      cerr << "Data is incorrect\n"; return EXIT_FAILURE;
    }
    delete[] peak_buff; //return all resources to where they belong
  }

  delete[] num_buff; //return all resources to where they belong
  is.close();
  return EXIT_SUCCESS;
}

  /* old snippets for direct comparison to matlab-shi#!
	  //input for id
  //~ char* name = "stars_indices.bin";
  //~ vector< vector<TT> > inread;
  //~ int res = loadBinArray<TT>(name, inread);
  //~ if(!res)
  //~ {
    //~ //test output
    //~ cout << "inread.size(): " << inread.size() << endl;
    //~ for(Size rows = 0; rows < inread.size(); ++rows)
    //~ {
      //~ cout << "inread["<< rows <<"].size(): " << inread[rows].size() << endl;
      //~ for(Size cols = 0; cols < inread[rows].size(); ++cols)
      //~ {
         //~ // cout << inread[rows][cols] << "-";
      //~ }
    //~ }
    //~ cout << endl;
  //~ }

  //~ TextFile pep_annots("spectra/peptides_tagsearch.txt", true); // todo try-catch
  //~ for(TextFile::iterator it = pep_annots.begin()+1; it != pep_annots.end(); ++it)
  //~ {
    //~ vector<String> tmp;
    //~ it->split('\t',tmp);

    //~ Size index, protein_index;
    //~ DoubleReal peptide_tag_score, perc_peptide_score;

    //~ if((tmp.size() == 6))
    //~ {
      //~ index = tmp[0].toInt()-1;// todo try-catch
      //~ peptide_tag_score = tmp[2].toDouble();// todo try-catch
      //~ perc_peptide_score = tmp[3].toDouble();// todo try-catch
      //~ protein_index = tmp[4].toInt();// todo try-catch
      //~ Size true_index = inread[0][index]-1;
      //~ if(true_index < experiment.size())
      //~ {
        //~ PeptideHit actual_hit(perc_peptide_score,0,-1,tmp[1].trim());
        //~ actual_hit.setProteinAccessions(vector<String>(1,tmp[5]));
        //~ PeptideIdentification pid;
        //~ pid.setScoreType("perc_peptide_score");
        //~ pid.insertHit(actual_hit);
        //~ experiment[true_index].getPeptideIdentifications().push_back(pid);
      //~ }
      //~ else
      //~ {
        //~ cout << "no good file " << endl; break;
      //~ }
    //~ }
    //~ else
    //~ {
      //~ cout << "no good file " << endl; break;
    //~ }
  //~ }

  //test output
  //~ for(Peakexperiment::iterator it = experiment.begin(); it != experiment.end(); ++it)
  //~ {
    //~ for(vector<PeptideIdentification>::const_iterator it_pid = it->getPeptideIdentifications().begin(); it_pid != it->getPeptideIdentifications().end(); ++it_pid)
    //~ {
      //~ for(vector<PeptideHit>::const_iterator it_hit = it_pid->getHits().begin(); it_hit != it_pid->getHits().end(); ++it_hit)
      //~ {
        //~ cout << it_hit->getSequence() << " ";
      //~ }
      //~ cout << " - ";
    //~ }
    //~ cout << " \'\' " << endl;
  //~ }

  //~ vector<PeakSpectrum> blists(experiment.size());
  //~ vector<Real> ei(experiment.size()), tp(experiment.size());

      vector<bool> matches_B(10,false);
      matches_B[0]=true;matches_B[4]=true;matches_B[9]=true;
      vector<Size> matched_indices;
      for(Size s = 0; s < matches_B.size(); ++s)
      {
        if(matches_B[s] == true){ matched_indices.push_back(s); }
      }
      for(Size s = 0; s < matches_B.size(); ++s)
      {
        cout << matches_B[s];
      }
      cout << endl;

      vector<Size> matched_indices_(matched_indices.size()-1);
      for(Size s = 0; s < matched_indices_.size(); ++s)
      {
        matched_indices_[s] = matched_indices[s+1] - matched_indices[s];
      }

      for(Size s = 0; s < matched_indices_.size(); ++s)
      {
        cout << matched_indices_[s];
      }
      cout << endl;

      vector<Size> ip_count_B(6,0); //n-te zelle enthält #(abstand consekutiver matches == n)
      for(Size s = 0; s < matched_indices_.size(); ++s)
      {
        cout << matched_indices_[s] << endl;
        switch(matched_indices_[s])
        {
          case 0: cout << "err" << endl; break;
          case 1: ++ip_count_B[0]; break;
          case 2: ++ip_count_B[1]; break;
          case 3: ++ip_count_B[2]; break;
          case 4: ++ip_count_B[3]; break;
          case 5: ++ip_count_B[4]; break;
          default: ++ip_count_B[5]; break;
        }
      }

  //read a binarray int32_t
  char* name_1 = "kept_indices.bin";
  vector< vector<int32_t> > inread_1;
  int res_ = loadBinArray<int32_t>(name_1, inread_1);
  if(!res_)
  {
    cout << name_1 << endl;
  }

  //input test
  char* name_2 = "aligns/pairs_modpos.bin";
  vector< vector<sup> > inread_2;
  res_ = loadBinArray<sup>(name_2, inread_2);
  if(!res_)
  {
    cout << name_2 << endl;
    //test output
    cout << "inread.size(): " << inread_2.size() << endl;
    for(Size rows = 0; rows < inread_2.size(); ++rows)
    {
      cout << "inread["<< rows <<"].size(): " << inread_2[rows].size() << endl;
      for(Size cols = 0; cols < inread_2[rows].size(); ++cols)
      {
         cout << inread[rows][cols] << "-";
      }
    }
    cout << endl;
  }

  char* pairs_file = "aligns/pairs.pklbin";
  loadPklBin(pairs_file, pairs_experiment);
  */

Int main()
{

	DoubleReal max_pm;
	DoubleReal pm_tol;
	DoubleReal min_intensity_coverage_ratio;
	DoubleReal p_value;

	  //input
	MzMLFile mzml;
	PeakMap experiment; //MSExperiment<Peak1D> see StandardTypes
	mzml.load("test.mzML",experiment);

  //~ // maybe a mgf too
  //~ //load a mgf
  //~ String filename_pairs("pairs_pklbin.mgf");
  //~ FileHandler fh;
  //~ FileTypes::Type file_type = fh.getType(filename_pairs);
  //~ if (file_type!=FileTypes::MGF)
  //~ {
    //~ cout << "Open file error: Could not determine file type of '"<<filename_pairs<<"'!";
    //~ return EXIT_FAILURE;
  //~ }
  //~ fh.loadExperiment(filename_pairs, pairs_experiment, file_type);

	std::vector<Size> selected_specs; //where indices to not-to-bad-to-be-filtered-out-spectra's (or ms1) indices come; only spectra to those will be considered
	std::vector< DoubleReal > total_intensity; // total intensity of the selected spectra (for intensity-coverage ratio with matches)
	std::vector< std::pair<Size,Size> > pot_pairs; // name is program; will be reduced in run of the program

	for(Size i = 0; i < experiment.size(); ++i)
	{
		//select specs
			//filter criteria: has no valid front() precursor in getPrecursors(), is MS1, ...
      if( (!experiment[i].getPrecursors().empty()) and (experiment[i].getMSLevel()==2) )
      {
        selected_specs.push_back(i);
      }
	}

	//~ for(Size i = 0; i < experiment.size(); ++i)
	//~ {
    //modify selected specs
			// remove isotope peaks, remove low peaks, remove peaks > .getPrecursors().front().getMZ(),
	//~ }

  total_intensity.reserve(selected_specs.size());
	for(Size i = 0; i < selected_specs.size()-1; ++i)
	{
    //fill pot_pairs
		for(Size j = i+1; j < selected_specs.size(); ++j)
		{
			if( !(abs(experiment[j].getPrecursors().front().getMZ()-experiment[i].getPrecursors().front().getMZ()) > max_pm+pm_tol) )
			{
				pot_pairs.push_back(std::pair<Size,Size>(i,j));
			}
		}
    //fill total_intensity
    DoubleReal intens;
    for(Size j = 0; j < experiment[i].size(); ++j)
    {
      intens += experiment[i][j].getIntensity();
    }
	}

	std::vector< std::pair<DoubleReal,DoubleReal> > pot_pairs_xcs; // the pairs xcorrelation scores
	std::vector< Size > stats_num(selected_specs.size()); // the number of scores for respective spec (for online calc of the var)
	std::vector< DoubleReal > means(selected_specs.size()), variances(selected_specs.size()), stddevs(selected_specs.size()); // the selected specs xcorrelation means, variation and std. deviation (unbiased gaussian)

  std::cout << pot_pairs.size() << std::endl;

  pot_pairs_xcs.reserve(pot_pairs.size());
	XCorrelation<Peak1D> x_corr;
	for(Size i = 0; i < pot_pairs.size(); ++i)
	{
		DoubleReal best_score1, best_score2, best_shift;
		std::list<std::pair<Size,Size> > best_matches;
		//~ x_corr.getXCorrelation(experiment[pot_pairs[i].first], experiment[pot_pairs[i].second], best_score1, best_score2, best_shift, best_matches);
		//calc mean and variation

			/* means[spec1][0]=(means[spec1][0]*means[spec1][1]+curResult.score1)/(means[spec1][1]+1);
			means[spec2][0]=(means[spec2][0]*means[spec2][1]+curResult.score2)/(means[spec2][1]+1);
			means[spec1][1]++;    means[spec2][1]++;
			float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
			varTerms[spec1] = updateRatio*varTerms[spec1]+(curResult.score1*curResult.score1)/means[spec1][1];
			updateRatio = (means[spec2][1]-1)/means[spec2][1];
			varTerms[spec2] = updateRatio*varTerms[spec2]+(curResult.score2*curResult.score2)/means[spec2][1]; */

    means[pot_pairs[i].first] = ( means[pot_pairs[i].first]*stats_num[pot_pairs[i].first] + best_score1 )/ (1+stats_num[pot_pairs[i].first]);
    means[pot_pairs[i].second] = ( means[pot_pairs[i].second]*stats_num[pot_pairs[i].second] + best_score2 )/ (1+stats_num[pot_pairs[i].second]);
    variances[pot_pairs[i].first] = ( variances[pot_pairs[i].first]*stats_num[pot_pairs[i].first] + (best_score1*best_score1) )/ (1+stats_num[pot_pairs[i].first]);
    variances[pot_pairs[i].second] = ( variances[pot_pairs[i].second]*stats_num[pot_pairs[i].second] + (best_score1*best_score2) )/ (1+stats_num[pot_pairs[i].second]);
    ++stats_num[pot_pairs[i].first];
    ++stats_num[pot_pairs[i].second];
	}

	//calc std.dev. sigma as sqrt of variance
	for(Size i = 0; i < selected_specs.size(); ++i)
	{
    stddevs[i]= sqrt( (means[i]/(means[i]-1)) * (variances[i]-means[i]*means[i]) );
  }

	//filter pairs not fitting ratio or gcdf
    std::vector< std::pair<Size,Size> >::iterator it_pairs = pot_pairs.begin();
    std::vector< std::pair<DoubleReal,DoubleReal> >::iterator it_pairs_xcs = pot_pairs_xcs.begin();
		while(it_pairs != pot_pairs.end())
    {
			if( boost::math::cdf(boost::math::normal(means[it_pairs->first],stddevs[it_pairs->first]), it_pairs_xcs->first) >= 1-p_value and
          boost::math::cdf(boost::math::normal(means[it_pairs->second],stddevs[it_pairs->second]), it_pairs_xcs->second) >= 1-p_value and
          it_pairs_xcs->first/total_intensity[it_pairs->first] >= min_intensity_coverage_ratio and
          it_pairs_xcs->second/total_intensity[it_pairs->second] >= min_intensity_coverage_ratio)
			{
        ++it_pairs;
        ++it_pairs_xcs;
      }
			else
      {
        it_pairs = pot_pairs.erase(it_pairs);
        it_pairs_xcs = pot_pairs_xcs.erase(it_pairs_xcs);
      }
		}

    std::cout << pot_pairs.size() << std::endl;


  return EXIT_SUCCESS;
} //end of main
