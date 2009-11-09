#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/XCorrelation.h>
#include <OpenMS/COMPARISON/SPECTRA/AntisymetricAlignment.h>
#include <OpenMS/COMPARISON/CLUSTERING/StarClusters.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <boost/math/distributions/normal.hpp> // for binomial distribution
#include <iostream>
#include <cmath>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@brief topp tool test in workbench

	<B>The command line parameters of this tool are:</B>
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

	class SpectralNetworking
	: public TOPPBase
	{
	 public:
		SpectralNetworking()
			: TOPPBase("SpectralNetworking","pipeline for spectral networking", false, "0.3beta")
		{

		}

 protected:

	virtual void registerOptionsAndFlags_()
	{
		//~ file section
		registerInputFile_("in","<file>","","Input ms/ms data in mzML format",true);
		setValidFormats_("in",StringList::create("mzML"));
		registerInputFile_("ids","<file>","","Input of ids to consensuses for propagation",false);
		setValidFormats_("ids",StringList::create("idXML"));
		registerInputFile_("network","<file>","","Input of network edges for propagation",false);
		registerStringOption_("directory", "output directory", "", "output directory", false);
		registerStringOption_("prefix", "prefix for the outputfiles", "", "prefix for the outputfiles", true);

		//~ internal parameters section
		registerDoubleOption_("peak_tolerance", "<tolerance>", 0.3, "Defines the absolut (in Da) peak tolerance", false);
		registerDoubleOption_("parentmass_tolerance", "<tolerance>", 3.0, "Defines the absolut (in Da) peak tolerance", false);
		registerDoubleOption_("min_dist", "<tolerance>", 57.0214637230, "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching", false);
		registerDoubleOption_("scan_resolution", "<tolerance>", 0.1, "Defines the resolution of the instrument used to create the spectra", false);
		registerDoubleOption_("max_pm_diff", "<limit>", 150.0, "Defines the maximal absolut (in Da) shift between the two spectra", false);
		registerDoubleOption_("consensus_sz", "<limit>", 0.3, "Defines the bin size (in Da) while makeing a consensus from two or more spectra", false);
		registerDoubleOption_("pair_filter_p_value", "<limit>", 0.05, "Defines the p_value while filtering optential spectral pairs", false);
		registerDoubleOption_("pair_filter_min_ratio", "<limit>", 0.4, "Defines the minimal matchratio while filtering optential spectral pairs", false);
		registerIntOption_("propagation_max_hops", "<limit>", 3, "Defines the maximal number of hops a identification may be propagated (saves time and space)", false);

		//~ addEmptyLine_();
		//~ addText_("Parameters for the filter can only be fiven in the INI file.");

	}
	///

	//~ builds a spectral network ready for initial identification and subsequent propagation
	void spanNetwork_(MSExperiment<Peak1D>& experiment, Param& param, String& outputfile_name_consensus, String& outputfile_name_edges)
	{
		//~ parameters
		DoubleReal pm_tol = param.getValue("peak_tolerance");
		DoubleReal max_pm_diff = param.getValue("max_pm_diff_diff");
		DoubleReal pairfilter_p_value = param.getValue(" pairfilter_p_value");
		DoubleReal pairfilter_min_ratio = param.getValue("pairfilter_min_ratio");
		DoubleReal scan_resolution = param.getValue("scan_resolution");
		Real consensus_sz = param.getValue("consensus_sz"); /* e.g. peak_tolerance */;
		/// @improvement is floor or round more appropriate?
		UInt consensus_sp = (UInt)ceil(consensus_sz/scan_resolution);

		StopWatch w;
		w.start();
		writeLog_(String("Determining potential network edges ..") );

		std::vector< DoubleReal > total_intensity; // total intensity of the selected spectra (for intensity-coverage ratio with matches)
		std::vector< std::pair<Size,Size> > pot_pairs; // name is program; will be reduced in run of the program
		total_intensity.reserve(experiment.size());
		pot_pairs.reserve((Size)(experiment.size()*(experiment.size()+1)/2));
		for(Size i = 0; i < experiment.size(); ++i) //can run up to size-1 because 2nd for-loop wont be exec so no harm done but total intensity filled for all
		{
			//fill pot_pairs
			for(Size j = i+1; j < experiment.size(); ++j)
			{
				if(fabs(experiment[j].getPrecursors().front().getMZ()-experiment[i].getPrecursors().front().getMZ()) < max_pm_diff+pm_tol+0.00001 )
				{
					if(fabs(experiment[i].getPrecursors().front().getMZ() > experiment[j].getPrecursors().front().getMZ()))
					{
						pot_pairs.push_back(std::pair<Size,Size>(j,i));
					}
					else
					{
						pot_pairs.push_back(std::pair<Size,Size>(i,j));
					}
				}
			}
			//fill total_intensity
			DoubleReal intens;
			for(Size j = 0; j < experiment[i].size(); ++j)
			{
				intens += experiment[i][j].getIntensity();
			}
			total_intensity.push_back(intens);
		}

		writeLog_(String(" .. preselected ") + pot_pairs.size() + String(" edges ..") );

		std::vector< std::pair<DoubleReal,DoubleReal> > pot_pairs_xcs; // the pairs xcorrelation scores
		std::vector< Size > stats_num(experiment.size()); // the number of scores for respective spec (for online calc of the var)
		std::vector< DoubleReal > means(experiment.size()), variances(experiment.size()), stddevs(experiment.size()); // the selected specs xcorrelation means, variation and std. deviation (unbiased gaussian)

		pot_pairs_xcs.reserve(pot_pairs.size());
		XCorrelation<Peak1D> x_corr;
		x_corr.setParameters(param);
		/// @todo set parameters for XCorr (and AntiSymetricAlignment and StarClusters?)
		for(Size i = 0; i < pot_pairs.size(); ++i)
		{
			DoubleReal best_score1, best_score2, best_shift;
			std::list<std::pair<Size,Size> > best_matches;
			x_corr.getXCorrelation(experiment[pot_pairs[i].first], experiment[pot_pairs[i].second], best_score1, best_score2, best_shift, best_matches);

			//calc mean and variation
			means[pot_pairs[i].first] = ( means[pot_pairs[i].first]*stats_num[pot_pairs[i].first] + best_score1 )/ (1+stats_num[pot_pairs[i].first]);
			means[pot_pairs[i].second] = ( means[pot_pairs[i].second]*stats_num[pot_pairs[i].second] + best_score2 )/ (1+stats_num[pot_pairs[i].second]);
			variances[pot_pairs[i].first] = ( variances[pot_pairs[i].first]*stats_num[pot_pairs[i].first] + (best_score1*best_score1) )/ (1+stats_num[pot_pairs[i].first]);
			variances[pot_pairs[i].second] = ( variances[pot_pairs[i].second]*stats_num[pot_pairs[i].second] + (best_score1*best_score2) )/ (1+stats_num[pot_pairs[i].second]);
			++stats_num[pot_pairs[i].first];
			++stats_num[pot_pairs[i].second];
		}

		//calc std.dev. sigma as sqrt of variance
		for(Size i = 0; i < experiment.size(); ++i)
		{
			stddevs[i]= sqrt( (means[i]/(means[i]-1)) * (variances[i]-means[i]*means[i]) );
		}

		//filter pairs not fitting ratio or gcdf
		std::vector< std::pair<Size,Size> >::iterator it_pairs = pot_pairs.begin();
		std::vector< std::pair<DoubleReal,DoubleReal> >::iterator it_pairs_xcs = pot_pairs_xcs.begin();
		while(it_pairs != pot_pairs.end())
		{
			if(  boost::math::cdf(boost::math::normal(means[it_pairs->first],stddevs[it_pairs->first]), it_pairs_xcs->first) >= 1-pairfilter_p_value and
					 boost::math::cdf(boost::math::normal(means[it_pairs->second],stddevs[it_pairs->second]), it_pairs_xcs->second) >= 1-pairfilter_p_value and
					 it_pairs_xcs->first/total_intensity[it_pairs->first] >= pairfilter_min_ratio and
					 it_pairs_xcs->second/total_intensity[it_pairs->second] >= pairfilter_min_ratio )
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

		writeLog_(String(" .. final selected ") + pot_pairs.size() + String(" edges ..") );


		w.stop();
		writeLog_(String(" done .. took") + String(w.getClockTime()) + String(" seconds"));
		w.reset();

		/*-------------------/
		// build alignments
		/-------------------*/
		w.start();
		writeLog_(String("Defining network edges ..") );

		AntisymetricAlignment<Peak1D> asa;
		asa.setParameters(param);

		PeakMap aligned_spectra; // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		std::vector<DoubleReal> mod_positions; // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		std::vector< std::pair<Size,Size> >& aligned_pairs = pot_pairs; // indices of aligned_spectra in original experiment
		for(Size i = 0; i < aligned_pairs.size(); ++i)
		{
			MSSpectrum<Peak1D> res_1, res_2;
			DoubleReal mod_pos, score;
			asa.getAntisymetricAlignment(res_1, res_2, score, mod_pos, experiment[aligned_pairs[i].first], experiment[aligned_pairs[i].second]);
			res_1.setMetaValue("original index", aligned_pairs[i].first);
			res_1.setMetaValue("paired with index", aligned_pairs[i].second);
			res_2.setMetaValue("original index", aligned_pairs[i].second);
			res_2.setMetaValue("paired with index", aligned_pairs[i].first);
			aligned_spectra.push_back(res_1);
			aligned_spectra.push_back(res_2);
			mod_positions.push_back(mod_pos);
			/// @improvement scoring? what to do with scoring? thin out network??
		}

		w.stop();
		writeLog_(String(" done .. took") + String(w.getClockTime()) + String(" seconds, defined ") + String(aligned_pairs.size()) + String(" edges"));
		w.reset();

		/*-------------------------------------------------------/
		// build consensuses i.e. network nodes i.e. star center
		/-------------------------------------------------------*/
		w.start();
		writeLog_(String("Defining star cluster and create consensuses in centers") );

		StarClusters<Peak1D> stars;
		stars.setParameters(param);

		std::map< Size, std::set<Size> > indices_of_i_in_aligned_spectra; // map key is index i to spec s in original experiment, map value is all indices to s' in aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		for(Size i = 0; i < aligned_pairs.size(); ++i)
		{
			indices_of_i_in_aligned_spectra[aligned_pairs[i].first].insert(i*2);
			indices_of_i_in_aligned_spectra[aligned_pairs[i].second].insert((i*2)+1);
		}

		//orientate the specs - only neccessary if the spectra in pairs came from antisymmetric path alignment
		for(std::map< Size, std::set<Size> >::const_iterator it = indices_of_i_in_aligned_spectra.begin(); it != indices_of_i_in_aligned_spectra.end(); ++it)
		{
			std::vector< MSSpectrum<Peak1D> > unmerged;
			unmerged.push_back(aligned_spectra[*(it->second.begin())]);

			MSSpectrum<Peak1D> base(aligned_spectra[*(it->second.begin())]), rev_base(aligned_spectra[*(it->second.begin())]);
			stars.reverseSpectrum(rev_base);

			std::set<Size>::const_iterator others = it->second.begin();
			++others;
			for(; others != it->second.end(); ++others)
			{
				MSSpectrum<Peak1D> current(aligned_spectra[*others]);
				//find matches to base in -parentmass_tolerance:peakmass_tolerance:+parentmass_tolerance shifts
				DoubleReal best_score_base = x_corr.bestMatchIntensity(base, current);

				//find matches to rev_base in -parentmass_tolerance:peakmass_tolerance:+parentmass_tolerance shifts
				DoubleReal best_score_rev_base = x_corr.bestMatchIntensity(rev_base, current);

				//take best, possibly reverse
				if(best_score_base < best_score_rev_base)
				{
					stars.reverseSpectrum(current);
				}
				unmerged.push_back(current);
			}

			//consensus building
			experiment[it->first] = BinnedSpectrum<Peak1D>(consensus_sz, consensus_sp,unmerged);
			experiment[it->first].setMetaValue("consensus", DataValue::EMPTY_VALUE);
		}

		w.stop();
		writeLog_(String(" done .. took") + String(w.getClockTime()) + String(" seconds, made ") + String(experiment.size()) + String(" consensuses"));
		writeLog_(String("all done"));
		w.reset();

		//~ save consensuses
		writeLog_(String("writing consensuses"));
		MzMLFile mzml;
		mzml.store(outputfile_name_consensus, experiment);

		/*-----------------------------/
		// save star edges and modpos
		/-----------------------------*/
		TextFile stars_txt;
		writeLog_(String("writing edges"));

		for(Size i = 0; i < aligned_pairs.size()-1; ++i)
		{
			stars_txt.push_back(String(aligned_pairs[i].first) +String('-') +String(aligned_pairs[i].second) +String('-') +String(mod_positions[i]) +String('\n'));
		}
		stars_txt.push_back(String(aligned_pairs[aligned_pairs.size()-1].first) +String('-') +String(aligned_pairs[aligned_pairs.size()-1].second) +String('-') +String(mod_positions[aligned_pairs.size()-1]) +String("\r\n"));

		stars_txt.store(outputfile_name_edges);
		//~ reading see TextFile_test
		return;
	}
	///

	void propagateNetwork_(MSExperiment<Peak1D>& experiment, Param& param, ConsensusMap& ids, std::vector< std::pair<Size,Size> >& aligned_pairs, std::vector< DoubleReal > mod_positions)
	{

		/*----------------------/
		// calc. max hop number
		/----------------------*/
		int max_hops = param.getValue("propagation_max_hops");

		/*--------------------/
		// build propagation
		/--------------------*/
		StarClusters<Peak1D> sc;
		sc.setParameters(param);
		sc.propagateNetwork(experiment, ids, aligned_pairs, mod_positions, max_hops);

		/*---------------------------/
		// score and export network
		/---------------------------*/
		/// @improvement dont punish to hard if spectrum does cover only a part of/or more than the sequence - e.g maybe mod that shortened the pep


		return;
	}
	///

	ExitCodes main_(int, const char**)
	{
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------

		String inputfile_name = getStringOption_("in");
		inputFileReadable_(inputfile_name);

		String ids_file = getStringOption_("ids");
		String network_file = getStringOption_("network");

		String dir = getStringOption_("directory");
		String prefix = getStringOption_("prefix");

		//~ parameter file for internal stuff read or default setting
		Param param= getParam_();

		if(ids_file.empty() and network_file.empty())
		{
			//~ make consensus case
			String outputfile_name_consensus = dir+prefix+String("-consensus.mzML");
			outputFileWritable_(outputfile_name_consensus);
			String outputfile_name_edges = dir+prefix+String("-edges.txt");
			outputFileWritable_(outputfile_name_edges);

			StopWatch w;
			w.start();

			writeLog_(String("Loading dataset from ") + inputfile_name +  String(" .. ") );
			MzMLFile mzml;
			MSExperiment<Peak1D> experiment; //MSExperiment<Peak1D>= PeakMap see StandardTypes
			mzml.load(inputfile_name,experiment);
			experiment.sortSpectra();
			w.stop();
			writeLog_(String("done  .. took") + String(w.getClockTime()) + String(" seconds, loaded ") + String(experiment.size()) + String(" spectra"));
			w.reset();

			//~ start consensus making
			spanNetwork_(experiment, param, outputfile_name_consensus, outputfile_name_edges);
		}
		else if(!ids_file.empty() and !network_file.empty())
		{
			//~ make propagation case
			inputFileReadable_(ids_file);
			inputFileReadable_(network_file);

			String outputfile_name_specnet = dir+prefix+String("-specnet.consensusXML");
			outputFileWritable_(outputfile_name_specnet);

			StopWatch w;
			w.start();

			writeLog_(String("Loading nodes from ") + inputfile_name +  String(" .. ") );
			MzMLFile mzml;
			MSExperiment<Peak1D> experiment; //MSExperiment<Peak1D>= PeakMap see StandardTypes
			mzml.load(inputfile_name,experiment);
			w.stop();
			writeLog_(String("done  .. took") + String(w.getClockTime()) + String(" seconds, loaded ") + String(experiment.size()) + String(" spectra"));
			w.reset();

			w.start();
			writeLog_(String("Loading network of stars from ") + network_file +  String(" ..") );

			std::vector< std::pair<Size,Size> > aligned_pairs;
			std::vector< DoubleReal > mod_positions;
			TextFile stars_txt;
			stars_txt.store(network_file);
			for(Size i = 0; i < stars_txt.size(); ++i)
			{
				std::vector< String > splits;
				stars_txt[i].split('-',splits);
				if(splits.size()!=3)
				{
					writeLog_(String("Invalid input, aborting!"));
					return PARSE_ERROR;
				}
				try
				{
					aligned_pairs.push_back(std::make_pair<Size,Size>(splits[0].toDouble(),splits[1].toDouble()));
					mod_positions.push_back(splits[2].toDouble());
				}
				catch(...)
				{
					writeLog_(String("Invalid input, aborting!"));
					return PARSE_ERROR;
				}
			}

			w.stop();
			writeLog_(String(" done  .. took") + String(w.getClockTime()) + String(" seconds, loaded ") + String(aligned_pairs.size()) + String(" edges"));
			w.reset();

			w.start();
			writeLog_(String("Loading identifications to nodes from ") + ids_file +  String(" ..") );

			//~ adopt identification
			ConsensusMap ids;
			try
			{
				IdXMLFile idxml;
				IDMapper idm;
				std::vector< ProteinIdentification > protein_ids;
				std::vector< PeptideIdentification > peptide_ids;
				idxml.load(ids_file, protein_ids, peptide_ids);
				idm.annotate(ids,peptide_ids,protein_ids);
			}
			catch(...)
			{
				writeLog_(String("Invalid input, aborting!"));
				return PARSE_ERROR;
			}

			w.stop();
			writeLog_(String(" done  .. took") + String(w.getClockTime()) + String(" seconds, loaded ") + String(ids.size()) + String(" ids"));
			w.reset();


			//~ start propagation
			propagateNetwork_(experiment, param, ids, aligned_pairs, mod_positions);
			ConsensusXMLFile cxml;
			cxml.store(outputfile_name_specnet, ids);
		}
		else
		{
			//~ invalid input case
			return PARSE_ERROR;
		}

		return EXECUTION_OK;
	}
	///
};

int main( int argc, const char** argv )
{
	SpectralNetworking tool;
	return tool.main(argc,argv);
}
/// @endcond
