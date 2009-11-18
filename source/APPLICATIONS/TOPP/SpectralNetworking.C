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
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <boost/math/distributions/normal.hpp> // for binomial distribution
#include <iostream>
#include <cmath>
#include <exception>

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
			: TOPPBase("SpectralNetworking","pipeline for spectral networking", false, "0.4beta")
		{

		}

 protected:

	Param getSubsectionDefaults_(const String& section) const
	{
		if (section == "general")
		{
			Param p(StarClusters<Peak1D>().getDefaults());
			//~ remove all starclusters specific i.e. all but parentmass and peak tolerance
			p.removeAll("m");
			return p;
		}
		if (section == "filter")
		{
			Param p;
			p.setValue("scan_resolution", 0.1, "Defines the resolution of the instrument used to create the spectra");
			p.setValue("max_pm_diff", 150.0, "Defines the maximal absolut (in Da) shift between the two spectra");
			p.setValue("consensus_bins", 0.3, "Defines the bin size (in Da) while makeing a consensus from two or more spectra");
			p.setValue("pairfilter_p_value", 0.05, "Defines the p_value while filtering optential spectral pairs");
			p.setValue("pairfilter_min_ratio", 0.4, "Defines the minimal matchratio while filtering optential spectral pairs");
			return p;
		}
		if (section == "spectra")
		{
			Param p(AntisymetricAlignment<Peak1D>().getDefaults());
			//~ remove all AntisymetricAlignment unspecific i.e. all but parentmass and peak tolerance
			p.removeAll("p");
			return p;
		}
		if (section == "network")
		{
			Param p(StarClusters<Peak1D>().getDefaults());
			//~ remove all starclusters unspecific i.e. parentmass and peak tolerance
			p.removeAll("p");
			return p;
		}

		//~ std::cout << section << std::endl;
		Param p;
		return p;
	}
	///

	void registerOptionsAndFlags_()
	{
		//~ file section
		registerInputFile_("in","<file>","","Input ms/ms data in mzML format",true);
		setValidFormats_("in",StringList::create("mzML"));
		registerStringOption_("prefix", "prefix for the outputfiles", "", "prefix for the outputfiles", true);
		registerStringOption_("directory", "output directory", "", "output directory", false);

		registerInputFile_("ids","<file>","","Input of ids to consensuses for propagation",false);
		setValidFormats_("ids",StringList::create("idXML"));
		registerInputFile_("network","<file>","","Input of network edges for propagation",false);

		registerSubsection_("general","General section");
		registerSubsection_("filter","Filter section");
		registerSubsection_("spectra","Spectra comparison section");
		registerSubsection_("network","Network section");

		//~ addEmptyLine_();
		//~ addText_("Parameters for the sections can only be given in the INI file.");
	}
	///

	//~ builds a spectral network ready for initial identification and subsequent propagation
	void spanNetwork_(MSExperiment<Peak1D>& experiment, String& outputfile_name_consensus, String& outputfile_name_edges)
	{
		//~ parameters
		DoubleReal pm_tol = getParam_().getValue("general:parentmass_tolerance");
		DoubleReal max_pm_diff = getParam_().getValue("filter:max_pm_diff");
		DoubleReal pairfilter_p_value = getParam_().getValue("filter:pairfilter_p_value");
		DoubleReal pairfilter_min_ratio = getParam_().getValue("filter:pairfilter_min_ratio");
		DoubleReal scan_resolution = getParam_().getValue("filter:scan_resolution");
		Real consensus_sz = getParam_().getValue("filter:consensus_bins"); /* e.g. peak_tolerance */;
		/// @improvement is floor or round more appropriate?
		UInt consensus_sp = (UInt)ceil(consensus_sz/scan_resolution);

		writeLog_(String("Determining network edges .. ") );

		StopWatch w;
		w.start();

		ProgressLogger logger;
		logger.setLogType(log_type_);

		std::vector< DoubleReal > total_intensity; // total intensity of the selected spectra (for intensity-coverage ratio with matches)
		std::vector< std::pair<Size,Size> > pot_pairs; // name is program; will be reduced in run of the program
		total_intensity.reserve(experiment.size());
		pot_pairs.reserve((Size)(experiment.size()*(experiment.size()+1)/2));
		logger.startProgress(0,experiment.size(),"preselecting edges");
		for(Size i = 0; i < experiment.size(); ++i) //can run up to size-1 because 2nd for-loop wont be exec so no harm done but total intensity filled for all
		{
			logger.setProgress(i);
			//fill pot_pairs
			for(Size j = i+1; j < experiment.size(); ++j)
			{
				DoubleReal pm_i = experiment[i].getPrecursors().front().getMZ();
				int c_i= experiment[i].getPrecursors().front().getCharge();
				DoubleReal pm_j = experiment[j].getPrecursors().front().getMZ();
				int c_j = experiment[j].getPrecursors().front().getCharge();
				/// @attention if the precursor charge is unknown, i.e. 0 best guess is its doubly charged
				(c_i==0)?c_i=2:c_i=c_i;
				(c_j==0)?c_j=2:c_j=c_j;
			/// @attention singly charged mass difference!
				Real pm_diff = (pm_j*c_j + (c_j-1)*Constants::PROTON_MASS_U)-(pm_i*c_i + (c_i-1)*Constants::PROTON_MASS_U);
				if(fabs(pm_diff) < max_pm_diff+pm_tol+0.00001 )
				{
					if(pm_diff < 0)
					{
						//~ pm_i > pm_j
						pot_pairs.push_back(std::pair<Size,Size>(j,i));
					}
					else
					{
						//~ pm_i <= pm_j
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
		logger.endProgress();
		writeLog_(String(" .. preselected ") + pot_pairs.size() + String(" edges ..") );

		std::vector< std::pair<DoubleReal,DoubleReal> > pot_pairs_xcs; // the pairs xcorrelation scores
		std::vector< Size > stats_num(experiment.size()); // the number of scores for respective spec (for online calc of the var)
		std::vector< DoubleReal > means(experiment.size()), squared_value_means(experiment.size()), std_devs(experiment.size()); // the selected specs xcorrelation means, variation and std. deviation (unbiased gaussian)

		pot_pairs_xcs.reserve(pot_pairs.size());

		std::vector<String> prefixes;
		prefixes.push_back("general:");
		prefixes.push_back("spectra:");
		Param xcorr_param;
		xcorr_param.insert("",getParam_().copy(prefixes[0],true));
		xcorr_param.insert("",getParam_().copy(prefixes[1],true));
		XCorrelation<Peak1D> x_corr;
		x_corr.setParameters(xcorr_param);
		logger.startProgress(0,pot_pairs.size()+experiment.size(),"evaluating selected edges");
		for(Size i = 0; i < pot_pairs.size(); ++i)
		{
			logger.setProgress(i);
			DoubleReal best_score1, best_score2, best_shift;
			std::list<std::pair<Size,Size> > best_matches;
			x_corr.getXCorrelation(experiment[pot_pairs[i].first], experiment[pot_pairs[i].second], best_score1, best_score2, best_shift, best_matches);

			//calc mean and variation
			means[pot_pairs[i].first] = ( means[pot_pairs[i].first]*stats_num[pot_pairs[i].first] + best_score1 )/ (1+stats_num[pot_pairs[i].first]);
			means[pot_pairs[i].second] = ( means[pot_pairs[i].second]*stats_num[pot_pairs[i].second] + best_score2 )/ (1+stats_num[pot_pairs[i].second]);
			squared_value_means[pot_pairs[i].first] = ( squared_value_means[pot_pairs[i].first]*stats_num[pot_pairs[i].first] + (best_score1*best_score1) )/ (1+stats_num[pot_pairs[i].first]);
			squared_value_means[pot_pairs[i].second] = ( squared_value_means[pot_pairs[i].second]*stats_num[pot_pairs[i].second] + (best_score1*best_score2) )/ (1+stats_num[pot_pairs[i].second]);
			++stats_num[pot_pairs[i].first];
			++stats_num[pot_pairs[i].second];

			pot_pairs_xcs.push_back(std::make_pair<DoubleReal,DoubleReal>(best_score1,best_score2));
		}

		//calc std.dev. sigma as sqrt of variance
		for(Size i = 0; i < experiment.size(); ++i)
		{
			logger.setProgress(pot_pairs.size()+i);
			/* debug std::cout << ((stats_num[i]/(stats_num[i]-1)) * (squared_value_means[i]-means[i]*means[i])) << std::endl; */
			std_devs[i]= sqrt( ((stats_num[i]/(stats_num[i]-1)) * (squared_value_means[i]-means[i]*means[i])) );
			/* debug std::cout << std_devs[i] << std::endl; */
		}
		logger.endProgress();

		//filter pairs not fitting ratio or gcdf
		std::vector< std::pair<Size,Size> >::iterator it_pairs = pot_pairs.begin();
		std::vector< std::pair<DoubleReal,DoubleReal> >::iterator it_pairs_xcs = pot_pairs_xcs.begin();
		logger.startProgress(0,pot_pairs.size(),"refining edge selection");
		while(it_pairs != pot_pairs.end())
		{
			logger.setProgress(it_pairs - pot_pairs.begin());
			/// @improvement : not erase but copy and reassign
			if(  boost::math::cdf(boost::math::normal(means[it_pairs->first],std_devs[it_pairs->first]), it_pairs_xcs->first) >= 1-pairfilter_p_value and
					 boost::math::cdf(boost::math::normal(means[it_pairs->second],std_devs[it_pairs->second]), it_pairs_xcs->second) >= 1-pairfilter_p_value and
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
		logger.endProgress();

		w.stop();
		writeLog_(String(" .. done, took ") + String(w.getClockTime()) + String(" seconds, selected ") + pot_pairs.size() + String(" edges.") );
		w.reset();

		/*-------------------/
		// build alignments
		/-------------------*/
		w.start();
		writeLog_(String("Defining network edges ..") );

		prefixes.clear();
		prefixes.push_back("general:");
		prefixes.push_back("spectra:");
		Param asa_param;
		asa_param.insert("",getParam_().copy(prefixes[0],true));
		asa_param.insert("",getParam_().copy(prefixes[1],true));
		AntisymetricAlignment<Peak1D> asa;
		asa.setParameters(asa_param);

		PeakMap aligned_spectra; // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		std::vector<DoubleReal> mod_positions; // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		std::vector< std::pair<Size,Size> >& aligned_pairs = pot_pairs; // indices of aligned_spectra in original experiment
		logger.startProgress(0,aligned_pairs.size(),"finding putative modification positions");
		for(Size i = 0; i < aligned_pairs.size(); ++i)
		{
			logger.setProgress(i);
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
		logger.endProgress();


		w.stop();
		writeLog_(String(" .. done, took") + String(w.getClockTime()) + String(" seconds, defined ") + String(aligned_pairs.size()) + String(" edges."));
		w.reset();

		/*-------------------------------------------------------/
		// build consensuses i.e. network nodes i.e. star center
		/-------------------------------------------------------*/
		w.start();
		writeLog_(String("Defining star cluster and create consensuses in centers ..") );


		prefixes.clear();
		prefixes.push_back("general:");
		prefixes.push_back("networks:");
		Param star_param;
		star_param.insert("",getParam_().copy(prefixes[0],true));
		star_param.insert("",getParam_().copy(prefixes[1],true));
		StarClusters<Peak1D> stars;
		stars.setParameters(star_param);

		std::map< Size, std::set<Size> > indices_of_i_in_aligned_spectra; // map key is index i to spec s in original experiment, map value is all indices to s' in aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		for(Size i = 0; i < aligned_pairs.size(); ++i)
		{
			indices_of_i_in_aligned_spectra[aligned_pairs[i].first].insert(i*2);
			indices_of_i_in_aligned_spectra[aligned_pairs[i].second].insert((i*2)+1);
		}

		//orientate the specs - only neccessary if the spectra in pairs came from antisymmetric path alignment
		logger.startProgress(0,indices_of_i_in_aligned_spectra.size(),"building consensuses");
		Size progr(0);
		for(std::map< Size, std::set<Size> >::const_iterator it = indices_of_i_in_aligned_spectra.begin(); it != indices_of_i_in_aligned_spectra.end(); ++it)
		{
			logger.setProgress(progr);
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
				DoubleReal best_score_base = x_corr.bestMatchIntensity(base, current); /// @improvement replace bestmatchintensity by something faster and cleverer

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
			++progr;
		}
		logger.endProgress();

		w.stop();
		writeLog_(String(" .. done, took ") + String(w.getClockTime()) + String(" seconds, made ") + String(experiment.size()) + String(" consensuses."));
		w.reset();

		//~ save consensuses
		writeLog_(String("Writing consensuses."));
		MzMLFile mzml;
		mzml.store(outputfile_name_consensus, experiment);

		/*-----------------------------/
		// save star edges and modpos
		/-----------------------------*/
		TextFile stars_txt;
		writeLog_(String("Writing edges."));

		for(Size i = 0; i < aligned_pairs.size()-1; ++i)
		{
			stars_txt.push_back(String(aligned_pairs[i].first) +String('-') +String(aligned_pairs[i].second) +String('-') +String(mod_positions[i]) +String('\n'));
		}
		stars_txt.push_back(String(aligned_pairs[aligned_pairs.size()-1].first) +String('-') +String(aligned_pairs[aligned_pairs.size()-1].second) +String('-') +String(mod_positions[aligned_pairs.size()-1]) +String("\r\n"));

		stars_txt.store(outputfile_name_edges);

		writeLog_(String("All done. Please proceed with SpectralNetworking after identification of the consensuses."));
		return;
	}
	///

	void propagateNetwork_(MSExperiment<Peak1D>& experiment, ConsensusMap& ids, std::vector< std::pair<Size,Size> >& aligned_pairs, std::vector< DoubleReal > mod_positions)
	{

		/*----------------------/
		// calc. max hop number
		/----------------------*/
		int max_hops = 3 /* param.getValue("propagation_max_hops") */;

		/*--------------------/
		// build propagation
		/--------------------*/
		std::vector<String> prefixes;
		prefixes.push_back("general:");
		prefixes.push_back("network:");
		Param star_param;
		star_param.insert("",getParam_().copy(prefixes[0],true));
		star_param.insert("",getParam_().copy(prefixes[1],true));

		StarClusters<Peak1D> sc;
		sc.setParameters(star_param);
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

		String prefix = getStringOption_("prefix");
		String dir = getStringOption_("directory");

		String ids_file = getStringOption_("ids");
		String network_file = getStringOption_("network");

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
			//~ MSExperiment<Peak1D>= PeakMap see StandardTypes
			MSExperiment<Peak1D> experiment;
			mzml.load(inputfile_name,experiment);
			//~ experiment.sortSpectra();

			//~ for scoring reasons this is mandatory
			Normalizer normalizer;
			Param p(normalizer.getParameters());
			p.setValue("method", "to_one");
			normalizer.setParameters(p);
			for(Size i = 0; i < experiment.size(); ++i)
			{
				normalizer.filterSpectrum(experiment[i]);
			}

			w.stop();
			writeLog_(String("done  .. took ") + String(w.getClockTime()) + String(" seconds, loaded ") + String(experiment.size()) + String(" spectra"));
			w.reset();

			//~ start consensus making
			spanNetwork_(experiment, outputfile_name_consensus, outputfile_name_edges);
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
			//~ MSExperiment<Peak1D>= PeakMap see StandardTypes
			MSExperiment<Peak1D> experiment;
			mzml.load(inputfile_name,experiment);

				//~ for scoring reasons this is mandatory
			Normalizer normalizer;
			Param p(normalizer.getParameters());
			p.setValue("method", "to_one");
			normalizer.setParameters(p);
			for(Size i = 0; i < experiment.size(); ++i)
			{
				normalizer.filterSpectrum(experiment[i]);
			}

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
			propagateNetwork_(experiment, ids, aligned_pairs, mod_positions);
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
