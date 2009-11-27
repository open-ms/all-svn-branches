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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <boost/math/distributions/normal.hpp> // for binomial distribution
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

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
			: TOPPBase("SpectralNetworking","pipeline for spectral networking", false, "1.0beta")
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
			p.setValue("pairs_p_value", 0.05, "Defines the p_value while filtering optential spectral pairs");
			p.setValue("pairs_min_ratio", 0.4, "Defines the minimal match ratio while filtering potential spectral pairs");
			p.setValue("pairs_min_matches", 4, "Defines the minimum number of matched peaks while filtering potential spectral pairs");
			return p;
		}
		if (section == "consensus")
		{
			Param p;
			p.setValue("bin_size", 0.3, "Defines the bin size (in Da) while makeing a consensus from two or more spectra");
			p.setValue("set_min_size", 3, "Defines the minimum size of a spectra set to make a consensus from");
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
		registerSubsection_("filter","Filtering spectra pairs section");
		registerSubsection_("spectra","Spectra comparison section");
		registerSubsection_("network","Network section");
		registerSubsection_("consensus","Consensus makeing section");

		//~ addEmptyLine_();
		//~ addText_("Parameters for the sections can only be given in the INI file.");
	}
	///

	//~ builds a spectral network ready for initial identification and subsequent propagation
	void spanNetwork_(MSExperiment<Peak1D>& experiment, String& outputfile_name_consensus, String& outputfile_name_edges)
	{
		//~ parameters
		DoubleReal pm_tol = getParam_().getValue("general:parentmass_tolerance");
		//~ DoubleReal peak_tol = getParam_().getValue("general:peak_tolerance");
		DoubleReal max_pm_diff = getParam_().getValue("filter:max_pm_diff");
		DoubleReal pairs_p_value = getParam_().getValue("filter:pairs_p_value");
		DoubleReal pairs_min_ratio = getParam_().getValue("filter:pairs_min_ratio");
		Size pairs_min_matches = getParam_().getValue("filter:pairs_min_matches");
		DoubleReal scan_resolution = getParam_().getValue("filter:scan_resolution");
		Real consensus_sz = getParam_().getValue("consensus:bin_size"); /* e.g. peak_tolerance */;
		/// @improvement is floor or round more appropriate?
		UInt consensus_sp = (UInt)ceil(consensus_sz/scan_resolution);
		Size consensus_min_size = getParam_().getValue("consensus:set_min_size");

		writeLog_(String("Determining network edges .. ") );

		MzMLFile mzml;
		StopWatch w;
		w.start();

		ProgressLogger logger;
		logger.setLogType(log_type_);

		std::vector< DoubleReal > total_intensity; // total intensity of the selected spectra (for intensity-coverage ratio with matches)
		std::vector< std::pair<Size,Size> > pot_pairs; // name is program; will be reduced in run of the program
		total_intensity.reserve(experiment.size());
		pot_pairs.reserve((Size)(experiment.size()*(experiment.size()+1)/2));
		DoubleReal average_peak_num(0);
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
					if(fabs(pm_diff)<pm_tol)
					{
						//~ pm_i = pm_j
						pot_pairs.push_back(std::pair<Size,Size>(i,j));
					}
					else if(pm_diff < 0)
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
			average_peak_num += experiment[i].size();
		}
		average_peak_num /= (DoubleReal)experiment.size();
		logger.endProgress();
		writeLog_(String(".. dataset with average peak number of ") + String(average_peak_num) + String(" per spectrum ..") );
		writeLog_(String(".. preselected ") + pot_pairs.size() + String(" edges ..") );

		std::vector< std::pair<DoubleReal,DoubleReal> > pot_pairs_xcs; // the pairs xcorrelation scores
		pot_pairs_xcs.reserve(pot_pairs.size());

		std::vector<String> prefixes;
		prefixes.push_back("general:");
		prefixes.push_back("spectra:");
		Param xcorr_param;
		xcorr_param.insert("",getParam_().copy(prefixes[0],true));
		xcorr_param.insert("",getParam_().copy(prefixes[1],true));
		XCorrelation<Peak1D> x_corr;
		x_corr.setParameters(xcorr_param);

		std::vector< boost::accumulators::accumulator_set<DoubleReal, boost::accumulators::stats<boost::accumulators::tag::variance> > > xcorr_accumulators(experiment.size());
		std::vector<Size> edge_selection;
		logger.startProgress(0,pot_pairs.size(),"evaluating selected edges");
		DoubleReal average_best_matches(0);
		for(Size i = 0; i < pot_pairs.size(); ++i)
		{
			logger.setProgress(i);
			DoubleReal best_score1, best_score2, best_shift;
			std::list<std::pair<Size,Size> > best_matches;
			/// @improvement use xcorr with AND without shift, accept only those with sufficient combined match numbers
			x_corr.getXCorrelation(experiment[pot_pairs[i].first], experiment[pot_pairs[i].second], best_score1, best_score2, best_shift, best_matches);
			average_best_matches += (DoubleReal)best_matches.size();
			//~ half the peaks number of average size of the two spectra times the given ratio is the minimum match number
			/// @improvement make pairs_min_ratio an advanced parameter
			if((DoubleReal)best_matches.size()/((DoubleReal)(experiment[pot_pairs[i].first].size()+experiment[pot_pairs[i].second].size())/DoubleReal(2)) >=  pairs_min_ratio)
			{
				xcorr_accumulators[pot_pairs[i].first](best_score1);
				xcorr_accumulators[pot_pairs[i].second](best_score2);
				pot_pairs_xcs.push_back(std::make_pair<DoubleReal,DoubleReal>(best_score1,best_score2));
				edge_selection.push_back(i);
			}
		}
		for(Size i = 0; i < edge_selection.size(); ++i)
		{
			//~ edge_selection indices are alwas >= the ones to pot_pairs so no collision expected
			pot_pairs[i] = pot_pairs[edge_selection[i]];
		}
		average_best_matches /= (DoubleReal)pot_pairs.size();
		pot_pairs.resize(edge_selection.size());
		logger.endProgress();
		writeLog_(String(".. preselected edges having a average match size of ") + String(average_best_matches) + String(" ..") );
		writeLog_(String(".. ") + pot_pairs.size() + String(" edges passed minimum ratio evaluation ..") );

		//filter pairs not gcdf
		logger.startProgress(0,pot_pairs.size() + experiment.size(),"refining edge selection");
		std::vector<DoubleReal> means(experiment.size()), stddevs(experiment.size());
		for(Size i = 0; i < experiment.size(); ++i)
		{
			logger.setProgress(i);
			/// @improvement make nicer than try/catch!
			means[i] = boost::accumulators::mean(xcorr_accumulators[i]);
			stddevs[i] = sqrt(boost::accumulators::variance(xcorr_accumulators[i]));
		}
		//~ std::vector< std::pair<Size,Size> >::iterator it_pairs = pot_pairs.begin();
		//~ std::vector< std::pair<DoubleReal,DoubleReal> >::iterator it_pairs_xcs = pot_pairs_xcs.begin();
		edge_selection.clear();
		Size edge_without_stddev(0);
		for(Size i = 0; i < pot_pairs.size(); ++i)
		{
			logger.setProgress(i+experiment.size());
			try{
				if(  boost::math::cdf(boost::math::normal(means[pot_pairs[i].first] ,stddevs[pot_pairs[i].first]), pot_pairs_xcs[i].first) >= 1-pairs_p_value and
						 boost::math::cdf(boost::math::normal(means[pot_pairs[i].second] ,stddevs[pot_pairs[i].second]), pot_pairs_xcs[i].second) >= 1-pairs_p_value )
				{
					edge_selection.push_back(i);
				}
			}catch(...){ ++edge_without_stddev; }
		}
		writeLog_(String(edge_without_stddev) + String(" edges had a spec without stddev anyway ..") );
		for(Size i = 0; i < edge_selection.size(); ++i)
		{
			//~ edge_selection indices are alwas >= the ones to pot_pairs so no collision expected
			pot_pairs[i] = pot_pairs[edge_selection[i]];
		}
		pot_pairs.resize(edge_selection.size());
		logger.endProgress();
		w.stop();
		writeLog_(String(".. done, took ") + String(w.getClockTime()) + String(" seconds, selected ") + pot_pairs.size() + String(" edges.") );
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

		MSExperiment<Peak1D> aligned_spectra; // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		std::vector<DoubleReal> mod_positions; // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		std::vector< std::pair<Size,Size> >& aligned_pairs = pot_pairs; // indices of aligned_spectra in original experiment
		logger.startProgress(0,aligned_pairs.size(),"finding putative modification positions");
		edge_selection.clear();
		for(Size i = 0; i < aligned_pairs.size(); ++i)
		{
			logger.setProgress(i);
			MSSpectrum<Peak1D> res_1, res_2;
			DoubleReal mod_pos, score;
			asa.getAntisymetricAlignment(res_1, res_2, score, mod_pos, experiment[aligned_pairs[i].first], experiment[aligned_pairs[i].second]);
			/// @improvement make pairs_min_matches an advanced parameter
			if((DoubleReal)res_1.size()/((DoubleReal)(experiment[aligned_pairs[i].first].size()+experiment[aligned_pairs[i].second].size())/DoubleReal(2)) >=  pairs_min_ratio and res_1.size() >= pairs_min_matches)
			{
				edge_selection.push_back(i);
				/// @improvement pot. speedup by removing MetaValues ?!
				res_1.setMetaValue("original index", aligned_pairs[i].first);
				res_1.setMetaValue("paired with index", aligned_pairs[i].second);
				res_2.setMetaValue("original index", aligned_pairs[i].second);
				res_2.setMetaValue("paired with index", aligned_pairs[i].first);
				aligned_spectra.push_back(res_1);
				aligned_spectra.push_back(res_2);
				mod_positions.push_back(mod_pos);
			}
		}
		for(Size i = 0; i < edge_selection.size(); ++i)
		{
			//~ edge_selection indices are alwas >= the ones to aligned_pairs so no collision expected
			aligned_pairs[i] = aligned_pairs[edge_selection[i]];
		}
		aligned_pairs.resize(edge_selection.size());
		logger.endProgress();

		try{
		String alignments(outputfile_name_edges + String(".alignmentspecs.mzML"));
		mzml.store(alignments , aligned_spectra);
		}catch(...){writeLog_(String(".. aligned spectra not separately saved .."));}

		w.stop();
		writeLog_(String(".. done, took ") + String(w.getClockTime()) + String(" seconds, defined ") + String(aligned_pairs.size()) + String(" edges."));
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
		Size progr(0), made_consensuses(0);
		for(std::map< Size, std::set<Size> >::const_iterator it = indices_of_i_in_aligned_spectra.begin(); it != indices_of_i_in_aligned_spectra.end(); ++it)
		{
			logger.setProgress(progr);
			/// @improvement make consensus_min_size an advanced parameter
			if(it->second.size()>=consensus_min_size and it->second.size()>1)
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
					/// @improvement replace bestmatchintensity by something cleverer (and maybe faster? but this is hard!)
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
				String nat_id(experiment[it->first].getNativeID());
				experiment[it->first] = BinnedSpectrum<Peak1D>(consensus_sz, consensus_sp,unmerged);
				experiment[it->first].setMetaValue("consensus", DataValue::EMPTY_VALUE);
				experiment[it->first].setNativeID(nat_id);
				++made_consensuses;
			}
			++progr;
		}
		logger.endProgress();

		w.stop();
		writeLog_(String(".. done, took ") + String(w.getClockTime()) + String(" seconds, made ") + String(made_consensuses) + String(" consensuses."));
		w.reset();

		//~ save consensuses
		writeLog_(String("Writing consensuses."));
		mzml.store(outputfile_name_consensus, experiment);

		/*-----------------------------/
		// save star edges and modpos
		/-----------------------------*/
		TextFile stars_txt;

		if(aligned_pairs.size()>0)
		{
			writeLog_(String("Writing edges."));
			for(Size i = 0; i < aligned_pairs.size()-1; ++i)
			{
				/// @improvement save also the edges relation (i.e. equal/modified/ladder)
				stars_txt.push_back(String(aligned_pairs[i].first) +String('#') +String(aligned_pairs[i].second) +String('#') +String(mod_positions[i]) +String('\n'));
			}
			stars_txt.push_back(String(aligned_pairs[aligned_pairs.size()-1].first) +String('#') +String(aligned_pairs[aligned_pairs.size()-1].second) +String('#') +String(mod_positions[aligned_pairs.size()-1]) +String("\r\n"));

			stars_txt.store(outputfile_name_edges);
		}
		else
		{
			writeLog_(String("No edges."));
		}
		writeLog_(String("All done. Please proceed with SpectralNetworking after identification of the consensuses."));
		return;
	}
	///

	void propagateNetwork_(MSExperiment<Peak1D>& experiment, std::vector< PeptideIdentification > peptide_ids, std::vector< ProteinIdentification > protein_ids, std::vector< std::pair<Size,Size> >& aligned_pairs, std::vector< DoubleReal > mod_positions, ConsensusMap& ids )
	{
		std::vector<String> prefixes;
		prefixes.push_back("general:");
		prefixes.push_back("network:");
		Param star_param;
		star_param.insert("",getParam_().copy(prefixes[0],true));
		star_param.insert("",getParam_().copy(prefixes[1],true));

		/*----------------------/
		// calc. max hop number
		/----------------------*/
		/// @improvement calc max hop(from average connectivity) and set in star_param and make boolean use intelligent hop_calc an advanced parameter
		//~ int max_hops = param.getValue("network:max_hops");

		/*--------------------/
		// build propagation
		/--------------------*/
		StarClusters<Peak1D> sc;
		sc.setParameters(star_param);
		sc.propagateNetwork(experiment, peptide_ids, protein_ids, aligned_pairs, mod_positions, ids);

		/*------------/
		// score etc.
		/------------*/
		/// @improvement add featurehandles for those edges with fabs(pm_diff) < pm_tol and sufficient match and copy peptidehits if no id
		/// @improvement dont punish to hard if spectrum does cover only a part of/or more than the sequence - e.g maybe mod that shortened the pep
		/// @improvement txt output of peptide laddering sequences
		/// @improvement prepare featurehandles(metainfo?!) for colorization/different linestyle
		Size c(0);
		for(ConsensusMap::Iterator ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
		{
			if(ids_it->getPeptideIdentifications().front().empty())
			{
				ids_it->getPeptideIdentifications().clear();
				++c;
			}
		}
		writeLog_(String("unknowns ") + String(c));

		std::vector<String> colors;
		colors.push_back("#00FFFF");
		colors.push_back("#000000");
		colors.push_back("#0000FF");
		colors.push_back("#FF00FF");
		colors.push_back("#008000");
		colors.push_back("#808080");
		colors.push_back("#00FF00");
		colors.push_back("#800000");
		colors.push_back("#000080");
		colors.push_back("#808000");
		colors.push_back("#800080");
		colors.push_back("#FF0000");
		colors.push_back("#C0C0C0");
		colors.push_back("#008080");
		colors.push_back("#FFFF00");
		for(ConsensusMap::Iterator ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
		{
			if(ids_it->metaValueExists("network id"))
			{
				ids_it->setMetaValue("color", colors[(Size)ids_it->getMetaValue("network id")%colors.size()]);
				//~ set intensity for coloration
				ids_it->setIntensity(1000*((Size)ids_it->getMetaValue("network id"))%colors.size());
			}
		}

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
			writeLog_(String("done  .. took ") + String(w.getClockTime()) + String(" seconds, loaded ") + String(experiment.size()) + String(" spectra"));
			w.reset();

			w.start();
			writeLog_(String("Loading network of stars from ") + network_file +  String(" ..") );

			std::vector< std::pair<Size,Size> > aligned_pairs;
			std::vector< DoubleReal > mod_positions;
			TextFile stars_txt;
			stars_txt.load(network_file);
			for(Size i = 0; i < stars_txt.size(); ++i)
			{
				std::vector< String > splits;
				stars_txt[i].split('#',splits);
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
			writeLog_(String("done  .. took ") + String(w.getClockTime()) + String(" seconds, loaded ") + String(aligned_pairs.size()) + String(" edges"));
			w.reset();

			w.start();
			writeLog_(String("Loading identifications to nodes from ") + ids_file +  String(" ..") );

			IdXMLFile idxml;
			std::vector< ProteinIdentification > protein_ids;
			std::vector< PeptideIdentification > peptide_ids;
			try{
				idxml.load(ids_file, protein_ids, peptide_ids);
			}catch(...){
				writeLog_(String("Invalid input, aborting!"));
				return PARSE_ERROR;
			}

			w.stop();
			writeLog_(String(" done  .. took ") + String(w.getClockTime()) + String(" seconds, loaded ") + String(peptide_ids.size()) + String(" ids"));
			w.reset();

			w.start();

			//~ start propagation
			ConsensusMap ids;
			propagateNetwork_(experiment, peptide_ids, protein_ids, aligned_pairs, mod_positions, ids);
			ConsensusXMLFile cxml;
			writeLog_(String("Saving propagated network ..") );
			cxml.store(outputfile_name_specnet, ids);

			w.stop();
			writeLog_(String("All done. Saving took ") + String(w.getClockTime()) + String(" seconds. Ready to investigate in TOPPView."));
			w.reset();
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
