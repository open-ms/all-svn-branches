#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/XCorrelation.h>
#include <OpenMS/COMPARISON/SPECTRA/AntisymetricAlignment.h>
#include <OpenMS/COMPARISON/CLUSTERING/StarClusters.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

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
			: TOPPBase("SpectralNetworking","pipeline for spectral networking", false, "1.3beta")
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
			p.removeAll("f");
			return p;
		}
		if (section == "filter")
		{
			Param p;
			p.setValue("max_pm_diff", 150.0, "Defines the maximal absolut (in Da) shift between the two spectra");
			p.setValue("pairs_p_value", 0.05, "Defines the p_value while filtering potential spectral pairs");
			p.setValue("pairs_min_ratio", 0.4, "Defines the minimal match ratio while filtering potential spectral pairs");
			p.setValue("pairs_min_matches", 4, "Defines the minimum number of matched peaks while filtering potential spectral pairs");
			p.setValue("near_equal_min_match ratio", 1.1 , "Defines the minimum ratio of the number of matched peaks and intensities while clustering near equal spectral pairs");
			return p;
		}
		if (section == "consensus")
		{
			Param p;
			p.setValue("bin_size", 0.3, "Defines the bin size (in Da) while binning in makeing a consensus from two or more spectra");
			p.setValue("bin_spread", 0, "Defines the bin spread to both sides of a peak while binning in makeing a consensus from two or more spectra");
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

		registerInputFile_("consensus","<file>","","Processed input ms/ms data in mzML format",false);
		setValidFormats_("consensus",StringList::create("mzML"));
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

	void addZeroIonPeaks_(MSExperiment<Peak1D>& experiment)
	{
		DoubleReal peak_tol = getParam_().getValue("general:peak_tolerance");
		for(Size i = 0; i < experiment.size(); ++i)
		{
			/// @attention uncharged mass
			DoubleReal pm(experiment[i].getPrecursors().front().getUnchargedMass());

			//~ b0 = |H|; bk = pm - |OH|; y0 = |H3O|; yk = pm - |H|
			EmpiricalFormula f_H("H"), f_O("O");
			DoubleReal b0(f_H.getMonoWeight()+Constants::PROTON_MASS_U), bk(pm - f_H.getMonoWeight() - f_O.getMonoWeight()+Constants::PROTON_MASS_U), y0(3*f_H.getMonoWeight() + f_O.getMonoWeight()+Constants::PROTON_MASS_U), yk(pm - f_H.getMonoWeight()+Constants::PROTON_MASS_U);

			//~ b0
			MSSpectrum<Peak1D>::Iterator lo = experiment[i].MZBegin(b0-peak_tol);
			MSSpectrum<Peak1D>::Iterator hi = experiment[i].MZEnd(lo,b0+peak_tol,experiment[i].end());
			if(lo==hi)
			{
				Peak1D tmp;
				tmp.setPosition(b0);
				tmp.setIntensity(0.5);
				experiment[i].push_back(tmp);
			}
			//~ bk
			lo = experiment[i].MZBegin(bk-peak_tol);
			hi = experiment[i].MZEnd(lo,bk+peak_tol,experiment[i].end());
			if(lo==hi)
			{
				Peak1D tmp;
				tmp.setPosition(bk);
				tmp.setIntensity(0.5);
				experiment[i].push_back(tmp);
			}

			//~ y0
			lo = experiment[i].MZBegin(y0-peak_tol);
			hi = experiment[i].MZEnd(lo,y0+peak_tol,experiment[i].end());
			if(lo==hi)
			{
				Peak1D tmp;
				tmp.setPosition(y0);
				tmp.setIntensity(0.5);
				experiment[i].push_back(tmp);
			}
			//~ yk
			lo = experiment[i].MZBegin(yk-peak_tol);
			hi = experiment[i].MZEnd(lo,yk+peak_tol,experiment[i].end());
			if(lo==hi)
			{
				Peak1D tmp;
				tmp.setPosition(yk);
				tmp.setIntensity(0.5);
				experiment[i].push_back(tmp);
			}
		}
	}
	///

	void clusterInTolerance_(MSExperiment<Peak1D>& experiment, std::vector< std::set<Size> >& cluster)
	{
		Real consensus_sz = getParam_().getValue("consensus:bin_size");
		UInt consensus_sp = getParam_().getValue("consensus:bin_spread");
		DoubleReal pm_tol = getParam_().getValue("general:parentmass_tolerance");
		DoubleReal min_ratio = getParam_().getValue("filter:near_equal_min_match ratio");

		cluster.clear();
		if(min_ratio>1)
		{
			return;
		}

		std::vector<Size> cluster_mirror(experiment.size());
		for(Size i = 0; i < experiment.size(); ++i)
		{
			cluster_mirror[i] = i;
		}

		DoubleReal auto_si(0);
		//~ SpectrumAlignmentScore sas;
		SpectrumAlignment sa;
		for(Size i = 0; i < experiment.size(); ++i)
		{
			//~ auto_si == sas(experiment[i],experiment[i]);
			for(Size j = i+1; j < experiment.size(); ++j)
			{
				/// @attention uncharged masses
				DoubleReal pm_i(experiment[i].getPrecursors().front().getUnchargedMass());
				DoubleReal pm_j(experiment[j].getPrecursors().front().getUnchargedMass());
				Real pm_diff = (pm_j-pm_i);

				if(fabs(pm_diff) < pm_tol)
				{
					/// @todo cut off at right score
					DoubleReal auto_sj(0), sim(0), auto_mean(0), ratio(0), peak_mean(0);
					//~ auto_sj == sas(experiment[j],experiment[j]);
					//~ sim = sas(experiment[i],experiment[j]);
					//~ auto_mean = (auto_si+auto_sj)/2;
					peak_mean = ((DoubleReal)experiment[i].size()+(DoubleReal)experiment[j].size())/2;
					std::vector< std::pair< Size, Size > > alignment;
					sa.getSpectrumAlignment(alignment, experiment[i],experiment[j]);
					//~ ratio = sim/auto_mean;
					ratio = (DoubleReal)alignment.size()/peak_mean;
					std::cout << "ratio" << ratio << std::endl;
					if(ratio>min_ratio)
					{
						if(cluster_mirror[j]<cluster_mirror[i])
						{
							for(Size k = 0; k < cluster_mirror.size(); ++k)
							{
								if(cluster_mirror[k] == cluster_mirror[i])
								{
									cluster_mirror[k] = cluster_mirror[j];
								}
							}
						}
						else
						{
							cluster_mirror[j]=cluster_mirror[i];
						}
					}
				}
			}
		}

		/// @improvement remove clusters < 3 'cs of sufficient specs to consens on?
		cluster.clear();
		cluster.resize(experiment.size());
		for(Size k = 0; k < cluster_mirror.size(); ++k)
		{
			cluster[cluster_mirror[k]].insert(k);
		}
		std::vector< std::set<Size> >::iterator it = cluster.begin();
		while(it != cluster.end())
		{
			if(it->size()<2)
			{
				it = cluster.erase(it);
			}
			else
			{
				++it;
			}
		}
		std::cout << ".. " << cluster.size() << " consensuses"<< std::endl;

		//consensus building
		for(Size k = 0; k < cluster.size(); ++k)
		{
			std::vector< MSSpectrum<Peak1D> > unmerged;
			std::set<Size>::iterator it = cluster[k].begin();
			for(it; it != cluster[k].end(); ++it)
			{
				unmerged.push_back(experiment[*it]);
			}
			experiment[*(cluster[k].begin())] = BinnedSpectrum<Peak1D>(consensus_sz, consensus_sp, 0.1, 0.01, unmerged);
			experiment[*(cluster[k].begin())].setMetaValue("consensus", DataValue::EMPTY_VALUE);
		}
	}
	///

	//~ builds a spectral network ready for initial identification and subsequent propagation
	void spanNetwork_(MSExperiment<Peak1D>& experiment, String& outputfile_name_consensus, String& outputfile_name_edges, std::set<Size>& no_edges)
	{
		//~ parameters
		DoubleReal pm_tol = getParam_().getValue("general:parentmass_tolerance");
		//~ DoubleReal peak_tol = getParam_().getValue("general:peak_tolerance");
		DoubleReal max_pm_diff = getParam_().getValue("filter:max_pm_diff");
		DoubleReal pairs_p_value = getParam_().getValue("filter:pairs_p_value");
		DoubleReal pairs_min_ratio = getParam_().getValue("filter:pairs_min_ratio");
		Size pairs_min_matches = getParam_().getValue("filter:pairs_min_matches");
		Real consensus_sz = getParam_().getValue("consensus:bin_size"); /* e.g. peak_tolerance */;
		UInt consensus_sp = getParam_().getValue("consensus:bin_spread");
		Size consensus_min_size = getParam_().getValue("consensus:set_min_size");

		writeLog_(String("Determining network edges .. ") );

		MzMLFile mzml;
		StopWatch w;
		w.start();

		ProgressLogger logger;
		logger.setLogType(log_type_);

		std::vector< DoubleReal > total_intensity; // total intensity of the selected spectra (for intensity-coverage ratio with matches)
		Size clustered_size(experiment.size()-no_edges.size()); // if clustering was applied
		total_intensity.reserve(clustered_size);
		//~ pot_pairs.reserve((Size)((clustered_size)*(clustered_size+1)/2));
		std::vector< std::vector< std::pair<Size,Size> > > pot_pairs(experiment.size()); // name is program; will be reduced in run of the program

		DoubleReal average_peak_num(0);
		logger.startProgress(0,experiment.size(),"preselecting edges");
		for(Size i = 0; i < experiment.size(); ++i) //can run up to size-1 because 2nd for-loop wont be exec so no harm done but total intensity filled for all
		{
			logger.setProgress(i);

			std::set<Size>::iterator it_i = no_edges.find(i);
			if(it_i != no_edges.end())
			{
				continue;
			}

			//fill pot_pairs
			for(Size j = i+1; j < experiment.size(); ++j)
			{
				std::set<Size>::iterator it_j = no_edges.find(j);
				if(it_j != no_edges.end())
				{
					continue;
				}

				/// @attention uncharged masses
				DoubleReal pm_i(experiment[i].getPrecursors().front().getUnchargedMass());
				DoubleReal pm_j(experiment[j].getPrecursors().front().getUnchargedMass());
				Real pm_diff = (pm_j-pm_i);

				if(fabs(pm_diff) < max_pm_diff+pm_tol+0.00001 )
				{
					if(pm_diff < 0)
					{
						//~ pm_i > pm_j
						pot_pairs[i].push_back(std::pair<Size,Size>(j,i));
					}
					else
					{
						//~ pm_i <= pm_j
							pot_pairs[i].push_back(std::pair<Size,Size>(i,j));
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

 		average_peak_num /= (DoubleReal)clustered_size;
		logger.endProgress();
		Size edge_num(0);
		for(Size oi = 0; oi < pot_pairs.size(); ++oi)
		{
			edge_num+=pot_pairs[oi].size();
		}
		writeLog_(String(".. dataset with average peak number of ") + String(average_peak_num) + String(" per spectrum ..") );
		writeLog_(String(".. preselected ") + edge_num + String(" edges ..") );

		std::vector< std::vector< std::pair<DoubleReal,DoubleReal> > > pot_pairs_xcs(pot_pairs.size()); // the pairs xcorrelation scores

		std::vector<String> prefixes;
	prefixes.push_back("general:");
		prefixes.push_back("spectra:");
		Param xcorr_param;
		xcorr_param.insert("",getParam_().copy(prefixes[0],true));
		xcorr_param.insert("",getParam_().copy(prefixes[1],true));
		xcorr_param.remove("sv_penalty");  // these are antisymetric alignment only
		xcorr_param.remove("dif_penalty"); // these are antisymetric alignment only
		XCorrelation<Peak1D> x_corr;
		x_corr.setParameters(xcorr_param);
		std::vector<Size> edge_selection;

		std::vector< boost::accumulators::accumulator_set<DoubleReal, boost::accumulators::stats<boost::accumulators::tag::variance> > > xcorr_accumulators(experiment.size());
		logger.startProgress(0,pot_pairs.size(),"evaluating selected edges");
		DoubleReal average_best_matches(0);
		for(Size oi = 0; oi < pot_pairs.size(); ++oi)
		{
			logger.setProgress(oi);
			for(Size i = 0; i < pot_pairs[oi].size(); ++i)
			{
				DoubleReal best_score1_sf, best_score2_sf, best_score1_st, best_score2_st, best_shift;
				std::list<std::pair<Size,Size> > best_matches_sf, best_matches_st ;
				x_corr.getXCorrelation(experiment[pot_pairs[oi][i].first], experiment[pot_pairs[oi][i].second], best_score1_sf, best_score2_sf, best_shift, best_matches_sf,false);
				x_corr.getXCorrelation(experiment[pot_pairs[oi][i].first], experiment[pot_pairs[oi][i].second], best_score1_st, best_score2_st, best_shift, best_matches_st,true);
				// half the peaks number of average size of the two spectra times the given ratio is the minimum match number
				/// @improvement make pairs_min_ratio an advanced parameter
				if(best_matches_sf > best_matches_st and (DoubleReal)best_matches_sf.size()/((DoubleReal)(experiment[pot_pairs[oi][i].first].size()+experiment[pot_pairs[oi][i].second].size())/DoubleReal(2)) >=  pairs_min_ratio)
				{
					xcorr_accumulators[pot_pairs[oi][i].first](best_score1_sf);
					xcorr_accumulators[pot_pairs[oi][i].second](best_score2_sf);
					pot_pairs_xcs[oi].push_back(std::make_pair<DoubleReal,DoubleReal>(best_score1_sf*100,best_score2_sf*100));
					//~ pot_pairs_xcs[oi].push_back(std::make_pair<DoubleReal,DoubleReal>(best_matches_sf.size(),best_matches_sf.size()));
					edge_selection.push_back(i);
					average_best_matches += best_matches_sf.size();
				}
				else if((DoubleReal)best_matches_st.size()/((DoubleReal)(experiment[pot_pairs[oi][i].first].size()+experiment[pot_pairs[oi][i].second].size())/DoubleReal(2)) >=  pairs_min_ratio)
				{
					xcorr_accumulators[pot_pairs[oi][i].first](best_score1_st);
					xcorr_accumulators[pot_pairs[oi][i].second](best_score2_st);
					pot_pairs_xcs[oi].push_back(std::make_pair<DoubleReal,DoubleReal>(best_score1_st*100,best_score2_st*100));
					//~ pot_pairs_xcs[oi].push_back(std::make_pair<DoubleReal,DoubleReal>(best_matches_st.size(),best_matches_st.size()));
					edge_selection.push_back(i);
					average_best_matches += best_matches_st.size();
				}
			}
			for(Size i = 0; i < edge_selection.size(); ++i)
			{
				//~ edge_selection indices are alwas >= the ones to pot_pairs so no collision expected
				pot_pairs[oi][i] = pot_pairs[oi][edge_selection[i]];
			}
			pot_pairs[oi].resize(edge_selection.size());
			edge_selection.clear();
		}

		logger.endProgress();
		edge_num = 0;
		for(Size oi = 0; oi < pot_pairs.size(); ++oi)
		{
			edge_num+=pot_pairs[oi].size();
		}
		average_best_matches /= (DoubleReal)edge_num;
		writeLog_(String(".. preselected edges having a average match size of ") + String(average_best_matches) + String(" ..") );
		writeLog_(String(".. ") + edge_num + String(" edges passed minimum ratio evaluation ..") );

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
		for(Size oi = 0; oi < pot_pairs.size(); ++oi)
		{
			logger.setProgress(oi+experiment.size());
			for(Size i = 0; i < pot_pairs[oi].size(); ++i)
			{
				try{
					if(  boost::math::cdf(boost::math::normal(means[pot_pairs[oi][i].first] ,stddevs[pot_pairs[oi][i].first]), pot_pairs_xcs[oi][i].first) >= 1-pairs_p_value and
							 boost::math::cdf(boost::math::normal(means[pot_pairs[oi][i].second] ,stddevs[pot_pairs[oi][i].second]), pot_pairs_xcs[oi][i].second) >= 1-pairs_p_value )
					{
						edge_selection.push_back(i);
					}
				}catch(...){ ++edge_without_stddev; }
			}
			for(Size i = 0; i < edge_selection.size(); ++i)
			{
				//~ edge_selection indices are alwas >= the ones to pot_pairs so no collision expected
				pot_pairs[oi][i] = pot_pairs[oi][edge_selection[i]];
			}
			pot_pairs[oi].resize(edge_selection.size());
			edge_selection.clear();
		}
		writeLog_(String(edge_without_stddev) + String(" edges had a spec without stddev anyway ..") );
		logger.endProgress();
		w.stop();
		edge_num = 0;
		for(Size oi = 0; oi < pot_pairs.size(); ++oi)
		{
			edge_num+=pot_pairs[oi].size();
		}
		writeLog_(String(".. done, took ") + String(w.getClockTime()) + String(" seconds, selected ") + edge_num + String(" edges.") );
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
		std::vector< std::vector< std::pair<Size,Size> > >& aligned_pairs = pot_pairs; // indices of aligned_spectra in original experiment
		std::vector< std::vector<DoubleReal> > mod_positions(aligned_pairs.size()); // aligned_spectra (s' is a subset of peaks from s - the aligned ones )
		logger.startProgress(0,aligned_pairs.size(),"finding putative modification positions");
		edge_selection.clear();
		for(Size oi = 0; oi < aligned_pairs.size(); ++oi)
		{
			logger.setProgress(oi);
			for(Size i = 0; i < aligned_pairs[oi].size(); ++i)
			{
				MSSpectrum<Peak1D> res_1, res_2;
				DoubleReal mod_pos, score;
				asa.getAntisymetricAlignment(res_1, res_2, score, mod_pos, experiment[aligned_pairs[oi][i].first], experiment[aligned_pairs[oi][i].second]);
				/// @improvement make pairs_min_matches an advanced parameter
				/* debug std::cout << " with pair: "<< aligned_pairs[i].first << " + " << aligned_pairs[i].second << std::endl; */
				if((DoubleReal)res_1.size()/((DoubleReal)(experiment[aligned_pairs[oi][i].first].size()+experiment[aligned_pairs[oi][i].second].size())/DoubleReal(2)) >=  pairs_min_ratio/*  and res_1.size() >= pairs_min_matches */)
				{
					edge_selection.push_back(i);
					/// @improvement pot. speedup by removing MetaValues ?!
					res_1.setMetaValue("original index", aligned_pairs[oi][i].first);
					res_1.setMetaValue("paired with index", aligned_pairs[oi][i].second);
					res_2.setMetaValue("original index", aligned_pairs[oi][i].second);
					res_2.setMetaValue("paired with index", aligned_pairs[oi][i].first);
					aligned_spectra.push_back(res_1);
					aligned_spectra.push_back(res_2);
					mod_positions[oi].push_back(mod_pos);
				}
			}
			for(Size i = 0; i < edge_selection.size(); ++i)
			{
				//~ edge_selection indices are alwas >= the ones to aligned_pairs so no collision expected
				aligned_pairs[oi][i] = aligned_pairs[oi][edge_selection[i]];
			}
			aligned_pairs[oi].resize(edge_selection.size());
			edge_selection.clear();
		}

		logger.endProgress();

		//~ try{
		//~ String alignments(outputfile_name_edges + String(".alignmentspecs.mzML"));
		//~ mzml.store(alignments , aligned_spectra);
		//~ }catch(...){
			writeLog_(String(".. aligned spectra not separately saved .."));
		//~ }

		edge_num = 0;
		for(Size oi = 0; oi < aligned_pairs.size(); ++oi)
		{
			edge_num+=aligned_pairs[oi].size();
		}

		w.stop();
		writeLog_(String(".. done, took ") + String(w.getClockTime()) + String(" seconds, defined ") + edge_num + String(" edges."));
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
		for(Size oi = 0; oi < aligned_pairs.size(); ++oi)
		{
			for(Size i = 0; i < aligned_pairs[oi].size(); ++i)
			{
				indices_of_i_in_aligned_spectra[aligned_pairs[oi][i].first].insert(i*2);
				indices_of_i_in_aligned_spectra[aligned_pairs[oi][i].second].insert((i*2)+1);
			}
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
				experiment[it->first] = BinnedSpectrum<Peak1D>(consensus_sz, consensus_sp, 0.1, 0.01, unmerged);
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
			for(Size oi = 0; oi < aligned_pairs.size(); ++oi)
			{
				for(Size i = 0; i < aligned_pairs[oi].size(); ++i)
				{
					/// @improvement save also the edges relation (i.e. equal/modified/ladder)
					stars_txt.push_back(String(aligned_pairs[oi][i].first) +String('#') +String(aligned_pairs[oi][i].second) +String('#') +String(mod_positions[oi][i]) +String('\n'));
				}
			}
			stars_txt.push_back(String("\r\n"));
			stars_txt.store(outputfile_name_edges);
		}
		else
		{
			writeLog_(String("No edges."));
		}
		writeLog_(String("All done. Please proceed with SpectralNetworking after identification of the consensuses and merging the ids."));
		return;
	}
	///

	void propagateNetwork_(MSExperiment<Peak1D>& experiment, std::vector< std::pair<Size,Size> >& aligned_pairs, std::vector< DoubleReal > mod_positions, ConsensusMap& ids )
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
		sc.propagateNetwork(experiment, aligned_pairs, mod_positions, ids);

		/*------------/
		// score etc.
		/------------*/
		/// @improvement add featurehandles for those edges with fabs(pm_diff) < pm_tol and sufficient match and copy peptidehits if no id
		/// @improvement dont punish to hard if spectrum does cover only a part of/or more than the sequence - e.g maybe mod that shortened the pep
		/// @improvement txt output of peptide laddering sequences
		/// @improvement prepare featurehandles(metainfo?!) for colorization/different linestyle
		Size c(0);
		//~ ConsensusMap::Iterator zero_to_delete;
		for(ConsensusMap::Iterator ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
		{
			if(ids_it->getPeptideIdentifications().front().empty())
			{
				ids_it->getPeptideIdentifications().clear();
				++c;
			}
			//~ if(ids_it->getRT()==0 and ids_it->getMZ()==0)
			//~ {
				//~ zero_to_delete = ids_it;
			//~ }
		}
		writeLog_(String("nodes without id: ") + String(c));
		/// @improvement ahhrg!!! doeas not work - hardhack for the annoying 0,0 ConsensusFeature - find out where that comes from
		//~ ids.erase(zero_to_delete);

		std::vector<String> colors;
		//~ colors interpretable by QColor(QString) see http://www.w3.org/TR/SVG/types.html#ColorKeywords
		colors.push_back("blue				");
		colors.push_back("brown				");
		colors.push_back("cyan				");
		colors.push_back("gold				");
		colors.push_back("gray				");
		colors.push_back("grey				");
		colors.push_back("green				");
		colors.push_back("magenta			");
		colors.push_back("red					");
		colors.push_back("violet			");
		colors.push_back("darkblue		");
		colors.push_back("darkcyan		");
		colors.push_back("darkgray		");
		colors.push_back("darkgreen		");
		colors.push_back("darkgrey		");
		colors.push_back("darkmagenta	");
		colors.push_back("darkred			");
		colors.push_back("darkviolet	");

		std::map< Size, std::vector<String> > cropped_ids;

		for(ConsensusMap::Iterator ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
		{
			if(ids_it->metaValueExists("network id"))
			{
				ids_it->setMetaValue("color", colors[(Size)ids_it->getMetaValue("network id")%colors.size()]);
				//~ set intensity for coloration
				//~ ids_it->setIntensity(1000*((Size)ids_it->getMetaValue("network id"))%colors.size());
				cropped_ids[(Size)ids_it->getMetaValue("network id")].push_back(ids_it->getPeptideIdentifications().front().getHits().front().getSequence().toString());
			}
		}

		for(std::map< Size, std::vector<String> >::iterator it = cropped_ids.begin(); it != cropped_ids.end(); ++it)
		{
			for(Size i = 0; i < it->second.size(); ++i)
			{
				std::cout << it->second[i] << std::endl;
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

		String consensus_file = getStringOption_("consensus");
		String ids_file = getStringOption_("ids");
		String network_file = getStringOption_("network");

		StopWatch w;
		w.start();

		//-------------------------------------------------------------
		// read experiment from file, if a consensus
		// was given, it is preferred, else filtering
		// will be applied to unprocessed data
		//-------------------------------------------------------------
		MzMLFile mzml;
		//MSExperiment<Peak1D>= PeakMap see StandardTypes
		MSExperiment<Peak1D> experiment;
		if(!consensus_file.empty())
		{
			try{
				mzml.load(inputfile_name,experiment);
				writeLog_(String("Loading dataset from ") + consensus_file +  String(" .. ") );
			}catch(...){
				writeLog_(String("Invalid input, aborting!"));
				return PARSE_ERROR;
			}
			//~ experiment.sortSpectra();
		}
		else
		{
			try{
				mzml.load(inputfile_name,experiment);
				writeLog_(String("Loading dataset from ") + inputfile_name +  String(" .. ") );
			}catch(...){
				writeLog_(String("Invalid input, aborting!"));
				return PARSE_ERROR;
			}
			//~ experiment.sortSpectra();
		}

		w.stop();
		writeLog_(String("done  .. took ") + String(w.getClockTime()) + String(" seconds, loaded ") + String(experiment.size()) + String(" spectra"));
		w.reset();

		if(!ids_file.empty())
		{
			writeLog_(String("Loading identifications to dataset from ") + ids_file +  String(" ..") );

			IdXMLFile idxml;
			std::vector< ProteinIdentification > protein_ids;
			std::vector< PeptideIdentification > peptide_ids;
			try{
				idxml.load(ids_file, protein_ids, peptide_ids);
			}catch(...){
				writeLog_(String("Invalid input, aborting!"));
				return PARSE_ERROR;
			}

			IDMapper idm;
			idm.annotate(experiment, peptide_ids, protein_ids);

			StarClusters<Peak1D> scorer;
			for(Size i = 0; i < experiment.size(); ++i)
			{
				if(experiment[i].getPeptideIdentifications().size()>0)
				{
					scorer.idPrepare(experiment[i]);
					//~ take over the charge from the best DB hit
					experiment[i].getPrecursors().front().setCharge(experiment[i].getPeptideIdentifications().front().getHits().front().getCharge());
				}
				else
				{
					//~ by default a chargestate of 2 is assumed if not set
					if(experiment[i].getPrecursors().front().getCharge()<1)
					{
						/// @improvement apply a chargestate enumeration (1-2-3) in edgeselection and adopt best charge
						experiment[i].getPrecursors().front().setCharge(2);
					}
				}
			}

			writeLog_(String(" loaded ") + String(peptide_ids.size()) + String(" ids"));
		}

		//parent peak mower
		ParentPeakMower papemo;
		Param p(papemo.getParameters());
		p.setValue("remove_big_outliners", 1);
		papemo.setParameters(p);
		for(Size i = 0; i < experiment.size(); ++i)
		{
			papemo.filterSpectrum(experiment[i]);
			experiment[i].sortByPosition();
		}

		//for scoring reasons this is mandatory
		Normalizer normalizer;
		p = normalizer.getParameters();
		p.setValue("method", "to_TIC");
		normalizer.setParameters(p);
		for(Size i = 0; i < experiment.size(); ++i)
		{
			normalizer.filterSpectrum(experiment[i]);
			experiment[i].sortByPosition();
		}

		//threshold mower
		ThresholdMower thremowe;
		p = thremowe.getParameters();
		p.setValue("threshold", 0.01);
		thremowe.setParameters(p);
		for(Size i = 0; i < experiment.size(); ++i)
		{
			thremowe.filterSpectrum(experiment[i]);
			experiment[i].sortByPosition();
		}

		///@improvement find out if only the b0 ions or only the y0 ions suffice
		addZeroIonPeaks_(experiment);

		//~ ///@improvement find out if 'all_one' distorts the alignment
		//~ // normalizer to one
		//~ p = normalizer.getParameters();
		//~ p.setValue("method", "all_one");
		//~ normalizer.setParameters(p);
		//~ for(Size i = 0; i < experiment.size(); ++i)
		//~ {
			//~ normalizer.filterSpectrum(experiment[i]);
			//~ experiment[i].sortByPosition();
		//~ }

		if(network_file.empty())
		{
			String outputfile_name_consensus = dir+prefix+String("-consensus.mzML");
			outputFileWritable_(outputfile_name_consensus);
			String outputfile_name_edges = dir+prefix+String("-edges.txt");
			outputFileWritable_(outputfile_name_edges);

			//~ previous clustering of nearly equals
				///@attention if peak addition is applied this should be done _after_ clustering
			writeLog_(String("Cluster highly identical spectra .. ") );
			std::vector< std::set<Size> > clusters;
			clusterInTolerance_(experiment, clusters);
			std::set<Size> no_edges;
			for(Size i = 0; i < clusters.size(); ++i)
			{
				std::set<Size>::iterator it = clusters[i].begin();
				++it;
				for(it ; it != clusters[i].end(); ++it)
				{
					no_edges.insert(*it);
				}
			}
			writeLog_(String("Cut down ") + String(no_edges.size()) + String(" spectra, ") + String(experiment.size()-no_edges.size()) + String(" remain") );

			//~ start consensus making
			spanNetwork_(experiment, outputfile_name_consensus, outputfile_name_edges, no_edges);
		}
		else if(!ids_file.empty() and !network_file.empty())
		{
			//make propagation case
			inputFileReadable_(ids_file);
			inputFileReadable_(network_file);

			String outputfile_name_specnet = dir+prefix+String("-specnet.consensusXML");
			outputFileWritable_(outputfile_name_specnet);

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
					if(!splits.empty())
					{
						writeLog_(String("Invalid input, aborting!"));
						return PARSE_ERROR;
					}
					else
					{
						writeLog_(String("Ooops file contains empty lines, omitting!"));
						continue;
					}
				}
				try
				{
					aligned_pairs.push_back(std::make_pair<Size,Size>(splits[0].toDouble(),splits[1].toDouble()));
					mod_positions.push_back(splits[2].toDouble());
					/*debug if(splits[2].toDouble()>0){std::cout<<splits[2].toDouble()<<std::endl;}*/
				}
				catch(...)
				{
					writeLog_(String("Erroneous file for input, aborting!"));
					return PARSE_ERROR;
				}
			}

			w.stop();
			writeLog_(String("done  .. took ") + String(w.getClockTime()) + String(" seconds, loaded ") + String(aligned_pairs.size()) + String(" edges"));
			w.reset();

			w.start();

			//start propagation
			ConsensusMap ids;
			propagateNetwork_(experiment, aligned_pairs, mod_positions, ids);
			//~ ids.applyMemberFunction(&UniqueIdInterface::setUniqueId);
			ConsensusXMLFile cxml;
			writeLog_(String("Saving propagated network ..") );
			cxml.store(outputfile_name_specnet, ids);

			w.stop();
			writeLog_(String("All done. Saving took ") + String(w.getClockTime()) + String(" seconds. Ready to investigate in TOPPView."));
			w.reset();
		}
		else
		{
			//invalid input case
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
