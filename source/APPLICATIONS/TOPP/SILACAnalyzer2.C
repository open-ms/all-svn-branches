// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: $
// --------------------------------------------------------------------------

//OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

//clustering
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/COMPARISON/CLUSTERING/HashClustering.h>
#include <OpenMS/COMPARISON/CLUSTERING/CentroidLinkage.h>

//Contrib includes
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>

using namespace OpenMS;

typedef std::vector<BinaryTreeNode> Tree;
typedef std::vector<DataPoint*> Cluster;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SILACAnalyzer2 SILACAnalyzer2

  @brief Identifies peptide pairs in LC-MS data and determines their relative abundance.

  SILACAnalyzer2 is a tool for the fully automated analysis of quantitative proteomics data. It identifies pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we first explain the algorithm and then discuss the tuning of its parameters.

  <b>Algorithm</b>

  The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).

  We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
  - all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold and
  - the intensities within the first envelope (at m/z, m/z+0.5 and m/z+1) and second envelope (at m/z+3, m/z+3.5 and m/z+4) decrease successively.

  Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster centers around the unlabelled peptide of a pair. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

  @image html SILACAnalyzer_algorithm.png

  <b>Parameter Tuning</b>

  SILACAnalyzer2 can either search for SILAC pairs (as described in the above paragraph) or triplets:
  - type - double (for SILAC pairs), triple (for SILAC triplets)

  <i>input:</i>
  - in [*.mzML] - LC-MS dataset to be analyzed
  - ini [*.ini] - file containing all parameters (see discussion below)

  <i>standard output:</i>
  - out [*.consensusXML] - contains the list of identified peptide pairs (retention time and m/z of the lighter peptide, heavy-to-light ratio)
  - out_visual [*.featureXML] - contains the complete set of data points (retention time, m/z, intensity) of all peptide pairs

  The results of an analysis can easily visualized within TOPPView. Simply load *.consensusXML and *.featureXML as layers over the original *.mzML.

  <i>optional output:</i>
  @n If -silac_debug is enabled, SILACAnalyzer2 generates a number of further files:
  - [*.dat] - contains the list of identified peptide pairs in a simple text file, c.f. *.consensusXML
  - [*_clusters.dat] -  contains the complete set of data points of all peptide pairs in a simple text file, c.f. *.featureXML
  - [*.input] - gnuplot script for the visualization of the results. Running (gnuplot *.input) generates a number of *.eps plots. The range of clusters to be plotted can be specified by the parameters -cluster_min/max.

  The following parameters are straightforward:
  - mass_separation_light_medium - mass gap between light and medium isotopic envelopes [Da] (only relevant for the search for SILAC triplets, i.e. type triple)
  - mass_separation_light_heavy - mass gap between light and heavy isotopic envelopes [Da]
  - charge_min/max - range of charge states
  - mz_stepwidth - step width with which the interpolated spectrum, Fig. (b), is scanned. The step width should be of about the same order with which the raw data were recorded, see Fig. (a).

  The remaining parameters should be tuned in the following order:
  - intensity_cutoff - adjust the intensity cutoff such that the data points that pass the non-local filter (*.featureXML layer) form clear distinct clusters,  see Fig. (e). Ignore the coloring of the clusters at that stage.
  - rt_scaling - pick a representative cluster. rt_scaling = (width of the cluster in Da)/(height of the cluster in sec)
  - cluster_number_scaling - The clustering algorithm tries to determine the optimal number of clusters (i.e. the number of peptide pairs in the LC-MS data set). If neighboring clusters appear in the same color, the cluster number is too low. If a single cluster contains two colors, the cluster number is too high. The cluster number can be adjusted by this scaling factor.
  - optimal_silhouette_tolerance - The clustering algorithm tries to maximize the average-silhouette-width, see details in reference. The parameter specifies the relative tolerance (in %) by which the optimum can deviate from the maximum.

  <b>References:</b>
  @n L. Nilse, M. Sturm, D. Trudgian, M. Salek, P. Sims, K. Carroll, S. Hubbard, "SILACAnalyzer2 - a tool for differential quantitation of stable isotope derived data", unpublished.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SILACAnalyzer.cli
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

bool clusterCmp( DataPoint a, DataPoint b ) {
		if ( a.cluster_size == b.cluster_size && a.cluster_id > b.cluster_id) return true;
		if ( a.cluster_size > b.cluster_size ) return true;
		return false;
	}

class TOPPSILACAnalyzer2
: public TOPPBase
{
private:

	String type;
	DoubleReal mass_separation_light_medium;
	DoubleReal mass_separation_light_heavy;
	UInt charge_min;
	UInt charge_max;
	DoubleReal mz_stepwidth;
	DoubleReal intensity_cutoff;
	DoubleReal mz_threshold;
	DoubleReal rt_threshold;
	DoubleReal rt_scaling;
	DoubleReal optimal_silhouette_tolerance;
	DoubleReal cluster_number_scaling;
	int cluster_min;
	int cluster_max;
	String in;
	String out;
	String out_visual;
	ConsensusMap all_pairs;
	FeatureMap<> all_cluster_points;


	void handleParameters()
	{
		type = getStringOption_("type");
		if (type=="double") {
			mass_separation_light_medium = getParam_().getValue("algorithm:mass_separation_light_heavy");
		}
		else {
			mass_separation_light_medium = getParam_().getValue("algorithm:mass_separation_light_medium");
		}
		mass_separation_light_heavy = getParam_().getValue("algorithm:mass_separation_light_heavy");
		charge_min = getParam_().getValue("algorithm:charge_min");
		charge_max = getParam_().getValue("algorithm:charge_max");
		intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
		mz_threshold = getParam_().getValue("algorithm:mz_threshold");
		rt_threshold = getParam_().getValue("algorithm:rt_threshold");
		rt_scaling = getParam_().getValue("algorithm:rt_scaling");
		optimal_silhouette_tolerance = getParam_().getValue("algorithm:optimal_silhouette_tolerance");
		cluster_number_scaling = getParam_().getValue("algorithm:cluster_number_scaling");
		cluster_min = getParam_().getValue("algorithm:cluster_min");
		cluster_max = getParam_().getValue("algorithm:cluster_max");
		in = getStringOption_("in");
		out = getStringOption_("out");
		out_visual = getStringOption_("out_visual");
		//output variables
		all_pairs.getFileDescriptions()[0].filename = in;
		all_pairs.getFileDescriptions()[0].label = "light";
		all_pairs.getFileDescriptions()[1].filename = in;
		all_pairs.getFileDescriptions()[1].label = "medium";
		all_pairs.getFileDescriptions()[2].filename = in;
		all_pairs.getFileDescriptions()[2].label = "heavy";
		all_pairs.setExperimentType("silac");

	}



	std::vector<DataPoint> buildDataStructure(MSExperiment<Peak1D>& exp, int charge)
	{
		ProgressLogger logger_;
		std::vector<DataPoint> data;
		int id=0;

		logger_.setLogType(log_type_);
		logger_.startProgress(0,exp.size(),"reducing raw data");
		DoubleReal envelope_distance_light_medium = mass_separation_light_medium / (DoubleReal)charge;
		DoubleReal envelope_distance_light_heavy = mass_separation_light_heavy / (DoubleReal)charge;
		DoubleReal isotope_distance = 1.0 / (DoubleReal)charge;

		//Extract level 1 spectra
		IntList levels=IntList::create("1");
		exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<MSExperiment<Peak1D>::SpectrumType>(levels, true)), exp.end());

		//Estimate m/z step width
		UInt i=0;
		while(i<exp.size() && exp[i].size()<5) ++i;
		std::vector<Real> mz_spacing;
		std::vector<Real> rt_spacing;
		for (Size j=1; j<exp[i].size(); ++j)
		{
			mz_spacing.push_back(exp[i][j].getMZ()-exp[i][j-1].getMZ());
//			rt_spacing.push_back(exp[i][j].getRT()-exp[i][j-1].getRT());
		}
		std::sort(mz_spacing.begin(),mz_spacing.end());
		mz_stepwidth=mz_spacing[mz_spacing.size()/2];
//		std::sort(rt_spacing.begin(),rt_spacing.end());
//		rt_stepwidth=rt_spacing[rt_spacing.size()/2];

		// scan over the entire experiment and write to data structure
		for (MSExperiment<>::Iterator rt_it=exp.begin(); rt_it!=exp.end(); ++rt_it)
		{
			logger_.setProgress(rt_it-exp.begin());
			Size number_data_points = rt_it->size();
			// spectra with less than 10 data points are being ignored
			if (number_data_points>=10) { //filter MS1 spectra (
				// read one OpenMS spectrum into GSL structure
				std::vector<DoubleReal> mz_vec;
				std::vector<DoubleReal> intensity_vec;
				mz_vec.resize(number_data_points);
				intensity_vec.resize(number_data_points);
				Int j = 0;
				for (MSSpectrum<>::Iterator mz_it=rt_it->begin(); mz_it!=rt_it->end(); ++mz_it)
				{
					mz_vec[j] = mz_it->getMZ();
					intensity_vec[j] = mz_it->getIntensity();
					++j;
				}
				DoubleReal mz_min = mz_vec[0];
				DoubleReal mz_max = mz_vec[number_data_points-1];
				// linear interpolation
				// used for the detection of pairs (spline overestimates at noise level)
				gsl_interp_accel *acc = gsl_interp_accel_alloc();
				gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, number_data_points);
				gsl_spline_init(spline, &*mz_vec.begin(), &*intensity_vec.begin(), number_data_points);
				// spline interpolation
				// used for exact ratio calculation (more accurate when real peak pairs are present)
				gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
				gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, number_data_points);
				gsl_spline_init(spline2, &*mz_vec.begin(), &*intensity_vec.begin(), number_data_points);
				for (DoubleReal mz=mz_min+isotope_distance; mz<mz_max-envelope_distance_light_heavy-3*isotope_distance; mz+=mz_stepwidth)
				{
					DoubleReal int_lin1 = gsl_spline_eval (spline, mz, acc);
					DoubleReal int_lin2 = gsl_spline_eval (spline, mz+isotope_distance, acc);
					DoubleReal int_lin3 = gsl_spline_eval (spline, mz+2*isotope_distance, acc);
					DoubleReal int_lin4 = gsl_spline_eval (spline, mz+envelope_distance_light_medium, acc);
					DoubleReal int_lin5 = gsl_spline_eval (spline, mz+envelope_distance_light_medium+isotope_distance, acc);
					DoubleReal int_lin6 = gsl_spline_eval (spline, mz+envelope_distance_light_medium+2*isotope_distance, acc);
					DoubleReal int_lin7 = gsl_spline_eval (spline, mz+envelope_distance_light_heavy, acc);
					DoubleReal int_lin8 = gsl_spline_eval (spline, mz+envelope_distance_light_heavy+isotope_distance, acc);
					DoubleReal int_lin9 = gsl_spline_eval (spline, mz+envelope_distance_light_heavy+2*isotope_distance, acc);
					DoubleReal int_spline1 = gsl_spline_eval (spline2, mz, acc2);
					DoubleReal int_spline2 = gsl_spline_eval (spline2, mz+isotope_distance, acc2);
					DoubleReal int_spline3 = gsl_spline_eval (spline2, mz+2*isotope_distance, acc2);
					DoubleReal int_spline4 = gsl_spline_eval (spline2, mz+envelope_distance_light_medium, acc2);
					DoubleReal int_spline5 = gsl_spline_eval (spline2, mz+envelope_distance_light_medium+isotope_distance, acc2);
					DoubleReal int_spline6 = gsl_spline_eval (spline2, mz+envelope_distance_light_medium+2*isotope_distance, acc2);
					DoubleReal int_spline7 = gsl_spline_eval (spline2, mz+envelope_distance_light_heavy, acc2);
					DoubleReal int_spline8 = gsl_spline_eval (spline2, mz+envelope_distance_light_heavy+isotope_distance, acc2);
					DoubleReal int_spline9 = gsl_spline_eval (spline2, mz+envelope_distance_light_heavy+2*isotope_distance, acc2);

					bool condDouble1 = (int_lin1 >= intensity_cutoff) && (int_lin2 >= intensity_cutoff) && (int_lin3 >= intensity_cutoff) && (int_lin7 >= intensity_cutoff) && (int_lin8 >= intensity_cutoff) && (int_lin9 >= intensity_cutoff); // all six intensities peak simultaneously
					bool condDouble2 = (int_spline1 >= int_spline2) && (int_spline2 >= int_spline3) && (int_spline7 >= int_spline8) && (int_spline8 >= int_spline9); // isotopic peaks within one envelop decrease
					bool condTriple1 = (int_lin1 >= intensity_cutoff) && (int_lin2 >= intensity_cutoff) && (int_lin3 >= intensity_cutoff) && (int_lin4 >= intensity_cutoff) && (int_lin5 >= intensity_cutoff) && (int_lin6 >= intensity_cutoff) && (int_lin7 >= intensity_cutoff) && (int_lin8 >= intensity_cutoff) && (int_lin9 >= intensity_cutoff); // all nine intensities peak simultaneously
					bool condTriple2 = (int_spline1 >= int_spline2) && (int_spline2 >= int_spline3) && (int_spline4 >= int_spline5) && (int_spline5 >= int_spline6) && (int_spline7 >= int_spline8) && (int_spline8 >= int_spline9); // isotopic peaks within one envelop decrease
					if ((type=="double" && condDouble1 && condDouble2) || (type=="triple" && condTriple1 && condTriple2))
					{
						data.push_back(DataPoint(rt_it->getRT(),mz,int_spline1,id));
//						data.push_back(DataPoint(id,rt_it->getRT(),mz,int_spline1,int_spline2,int_spline3,int_spline4,int_spline5,int_spline6,int_spline7,int_spline8,int_spline9));
						++id;
					}
				}

				gsl_spline_free(spline);
				gsl_interp_accel_free(acc);
				gsl_spline_free(spline2);
				gsl_interp_accel_free(acc2);
			}
		}
		exp.clear(true);
		logger_.endProgress();
		return data;
	}



public:
	TOPPSILACAnalyzer2()
	: TOPPBase("SILACAnalyzer2","Determination of peak ratios in LC-MS data",true)
	{
	}

	void registerOptionsAndFlags_()
	{
		registerStringOption_("type","<name>","","SILAC experiment type\n",true);
		setValidStrings_("type", getToolList()[toolName_()] );

		registerInputFile_("in","<file>","","input file");
		setValidFormats_("in",StringList::create("mzML"));
		registerOutputFile_("out","<file>","","output file", false);
		setValidFormats_("out",StringList::create("consensusXML"));
		registerOutputFile_("out_visual","<file>","","output file containing cluster information",false);
		setValidFormats_("out_visual",StringList::create("featureXML"));

		registerFlag_("silac_debug","Enables writing of debug information",true);
		registerSubsection_("algorithm","Algorithm parameters section");

	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		Int SILAC_type = (type=="double" ?  2 : 3 );
		Param tmp;

		if (SILAC_type == 2)
		{
			tmp.setValue("mass_separation_light_heavy",6.0202, "mass gap between light and heavy isotopic envelopes, [Da]" );
		}
		else if (SILAC_type == 3)
		{
			tmp.setValue("mass_separation_light_medium" , 6.0202, "mass gap between light and medium isotopic envelopes, [Da]");
			tmp.setValue("mass_separation_light_heavy" , 8.0202, "mass gap between light and heavy isotopic envelopes, [Da]");
		}

		tmp.setValue("charge_min", 2, "Charge state range begin");
		tmp.setMinInt("charge_min", 1);

		tmp.setValue("charge_max", 3, "Charge state range end");
		tmp.setMinInt("charge_max",1);

		tmp.setValue("intensity_cutoff", 5000.0, "intensity cutoff");
		tmp.setMinFloat("intensity_cutoff", 0.0);

		tmp.setValue("mz_threshold", 50.0, "mz_threshold");
		tmp.setMinFloat("mz_threshold", 0.0);

		tmp.setValue("rt_threshold", 50.0, "rt_threshold");
		tmp.setMinFloat("rt_threshold", 0.0);

		tmp.setValue("rt_scaling",0.05,"scaling factor of retention times (Cluster height [s] an\ncluster width [Th] should be of the same order. The clustering algorithms work better for\nsymmetric clusters.)");
		tmp.setMinFloat("rt_scaling", 0.0);

		tmp.setValue("optimal_silhouette_tolerance",0.0,"The partition with most clusters is chosen, which deviates from the optimal silhouette width at most by this percentage.");
		tmp.setMinFloat("optimal_silhouette_tolerance",0.0);
		tmp.setMaxFloat("optimal_silhouette_tolerance",100.0);

		tmp.setValue("cluster_number_scaling", 1.0, "scaling factor of the number of clusters (The average-silhouette-width\nalgorithm returns an 'optimal' number of clusters. This number might need\nto be adjusted by this factor.)");
		tmp.setMinFloat("cluster_number_scaling", 0.0);

		tmp.setValue("cluster_min", 0, "Start of the clusters range to be plotted by the gnuplot script");
		tmp.setMinInt("cluster_min", 0);

		tmp.setValue("cluster_max", 2, "End of the clusters range to be plotted by the gnuplot script");
		tmp.setMinInt("cluster_max", 0);


		return tmp;
	}



	ExitCodes main_(int , const char**)
	{
		handleParameters();
		//--------------------------------------------------------------
		// determine file name for debug output
		//--------------------------------------------------------------
		String debug_trunk = in;
		if (in.has('.'))
		{
			debug_trunk = in.substr(0,-SignedSize(in.suffix('.').length())-1);
		}

		//iterate over all charge states
		for (UInt charge=charge_min; charge<=charge_max; ++charge)
		{
			std::cout << "charge state: " << charge << "+" << std::endl;


			// For each charge state run the experimental data (exp) are loaded again. Either the raw data (exp) or the distance matrix (distance_matrix) are in memory which keeps the memory footprint low.

			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			MzMLFile file;
			MSExperiment<Peak1D> exp;

			file.setLogType(log_type_);
			file.load(in,exp);

			//set input map size (only once)
			if (charge==charge_min)
			{
				exp.updateRanges();
				all_pairs.getFileDescriptions()[0].size = exp.getSize();
				all_pairs.getFileDescriptions()[1].size = exp.getSize();
				all_pairs.getFileDescriptions()[2].size = exp.getSize();
			}

			//-------------------------------------------------------------
			// build SILACData structure
			//-------------------------------------------------------------

			std::vector<DataPoint> data=buildDataStructure(exp, charge);

			//-------------------------------------------------------------
			// Perform clustering
			//-------------------------------------------------------------

			CentroidLinkage method(rt_scaling);
			HashClustering c(data,rt_threshold,mz_threshold,method);
			c.setLogType(log_type_);
			std::vector<Tree> subtrees;
			c.performClustering(subtrees);
			std::vector<std::vector<Real> > silhouettes;
			Size cluster_id=0;

			//run silhoutte optimization for all subtrees and find the appropriate best_n
			for ( std::vector<Tree>::iterator it=subtrees.begin();it!=subtrees.end();++it)
			{
				Tree tree=*it;
				std::vector< Real > asw = c.averageSilhouetteWidth(tree);
				silhouettes.push_back(asw);
				//Look only in the front area of the silhoutte values to avoid getting the wrong number
				std::vector< Real >::iterator max_el(max_element((asw.end()-((int)it->size()/10) ),asw.end()));
				Size best_n = (Size)tree.size();
				//Silhouette method can not create a single cluster from one subtree. If silhouette values in the front area are very low, take the subtree as one cluster
				if (*max_el< 0.7)
					best_n=1;
				else
				{
					for (Size i = 0; i < asw.size(); ++i)
					{
						if(std::fabs(asw[i]-(*max_el))<=0)
						{
							best_n = Size(it->size() - i);
							break;
						}
					}
				}

				best_n = Size(cluster_number_scaling * best_n); // slightly increase cluster number
				std::vector<Cluster> act_clusters;
				//Get clusters from the subtree
				c.cut(best_n,act_clusters,tree);
				//Run through all elements in each cluster and update cluster number
				for (std::vector<Cluster>::iterator cluster_it=act_clusters.begin();cluster_it!=act_clusters.end();++cluster_it)
				{
					Int rt=0;
					Int mz=0;
					Int intensity=0;
					for (Cluster::iterator element_it=cluster_it->begin();element_it!=cluster_it->end();++element_it)
					{
						//Calculate center of consensus feature for output
						if (out!="")
						{
							mz+=(*element_it)->mz;
							rt+=(*element_it)->rt;
							intensity+=(*element_it)->intensity;
						}
						(*element_it)->cluster_id=cluster_id;
						(*element_it)->cluster_size=cluster_it->size();
					}
					//Prepare output
					if (out!="")
					{
						mz/=cluster_it->size();
						rt/=cluster_it->size();
						intensity/=cluster_it->size();
						ConsensusFeature pair_light_heavy;
						pair_light_heavy.setRT(rt);
						pair_light_heavy.setMZ(mz);
						pair_light_heavy.setIntensity(intensity);
						pair_light_heavy.setCharge(charge);
						pair_light_heavy.setQuality(*max_el);
						FeatureHandle handle;
						handle.setRT(rt);
						handle.setMZ(mz);
						handle.setIntensity(intensity);
						handle.setCharge(charge);
						handle.setMapIndex(0);
						handle.setUniqueId(cluster_id);
						pair_light_heavy.insert(handle);
						handle.setRT(rt);
						DoubleReal envelope_distance_light_heavy = mass_separation_light_heavy / (DoubleReal)charge;
						handle.setMZ(mz+envelope_distance_light_heavy);
						handle.setIntensity(intensity);
						handle.setCharge(charge);
						handle.setMapIndex(2);
						handle.setUniqueId(cluster_id);
						pair_light_heavy.insert(handle);
						all_pairs.push_back(pair_light_heavy);
					}
					++cluster_id;
				}

			}
			//Sort data points: High cluster numbers first
			std::sort(data.begin(),data.end(),clusterCmp);

			//--------------------------------------------------------------
			//create features (for visualization)
			//--------------------------------------------------------------
			if (out_visual!="")
			{
				std::vector<String> colors;
				// 16 HTML colors
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

				for (std::vector<DataPoint>::iterator it=data.begin(); it!= data.end(); ++it)
				{
					//visualize the light variant
					Feature cluster_point;
					cluster_point.setRT(it->rt);
					cluster_point.setMZ(it->mz);
					cluster_point.setIntensity(it->intensity);
					cluster_point.setCharge(charge);
					cluster_point.setMetaValue("cluster_id",it->cluster_id);
					cluster_point.setMetaValue("color",colors[it->cluster_id%colors.size()]);
					all_cluster_points.push_back(cluster_point);
				}
				// required, as somehow the order of features on some datasets between Win & Linux is different and thus the TOPPtest might fail
				all_cluster_points.sortByPosition();
			}





//			//-------------------------------------------------------------
//			// generate debug output
//			//-------------------------------------------------------------
//			// strings repeatedly used in debug output
			String light_medium_string = String(0.01*floor(mass_separation_light_medium*100+0.5));
			String light_heavy_string = String(0.01*floor(mass_separation_light_heavy*100+0.5));


			if (getFlag_("silac_debug"))
			{
				String debug_suffix;
				if (type=="double") {
					debug_suffix = "_" + light_heavy_string + "Da_" + String(charge) +"+";
				}
				else {
					debug_suffix = "_" + light_medium_string + "Da_" + light_heavy_string + "Da_" + String(charge) +"+";
				}
				// names of dat files
				String debug_clusters_dat = debug_trunk + debug_suffix + "_clusters.dat";

				// write all cluster data points to *_clusters.dat

				std::ofstream stream_clusters(debug_clusters_dat.c_str());
				stream_clusters << "cluster_id cluster_size rt mz int" << std::endl;
				Int current_id = -1;
				for (std::vector<DataPoint>::iterator it=data.begin(); it!= data.end(); ++it)
				{
					if (it->cluster_id != current_id) stream_clusters << std::endl << std::endl;
					stream_clusters << it->cluster_id << " " << it->cluster_size << " "<< it->rt << " " << it->mz << " " << it->intensity << std::endl;
					current_id = it->cluster_id;
				}
				stream_clusters.close();

				String debug_silhouettes_dat = debug_trunk + debug_suffix + "_silhouettes.dat";
				std::ofstream stream_silhouettes(debug_silhouettes_dat.c_str());

				for (std::vector<std::vector<Real> >::iterator asw_vector_it=silhouettes.begin();asw_vector_it!=silhouettes.end();++asw_vector_it)
				{
					for (std::vector<Real>::iterator asw_it=asw_vector_it->begin();asw_it!=asw_vector_it->end();++asw_it)
					{
						stream_silhouettes << *asw_it << "\t";
					}
					stream_silhouettes << std::endl;
				}
				stream_silhouettes.close();

				String debug_silhouettes_r = debug_trunk + debug_suffix + "_silhouettes.R";
				std::ofstream stream_r(debug_silhouettes_r.c_str());
				stream_r << "con <- file(\"" << debug_silhouettes_dat << "\",\"r\")" << std::endl;
				stream_r << "lines <- readLines(con, n=-1)" << std::endl;
				String pdf_name= debug_trunk + debug_suffix + "_silhouettes.pdf";
				String pdf_title= debug_trunk + debug_suffix;
				stream_r << "pdf(\""<< pdf_name << "\", paper=\"special\",width=13,height=(2*length(lines)),title = \"Silhoutte Plots for "<< pdf_title << "\")" << std::endl;
				stream_r << "par(mfrow=c(ceiling(length(lines)/2),2))" << std::endl;
				stream_r << "for (i in 1:length(lines))" << std::endl;
				stream_r << "{" << std::endl;
				stream_r << "text_vec <- strsplit(lines[i],\"\\t\")" << std::endl;
				stream_r << "asw <- as.numeric(text_vec[[1]])" << std::endl;
				stream_r << "aswrange <- asw[(length(asw)):(0)]" << std::endl;
				stream_r << "plot (1:length(asw),aswrange,type=\"l\",xlab=\"# Cluster\",ylab=\"asw\",main=paste(\"Subtree \",i,\": max n \",which.max(aswrange)))" << std::endl;
				stream_r << "}" << std::endl;
				stream_r << "close(con)" << std::endl;
				stream_r << "dev.off()" << std::endl;
				stream_r.close();
//
//				// write ratios of all cluster to *.dat
//				std::ofstream stream_ratios(debug_dat.c_str());
//				if (type=="double") {
//					stream_ratios << "cluster_id cluster_size rt mz ratio_light_heavy intensity" << std::endl;
//				}
//				else {
//					stream_ratios << "cluster_id cluster_size rt mz ratio_light_medium ratio_light_heavy intensity" << std::endl;
//				}
//				for (Size i=0; i<best_clusters.size();++i)
//				{
//					DoubleReal rt = 0.0;
//					DoubleReal mz = 0.0;
//					DoubleReal int_l = 0.0;
//					DoubleReal int_m = 0.0;
//					DoubleReal int_h = 0.0;
//					// intensity vectors used for linear regression
//					std::vector<DoubleReal> i1(3*cluster_size[i]);
//					std::vector<DoubleReal> i2(3*cluster_size[i]);
//					std::vector<DoubleReal> i3(3*cluster_size[i]);
//					UInt j=0;
//					DoubleReal isotope_distance = 1.0 / (DoubleReal)charge;
//					for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
//					{
//						if ((Size)it->cluster_id==i)
//						{
//							i1[3*j] = it->int1;
//							i2[3*j] = it->int4;
//							i3[3*j] = it->int7;
//							i1[3*j+1] = it->int2;
//							i2[3*j+1] = it->int5;
//							i3[3*j+1] = it->int8;
//							i1[3*j+2] = it->int3;
//							i2[3*j+2] = it->int6;
//							i3[3*j+2] = it->int9;
//
//							rt += it->rt;
//							if (it->int1>int_l)
//							{
//								int_l = it->int1;
//								mz = it->mz;
//							}
//							if (it->int2>int_l)
//							{
//								int_l = it->int2;
//								mz = it->mz + isotope_distance;
//							}
//							if (it->int3>int_l)
//							{
//								int_l = it->int3;
//								mz = it->mz + 2.0 * isotope_distance;
//							}
//							if (it->int4>int_m)
//							{
//								int_m = it->int4;
//							}
//							if (it->int5>int_m)
//							{
//								int_m = it->int5;
//							}
//							if (it->int6>int_m)
//							{
//								int_m = it->int6;
//							}
//							if (it->int7>int_h)
//							{
//								int_h = it->int7;
//							}
//							if (it->int8>int_h)
//							{
//								int_h = it->int8;
//							}
//							if (it->int9>int_h)
//							{
//								int_h = it->int9;
//							}
//							++j;
//						}
//					}
//					rt /= (DoubleReal)(cluster_size[i]); // average retention time
//					Math::LinearRegression linear_reg_light_medium;
//					Math::LinearRegression linear_reg_light_heavy;
//					linear_reg_light_medium.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
//					linear_reg_light_heavy.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i3.begin());
//					if (type=="double") {
//						stream_ratios << i << " " << cluster_size[i] << " " << rt << " " << mz << " " << linear_reg_light_heavy.getSlope() << " " << *max_element(i1.begin(),i1.end()) + *max_element(i2.begin(),i2.end()) << std::endl;
//					}
//					else {
//						stream_ratios << i << " " << cluster_size[i] << " " << rt << " " << mz << " " << linear_reg_light_medium.getSlope() << " " << linear_reg_light_heavy.getSlope() << " " << *max_element(i1.begin(),i1.end()) + *max_element(i2.begin(),i2.end()) << std::endl;
//					}
//				}
//				stream_ratios.close();
			}

		} //end iterate over all charge states


		//--------------------------------------------------------------
		//Store output
		//--------------------------------------------------------------
		if (out!="")
		{

			// assign unique ids
			all_pairs.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			//annotate output with data processing info
			addDataProcessing_(all_pairs, getProcessingInfo_(DataProcessing::QUANTITATION));

			ConsensusXMLFile c_file;
			c_file.store(out,all_pairs);
		}

		if (out_visual!="")
		{
			// assign unique ids
			all_cluster_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			FeatureXMLFile f_file;
			f_file.store(out_visual,all_cluster_points);
		}
//
//		//--------------------------------------------------------------
//		//write gnuplot script
//		//--------------------------------------------------------------
//		// + silhoutte width in gnuplot
//		// strings repeatedly used in debug output
//		String light_medium_string = String(0.01*floor(mass_separation_light_medium*100+0.5));
//		String light_heavy_string = String(0.01*floor(mass_separation_light_heavy*100+0.5));
//		String rt_scaling_string = String(0.01*floor(rt_scaling*100+0.5));
//		String optimal_silhouette_tolerance_string = String(0.01*floor(optimal_silhouette_tolerance*100+0.5));
//		String cluster_number_scaling_string = String(0.01*floor(cluster_number_scaling*100+0.5));
//
//		if (getFlag_("silac_debug"))
//		{
//			// first lines of the gnuplot script
//			String r_script = debug_trunk + ".input";
//			std::ofstream stream_r_script(r_script.c_str());
//			stream_r_script << "dat <- read.table(file=\""<<debug_silhouettes_dat<<"\")" << std::endl;
//			stream_r_script << "for (k in 0:dim(dat)[1])" << std::endl;
//			stream_gnuplotscript << "set size 2.0, 2.0" << std::endl;
//			stream_gnuplotscript << "set size square" << std::endl << std::endl;
//			//iterate over all charge states
//			for (Size charge=charge_min; charge<=charge_max; ++charge)
//			{
//				String debug_suffix;
//				if (type=="double") {
//					debug_suffix = "_" + light_heavy_string + "Da_" + String(charge) +"+";
//				}
//				else {
//					debug_suffix = "_" + light_medium_string + "Da_" + light_heavy_string + "Da_" + String(charge) +"+";
//				}
//				// names of dat files
//				String debug_dat = debug_trunk + debug_suffix + ".dat";
//				String debug_clusters_dat = debug_trunk + debug_suffix + "_clusters.dat";
//				// names of postscript files
//				String debug_ratios_light_medium_intensity = debug_trunk + debug_suffix + "_ratios_light_medium_intensity.eps";
//				String debug_ratios_light_medium = debug_trunk + debug_suffix + "_ratios_light_medium.eps";
//				String debug_ratios_light_heavy = debug_trunk + debug_suffix + "_ratios_light_heavy.eps";
//				String debug_sizes = debug_trunk + debug_suffix + "_sizes.eps";
//				String debug_clusters = debug_trunk + debug_suffix + "_clusters.eps";
//				String debug_clustersIntLightMedium = debug_trunk + debug_suffix + "_clustersIntLightMedium.eps";
//				String debug_clustersIntLightHeavy = debug_trunk + debug_suffix + "_clustersIntLightHeavy.eps";
//
//				// write *_clusters.eps
//				stream_gnuplotscript << "set output \"" + debug_clusters + "\"" << std::endl;
//				if (type=="double") {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				else {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				stream_gnuplotscript << "set xlabel \'m/Z (Th)\'" << std::endl;
//				stream_gnuplotscript << "set ylabel \'RT (s)\'" << std::endl;
//				stream_gnuplotscript << "plot";
//				for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
//				{
//					if (i != cluster_min) stream_gnuplotscript << ",";
//					stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 4:3 title \"cluster " + String(i) + "\"";
//				}
//				stream_gnuplotscript << std::endl;
//
//				// write *_clustersIntLightMedium.eps
//				if (type=="triple") {
//					stream_gnuplotscript << "set output \"" + debug_clustersIntLightMedium + "\"" << std::endl;
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//					stream_gnuplotscript << "set xlabel \'intensity at m/Z\'" << std::endl;
//					stream_gnuplotscript << "set ylabel \'intensity at m/Z + " + light_medium_string + "Th\'" << std::endl;
//					stream_gnuplotscript << "plot";
//					for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
//					{
//						if (i != cluster_min) stream_gnuplotscript << ",";
//						stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:8 title \"cluster " + String(i) + "\"";
//					}
//					stream_gnuplotscript << std::endl;
//				}
//
//				// write *_clustersIntLightHeavy.eps
//				stream_gnuplotscript << "set output \"" + debug_clustersIntLightHeavy + "\"" << std::endl;
//				if (type=="double") {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				else {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				stream_gnuplotscript << "set xlabel \'intensity at m/Z\'" << std::endl;
//				stream_gnuplotscript << "set ylabel \'intensity at m/Z + " + light_heavy_string + "Th\'" << std::endl;
//				stream_gnuplotscript << "plot";
//				for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
//				{
//					if (i != cluster_min) stream_gnuplotscript << ",";
//					if (type=="double") {
//						stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:8 title \"cluster " + String(i) + "\"";
//					}
//					else {
//						stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:11 title \"cluster " + String(i) + "\"";
//					}
//				}
//				stream_gnuplotscript << std::endl;
//
//				// write *_ratios_light_heavy_intensity.eps
//				if (type=="double") {
//					stream_gnuplotscript << "set output \"" + debug_ratios_light_medium_intensity + "\"" << std::endl;
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//					stream_gnuplotscript << "set nokey" << std::endl;
//					stream_gnuplotscript << "set logscale x" << std::endl;
//					stream_gnuplotscript << "set logscale y" << std::endl;
//					stream_gnuplotscript << "set xlabel \'ratio light medium\'" << std::endl;
//					stream_gnuplotscript << "set ylabel \'intensity \'" << std::endl;
//					stream_gnuplotscript << "plot \'" + debug_dat + "\' using 5:6" << std::endl;
//					stream_gnuplotscript << "unset logscale x" << std::endl;
//					stream_gnuplotscript << "unset logscale y" << std::endl;
//					stream_gnuplotscript << std::endl;
//				}
//
//				// write *_ratios_light_medium.eps
//				if (type=="triple") {
//					stream_gnuplotscript << "set output \"" + debug_ratios_light_medium + "\"" << std::endl;
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//					stream_gnuplotscript << "set nokey" << std::endl;
//					stream_gnuplotscript << "set xlabel \'m/Z\'" << std::endl;
//					stream_gnuplotscript << "set ylabel \'ratio light medium\'" << std::endl;
//					stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:5";
//					stream_gnuplotscript << std::endl;
//				}
//
//				// write *_ratios_light_heavy.eps
//				stream_gnuplotscript << "set output \"" + debug_ratios_light_heavy + "\"" << std::endl;
//				if (type=="double") {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				else {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				stream_gnuplotscript << "set nokey" << std::endl;
//				stream_gnuplotscript << "set xlabel \'m/Z\'" << std::endl;
//				stream_gnuplotscript << "set ylabel \'ratio light heavy\'" << std::endl;
//				if (type=="double") {
//					stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:5";
//				}
//				else {
//					stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:6";
//				}
//				stream_gnuplotscript << std::endl;
//
//				// write *_sizes.eps
//				stream_gnuplotscript << "set output \"" + debug_sizes + "\"" << std::endl;
//				if (type=="double") {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				else {
//					stream_gnuplotscript << "set title \"SILACAnalyzer2 " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
//				}
//				stream_gnuplotscript << "set nokey" << std::endl;
//				stream_gnuplotscript << "set xlabel \'cluster ID\'" << std::endl;
//				stream_gnuplotscript << "set ylabel \'cluster size\'" << std::endl;
//				stream_gnuplotscript << "plot \'" + debug_dat + "\' using 1:2";
//				stream_gnuplotscript << std::endl;
//			}
//			stream_gnuplotscript.close();
//		}


		return EXECUTION_OK;
	}
};



//@endcond

int main(int argc, const char** argv )
{
	TOPPSILACAnalyzer2 tool;
	return tool.main(argc,argv);
}
