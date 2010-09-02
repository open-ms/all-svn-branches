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

//filtering
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>

//clustering
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/COMPARISON/CLUSTERING/HashClustering.h>
#include <OpenMS/COMPARISON/CLUSTERING/CentroidLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/QTClustering.h>

//Contrib includes
//#include <gsl/gsl_histogram.h>
//#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_fit.h>
//#include <gsl/gsl_interp.h>
//#include <gsl/gsl_spline.h>
//#include <gsl/gsl_errno.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>


//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>

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
  - type - double (for SILAC pairs), triple (for SILAC triplets), double_triple (for SILAC pairs and SILAC triplets)

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
  - [*_silhouettes.dat] - contains the silhouette values for every generated subtree
  - [*_silhouettes.R] - R script, which creates a silhouette plot of every subtree. It can be executed using <i>Rscript *_silhouettes.R </i>. If R is not installed see <a href="http://www.r-poject.org">http://www.r-poject.org</a>.
  - [*_subtrees.featureXML] - contains the complete set of data points (retention time, m/z, intensity) of all peptide pairs. Data poins contained in the same subtree have the same color.

  The following parameters are straightforward:
  - mass_separation_light_medium - mass gap between light and medium isotopic envelopes [Da] (only relevant for the search for SILAC triplets, i.e. type triple)
  - mass_separation_light_heavy - mass gap between light and heavy isotopic envelopes [Da]
  - charge_min/max - range of charge states

  The remaining parameters should be tuned in the following order:
  - mz_threshold - maximal m/z distance value of two data points belonging to one SILAC feature.
  - rt_threshold - maximal RT distance value of two data points belonging to one SILAC feature.
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
	return a.cluster_id > b.cluster_id;
}

bool clusterCompare( Cluster v1, Cluster v2) {
    return v1.size() > v2.size();
  }


class TOPPSILACAnalyzer2
: public TOPPBase
{
private:

	String type;
	Int charge_min;
	Int charge_max;
	DoubleReal mz_stepwidth;
	DoubleReal rt_stepwidth;
	DoubleReal intensity_cutoff;
	DoubleReal mz_threshold;
	DoubleReal rt_threshold;
	DoubleReal rt_scaling;
	DoubleReal model_deviation;
	DoubleReal label_selection;
	std::set<String> contained_labels;
	Int charge_selection;
	std::vector<std::set<DoubleReal> > filter_values;
	std::map<String,DoubleReal> label_identifiers;
	String in;
	String out;
	String out_visual;
	ConsensusMap all_pairs;
	ConsensusMap all_pairs1;
	FeatureMap<> all_cluster_points;
	FeatureMap<> all_cluster_points1;
	FeatureMap<> subtree_points;


	void handleParameters()
	{
		type = getStringOption_("type");
		String charge_string=getParam_().getValue("algorithm:charge");
		DoubleReal charge_min_temp,charge_max_temp;
		parseRange_(charge_string,charge_min_temp,charge_max_temp);
		charge_min=(Int)charge_min_temp;
		charge_max=(Int)charge_max_temp;
		intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
		mz_threshold = getParam_().getValue("algorithm:mz_threshold");
		rt_threshold = getParam_().getValue("algorithm:rt_threshold");
		rt_scaling = getParam_().getValue("algorithm:rt_scaling");
		model_deviation = getParam_().getValue("algorithm:maximum_model_deviation");
		label_identifiers.insert(std::make_pair("arg6",getParam_().getValue("labels:arg6")));
		label_identifiers.insert(std::make_pair("arg10",getParam_().getValue("labels:arg10")));
		label_identifiers.insert(std::make_pair("lys4",getParam_().getValue("labels:lys4")));
		label_identifiers.insert(std::make_pair("lys6",getParam_().getValue("labels:lys6")));
		label_identifiers.insert(std::make_pair("lys8",getParam_().getValue("labels:lys8")));
		String label_list_string=getParam_().getValue("algorithm:isotopic_labels");
		std::vector<String> label_list;
		boost::split( label_list, label_list_string , boost::is_any_of("[]") );
		 for (std::vector<String>::iterator it=label_list.begin();it!=label_list.end();++it)
		 {
			std::set<DoubleReal> act_filter_values;
			std::vector<String> act_labels;
			boost::split( act_labels, *it , boost::is_any_of(",") );
			for (std::vector<String>::iterator label_it=act_labels.begin();label_it!=act_labels.end() && *label_it!="";++label_it)
			{
				String act_label=*label_it;
				std::remove(act_label.begin(), act_label.end(), ' ');
				std::transform(act_label.begin(),act_label.end(),act_label.begin(),tolower);
				std::map<String,DoubleReal>::iterator pos = label_identifiers.find(act_label);
				if (pos==label_identifiers.end())
					throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,act_label);
				act_filter_values.insert(pos->second);
				contained_labels.insert(act_label);
			}
			if (act_filter_values.size()!=0)
				filter_values.push_back(act_filter_values);
		 }

		in = getStringOption_("in");
		out = getStringOption_("out");
		out_visual = getStringOption_("out_visual");
		String label_selection_string = getParam_().getValue("selected_search:label");
		std::map<String,DoubleReal>::iterator pos = label_identifiers.find(label_selection_string);
		if (label_selection_string=="")
			label_selection=-1;
		else if (pos==label_identifiers.end() || contained_labels.find(label_selection)==contained_labels.end())
			throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,label_selection_string);
		label_selection=pos->second;
		charge_selection = getParam_().getValue("selected_search:charge");
		if ((charge_selection < charge_min || charge_selection > charge_max) && charge_selection>0)
			throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Charge "+charge_selection);
		//output variables
		all_pairs.getFileDescriptions()[0].filename = in;
		all_pairs.getFileDescriptions()[0].label = "light";
		all_pairs.getFileDescriptions()[1].filename = in;
		all_pairs.getFileDescriptions()[1].label = "medium";
		all_pairs.getFileDescriptions()[2].filename = in;
		all_pairs.getFileDescriptions()[2].label = "heavy";
		all_pairs.setExperimentType("silac");
		all_pairs1.getFileDescriptions()[0].filename = in;
		all_pairs1.getFileDescriptions()[0].label = "light";
		all_pairs1.getFileDescriptions()[1].filename = in;
		all_pairs1.getFileDescriptions()[1].label = "medium";
		all_pairs1.getFileDescriptions()[2].filename = in;
		all_pairs1.getFileDescriptions()[2].label = "heavy";
		all_pairs1.setExperimentType("silac");

	}



	std::vector<std::vector<DataPoint> > buildDataStructure(MSExperiment<Peak1D>& exp)
	{

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
		}
		for (i=1; i<exp.size(); ++i)
		{
			rt_spacing.push_back(exp[i].getRT()-exp[i-1].getRT());
		}
		std::sort(mz_spacing.begin(),mz_spacing.end());
		mz_stepwidth=mz_spacing[mz_spacing.size()/2];
//		std::cout << mz_stepwidth << std::endl;
//		mz_stepwidth=0.0005;
		std::sort(rt_spacing.begin(),rt_spacing.end());
		rt_stepwidth=rt_spacing[rt_spacing.size()/2];
		std::list<SILACFilter> filters;
		for (Int charge=charge_min; charge<=charge_max; ++charge)
		{
			std::cout << charge << "+" << std::endl;;
			for (std::vector<std::set<DoubleReal> >::iterator value_vector_it=filter_values.begin();value_vector_it!=filter_values.end();++value_vector_it)
			{
				std::set<DoubleReal> value_set(value_vector_it->begin(),value_vector_it->end());
				for (std::set<DoubleReal>::iterator bla=value_set.begin();bla!=value_set.end();++bla)
				{
					std::cout << *bla << " ";
				}
				std::cout << std::endl;
				filters.push_back(SILACFilter(value_set,charge,model_deviation));
			}

		}

		SILACFiltering filtering(exp,mz_stepwidth,intensity_cutoff);
		filtering.setLogType(log_type_);
		for (std::list<SILACFilter>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
		{
			filtering.addFilter(*filter_it);
		}

		filtering.filterDataPoints(out.substr(0,in.find_first_of('.')));
		std::vector<std::vector<DataPoint> > data;
		for (std::list<SILACFilter>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
		{
			std::set<DoubleReal> envelope_distances=filter_it->getEnvelopeDistances();
			if (charge_selection <=0 || label_selection<=0 || (std::abs(*envelope_distances.rbegin()-label_selection)<std::numeric_limits<DoubleReal>::epsilon() && charge_selection==filter_it->getCharge()))
				data.push_back(filter_it->getElements());
		}
		exp.clear(true);
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

		registerSubsection_("labels","Label configuration");
		registerSubsection_("algorithm","Algorithm parameters section");
		registerSubsection_("selected_search","Label and charge selection for specific search");
	}

	Param getSubsectionDefaults_(const String& section) const
	{
		Param tmp;
		if (section == "labels")
		{

//			tmp.setValue("cleavage",StringList::create("Trypsin,Arg-C,Asp-N,Asp-N_ambic,Chymotrypsin,CNBr,CNBr+Trypsin,Formic_acid,Lys-C,Lys-C/P,PepsinA,Tryp-CNBr,TrypChymo,Trypsin/P,V8-DE,V8-E,semiTrypsin,LysC+AspN,None"),"CL");
			tmp.setValue("arg6", 6.0202, "Arg6 mass shift");
			tmp.setMinFloat("arg6", 0.0);
			tmp.setValue("arg10", 9.9304356, "Arg10 mass shift");
			tmp.setMinFloat("arg10", 0.0);
			tmp.setValue("lys4", 4.0202, "Lys4 mass shift");
			tmp.setMinFloat("lys4", 0.0);
			tmp.setValue("lys6", 6.0202, "Lys6 mass shift");
			tmp.setMinFloat("lys6", 0.0);
			tmp.setValue("lys8", 7.9427178, "Lys8 mass shift");
			tmp.setMinFloat("lys8", 0.0);
		}
		if (section == "algorithm")
		{
			tmp.setValue("isotopic_labels", "", "Isotopic labels. Add an entry for each SILAC type by using comma-seperated label identifiers (see \"labels\" section in advanced parameters)");

			tmp.setValue("charge",":","charge range");

			tmp.setValue("intensity_cutoff", 5000.0, "intensity cutoff");
			tmp.setMinFloat("intensity_cutoff", 0.0);

			tmp.setValue("mz_threshold", 50.0, "mz_threshold");
			tmp.setMinFloat("mz_threshold", 0.0);

			tmp.setValue("rt_threshold", 50.0, "rt_threshold");
			tmp.setMinFloat("rt_threshold", 0.0);

			tmp.setValue("rt_scaling",0.05,"scaling factor of retention times (Cluster height [s] an\ncluster width [Th] should be of the same order. The clustering algorithms work better for\nsymmetric clusters.)");
			tmp.setMinFloat("rt_scaling", 0.0);

			tmp.setValue("maximum_model_deviation",0.1,"Maximal value of which a predicted SILAC feature may deviate from the averagine model.");
			tmp.setMinFloat("maximum_model_deviation",0.0);

		}
		if (section == "selected_search")
		{
			tmp.setValue("label", "", "Label selection");
			tmp.setValue("charge", 0, "Charge selection. 0 for no charge selection. Charge must be in the selected charge range at the \"algorithm\" section.");
		}
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
			debug_trunk = in.substr(0,in.find_first_of('.'));
		}

		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------
		MzMLFile file;
		MSExperiment<Peak1D> exp;

		file.setLogType(log_type_);
		//			file.getOptions().setIntensityRange(DRange<1>(intensity_cutoff,std::numeric_limits<DoubleReal>::max()));
		file.load(in,exp);


		//set input map size (only once)

		exp.updateRanges();


		all_pairs.getFileDescriptions()[0].size = exp.getSize();
		all_pairs.getFileDescriptions()[1].size = exp.getSize();
		all_pairs.getFileDescriptions()[2].size = exp.getSize();

		all_pairs1.getFileDescriptions()[0].size = exp.getSize();
		all_pairs1.getFileDescriptions()[1].size = exp.getSize();
		all_pairs1.getFileDescriptions()[2].size = exp.getSize();


		//-------------------------------------------------------------
		// build SILACData structure
		//-------------------------------------------------------------

		std::vector<std::vector<DataPoint> > data=buildDataStructure(exp);

		//-------------------------------------------------------------
		// Perform clustering
		//-------------------------------------------------------------

		std::vector<std::vector<Real> > silhouettes;
		std::vector<Cluster> clusters;
		std::vector<Tree> subtrees;

		for (std::vector<std::vector<DataPoint> >::iterator data_it=data.begin();data_it!=data.end();++data_it)
		{
			CentroidLinkage method(rt_scaling);
			HashClustering c(*data_it,rt_threshold,mz_threshold,method);
			c.setLogType(log_type_);
			c.performClustering();
			std::vector<Tree> act_subtrees;
			c.getSubtrees(act_subtrees);
			subtrees.insert(subtrees.end(),act_subtrees.begin(),act_subtrees.end());
//			DoubleReal isotope_distance=1.000495/(DoubleReal)data_it->front().charge;
//			QTClustering c(*data_it,rt_threshold, mz_threshold,isotope_distance);
//			c.setLogType(log_type_);
//			std::vector<Cluster> act_clusters=c.performClustering();
			std::vector<Cluster> act_clusters;
			c.createClusters(act_clusters);
			clusters.insert(clusters.end(),act_clusters.begin(),act_clusters.end());
			const std::vector<std::vector<Real> >& act_silhouettes=c.getSilhouetteValues();
			silhouettes.insert(silhouettes.end(),act_silhouettes.begin(),act_silhouettes.end());


		}

		std::sort(clusters.begin(),clusters.end(),clusterCompare);


		if (getFlag_("silac_debug"))
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
			Size subtree_number=1;
			for (std::vector<std::vector<BinaryTreeNode> >::iterator subtree_it=subtrees.begin();subtree_it!=subtrees.end();++subtree_it)
			{
				std::set<DataPoint*> leafs;
				for (std::vector<BinaryTreeNode>::iterator tree_it=subtree_it->begin(); tree_it!= subtree_it->end(); ++tree_it)
				{
					leafs.insert(tree_it->data1);
					leafs.insert(tree_it->data2);
				}
				for (std::set<DataPoint*>::iterator leafs_it=leafs.begin();leafs_it!=leafs.end();++leafs_it)
				{
					//visualize the light variant
					Feature tree_point;
					tree_point.setRT((*leafs_it)->rt);
					tree_point.setMZ((*leafs_it)->mz);
					tree_point.setIntensity((*leafs_it)->intensities[0][0]);
					tree_point.setCharge((*leafs_it)->charge);
					tree_point.setMetaValue("subtree",subtree_number);
					tree_point.setMetaValue("color",colors[subtree_number%colors.size()]);
					subtree_points.push_back(tree_point);
				}
				++subtree_number;
			}

			// required, as somehow the order of features on some datasets between Win & Linux is different and thus the TOPPtest might fail
			subtree_points.sortByPosition();
		}

		std::map<Size,String> silac_types;
		silac_types.insert(std::make_pair(0,"Singlet"));
		silac_types.insert(std::make_pair(1,"Doublet"));
		silac_types.insert(std::make_pair(2,"Triplet"));
		silac_types.insert(std::make_pair(3,"Quadruplet"));

		if (out!="")
		{
			Size id=0;
			for(std::vector<Cluster>::iterator cluster_it=clusters.begin();cluster_it!=clusters.end();++cluster_it)
			{
				DoubleReal rt = 0.0;
				DoubleReal mz = 0.0;
				DoubleReal total_intensity = 0.0;
				Size mass_shifts_size=(*(cluster_it->begin()))->mass_shifts.size();
				std::vector<DoubleReal> max_intensities(mass_shifts_size,0.0);
				String silac_type=silac_types[mass_shifts_size];
				Int charge=(*(cluster_it->begin()))->charge;
				// intensity vectors used for linear regression
				std::vector<std::vector<DoubleReal> > intensities(mass_shifts_size);
				DoubleReal max_intensity=0.0;
				for (Cluster::iterator el_it=cluster_it->begin();el_it!=cluster_it->end();++el_it)
				{

					std::vector<std::vector<DoubleReal> >& element_intensities=(*el_it)->intensities;
					total_intensity+=element_intensities[0][0];
					rt += element_intensities[0][0]*(*el_it)->rt;

					for (Size i=0;i<mass_shifts_size;++i)
					{
						std::vector<DoubleReal>& act_intensities=element_intensities[i];
						intensities[i].insert(intensities[i].end(),act_intensities.begin(),act_intensities.end());

						std::vector<DoubleReal>::iterator max_position=std::max_element(act_intensities.begin(),act_intensities.end());
						if (*max_position>max_intensities[i])
						{
							max_intensities[i]=*max_position;
							mz = (*el_it)->mz;
						}
					}

				}
				rt /= total_intensity; // average retention time

				ConsensusFeature consensus_feature;
				consensus_feature.setRT(rt);
				consensus_feature.setMZ(mz);
				consensus_feature.setIntensity(max_intensities[0]);
				consensus_feature.setCharge(charge);
				consensus_feature.setQuality(0);

				for (Size j=0;j<mass_shifts_size;++j)
				{
					FeatureHandle handle;
					handle.setRT(rt);
					handle.setMZ(mz+(*(cluster_it->begin()))->mass_shifts[j]);
					handle.setIntensity(max_intensities[j]);
					handle.setCharge(charge);
					handle.setMapIndex(j);
					handle.setUniqueId(id);
					consensus_feature.insert(handle);
				}
				all_pairs.push_back(consensus_feature);
				++id;
			}
		}

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
			for (std::vector<std::vector<DataPoint> >::iterator data_it=data.begin(); data_it!= data.end(); ++data_it)
			{
				for (std::vector<DataPoint>::iterator it=data_it->begin(); it!= data_it->end(); ++it)
				{
					//visualize the light variant
					Feature cluster_point;
					cluster_point.setRT(it->rt);
					cluster_point.setMZ(it->mz);
					cluster_point.setIntensity(it->intensities[0][0]);
					cluster_point.setCharge(it->charge);
					cluster_point.setOverallQuality(it->quality);
					cluster_point.setQuality(0,it->quality);
					cluster_point.setMetaValue("SILAC type",silac_types[it->mass_shifts.size()-1]);
					String mass_shift="";
					for (std::vector<DoubleReal>::iterator shift_it=it->mass_shifts.begin()+1;shift_it!=it->mass_shifts.end();++shift_it)
					{
						mass_shift+=((String)*shift_it)+" ";
					}
					if (mass_shift!="")
						cluster_point.setMetaValue("Mass shift (l/h)",mass_shift);
					cluster_point.setMetaValue("Cluster id",it->cluster_id);
					cluster_point.setMetaValue("color",colors[it->cluster_id%colors.size()]);
					if (getFlag_("silac_debug"))
					{
						cluster_point.setMetaValue("Cluster size",it->cluster_size);
						cluster_point.setMetaValue("feature_id",it->feature_id);
					}

					all_cluster_points.push_back(cluster_point);
				}
			}
			// required, as somehow the order of features on some datasets between Win & Linux is different and thus the TOPPtest might fail
			all_cluster_points.sortByPosition();
		}




		//-------------------------------------------------------------
		// generate debug output
		//-------------------------------------------------------------
		// strings repeatedly used in debug output


		if (getFlag_("silac_debug"))
		{
//			// names of dat files
			String debug_clusters_dat = debug_trunk + "_cluster_sizes.dat";
//
//			// write all cluster data points to *_clusters.dat
//
			std::ofstream stream_clusters(debug_clusters_dat.c_str());
			for(std::vector<Cluster>::iterator cluster_it=clusters.begin();cluster_it!=clusters.end();++cluster_it)
			{
				stream_clusters << cluster_it->size() << std::endl;
			}
//			stream_clusters << "cluster_id rt mz int" << std::endl;
//			Int current_id = -1;
//			for (std::vector<std::vector<DataPoint> >::iterator data_it=data.begin(); data_it!= data.end(); ++data_it)
//			{
//				for (std::vector<DataPoint>::iterator it=data_it->begin(); it!= data_it->end(); ++it)
//				{
//					if (it->cluster_id != current_id) stream_clusters << std::endl << std::endl;
//					stream_clusters << it->cluster_id << " " << it->rt << " " << it->mz << " " << it->intensities[0] << std::endl;
//					current_id = it->cluster_id;
//				}
//			}

			stream_clusters.close();

			String debug_silhouettes_dat = debug_trunk + "_silhouettes.dat";
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

//			String debug_intensities_dat = debug_trunk + "_intensities.dat";
//			std::ofstream stream_intensities(debug_intensities_dat.c_str());
//			for (std::vector<std::vector<DataPoint> >::iterator data_it=data.begin(); data_it!= data.end(); ++data_it)
//			{
//				for (std::vector<DataPoint>::iterator it=data_it->begin(); it!= data_it->end(); ++it)
//				{
//					std::vector<DoubleReal> intensities=it->intensities;
//					stream_intensities << intensities[0] << "\t" << intensities[3] << std::endl;
//				}
//			}
//			for(std::vector<Cluster>::iterator cluster_it=clusters.begin();cluster_it!=clusters.end();++cluster_it)
//			{
//				for (Cluster::iterator el_it=cluster_it->begin();el_it!=cluster_it->end();++el_it)
//				{
//					std::vector<DoubleReal> intensities=(*el_it)->intensities;
//					stream_intensities << intensities[0] << "\t" << intensities[3] << std::endl;
//				}
//				stream_intensities << "*" << std::endl;
//			}
//			stream_intensities.close();
//
//			Size data_size=0;
//			for (std::vector<std::vector<DataPoint> >::iterator data_it=data.begin(); data_it!= data.end(); ++data_it)
//			{
//				data_size+=data_it->size();
//			}

			String debug_silhouettes_r = debug_trunk + "_silhouettes.R";
			std::ofstream stream_r(debug_silhouettes_r.c_str());
			stream_r << "# SILACAnalyzer debug script to run in R, which creates a pdf file containing a silhouette plot for every clustered subtree"<<std::endl;
			stream_r << "# This script can be executed from the command line using 'Rscript "<<debug_silhouettes_r<<"'"<<std::endl;
			stream_r << "con <- file(\"" << debug_silhouettes_dat << "\",\"r\")" << std::endl;
			stream_r << "lines <- readLines(con, n=-1)" << std::endl;
			String pdf_name= debug_trunk + "_silhouettes.pdf";
			String pdf_title= debug_trunk;
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
//			stream_r << "#The following part creates an intensity vs intensity plot" << std::endl;
//			stream_r << "data<-read.table(file=\"" << debug_intensities_dat << "\")" << std::endl;
//			pdf_name= debug_trunk + "_intensities.pdf";
//			stream_r << "pdf(file=\""<< pdf_name <<"\")" << std::endl;
//			stream_r << "plot(data,pch=20,xlab=\"Heavy intensity\",ylab=\"Light Intensity\",main=\"Intensity vs. Intensity plot ("<< data_size <<" features)\")" << std::endl;
//			stream_r << "dev.off()" << std::endl;
			stream_r.close();

			//				// write ratios of all cluster to *.dat
//			String debug_dat = debug_trunk + ".dat";
//			std::ofstream stream_ratios(debug_dat.c_str());
//			if (type=="double") {
//				stream_ratios << "cluster_id cluster_size rt mz ratio_light_heavy intensity" << std::endl;
//			}
//			else {
//				stream_ratios << "cluster_id cluster_size rt mz ratio_light_medium ratio_light_heavy intensity" << std::endl;
//			}
//			Size i=0;
//			for(std::vector<Cluster>::iterator cluster_it=clusters.begin();cluster_it!=clusters.end();++cluster_it)
//			{
//				DoubleReal rt = 0.0;
//				DoubleReal mz = 0.0;
//				DoubleReal int_l = 0.0;
//				DoubleReal int_m = 0.0;
//				DoubleReal int_h = 0.0;
//				Int silac_type=(*(cluster_it->begin()))->silac_type;
//				// intensity vectors used for linear regression
//				std::vector<DoubleReal> i1(3*cluster_it->size());
//				std::vector<DoubleReal> i2(3*cluster_it->size());
//				std::vector<DoubleReal> i3(3*cluster_it->size());
//				UInt j=0;
//				for (Cluster::iterator el_it=cluster_it->begin();el_it!=cluster_it->end();++el_it)
//				{
//					std::vector<DoubleReal> intensities=(*el_it)->intensities;
//					rt += (*el_it)->rt;
//					i1[3*j] = intensities[0];
//					i1[3*j+1] = intensities[1];
//					i1[3*j+2] = intensities[2];
//					i2[3*j] = intensities[3];
//					i2[3*j+1] = intensities[4];
//					i2[3*j+2] = intensities[5];
//					if ((*el_it)->silac_type == DataPoint::TRIPLE)
//					{
//						i3[3*j] = intensities[6];
//						i3[3*j+1] = intensities[7];
//						i3[3*j+2] = intensities[8];
//					}
//
//					std::pair<DoubleReal,DoubleReal> low_maximum=clusterMaximum(**el_it,0,1,2);
//					if (low_maximum.second > int_l)
//					{
//						int_l=low_maximum.second;
//						mz=low_maximum.first;
//					}
//
//					if ((*el_it)->silac_type == DataPoint::DOUBLE)
//					{
//						std::pair<DoubleReal,DoubleReal> heavy_maximum=clusterMaximum(**el_it,3,4,5);
//						if (low_maximum.second > int_l)
//						{
//							int_h=low_maximum.second;
//							mz=low_maximum.first;
//						}
//					}
//					else
//					{
//						std::pair<DoubleReal,DoubleReal> medium_maximum=clusterMaximum(**el_it,3,4,5);
//						if (medium_maximum.second > int_l)
//						{
//							int_m=medium_maximum.second;
//							mz=medium_maximum.first;
//						}
//						std::pair<DoubleReal,DoubleReal> heavy_maximum=clusterMaximum(**el_it,6,7,8);
//						if (heavy_maximum.second > int_l)
//						{
//							int_h=heavy_maximum.second;
//							mz=heavy_maximum.first;
//						}
//					}
//					++j;
//				}
//				rt /= (DoubleReal)(cluster_it->size()); // average retention time
//				if (silac_type==DataPoint::TRIPLE)
//				{
//					Math::LinearRegression linear_reg_light_heavy;
//					linear_reg_light_heavy.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i3.begin());
//					Math::LinearRegression linear_reg_light_medium;
//					linear_reg_light_medium.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
//					stream_ratios << i << " " << cluster_it->size() << " " << rt << " " << mz << " " << linear_reg_light_medium.getSlope() << " " << linear_reg_light_heavy.getSlope() << " " << *max_element(i1.begin(),i1.end()) + *max_element(i2.begin(),i2.end()) << std::endl;
//				}
//				else
//				{
//					Math::LinearRegression linear_reg_light_heavy;
//					linear_reg_light_heavy.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
//					stream_ratios << i << " " << cluster_it->size() << " " << rt << " " << mz << " " << linear_reg_light_heavy.getSlope() << " " << *max_element(i1.begin(),i1.end()) + *max_element(i2.begin(),i2.end()) << std::endl;
//				}
//				++i;
//			}
		}


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

//		String qt_out;
//		if (out.has('.'))
//		{
//			qt_out = out.substr(0,out.find_first_of('.'));
//		}
//
//		if (out!="")
//		{
//
//			// assign unique ids
//			all_pairs1.applyMemberFunction(&UniqueIdInterface::setUniqueId);
//
//			//annotate output with data processing info
//			addDataProcessing_(all_pairs1, getProcessingInfo_(DataProcessing::QUANTITATION));
//
//			ConsensusXMLFile c_file;
//			c_file.store(qt_out+"_qt.consensusXML",all_pairs1);
//		}

		if (out_visual!="")
		{
			// assign unique ids
			all_cluster_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			FeatureXMLFile f_file;
			f_file.store(out_visual,all_cluster_points);
		}
//		if (out_visual.has('.'))
//		{
//			qt_out = out.substr(0,out_visual.find_first_of('.'));
//		}
//		if (out_visual!="")
//		{
//			// assign unique ids
//			all_cluster_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);
//
//			FeatureXMLFile f_file;
//			f_file.store(qt_out+"_qt.featureXML",all_cluster_points1);
//		}
		if (getFlag_("silac_debug"))
		{
			// assign unique ids
			subtree_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			FeatureXMLFile t_file;
			t_file.store(debug_trunk+ "_subtrees.featureXML",subtree_points);
		}
		return EXECUTION_OK;
	}
};



//@endcond

int main(int argc, const char** argv )
{
	TOPPSILACAnalyzer2 tool;
	return tool.main(argc,argv);
}
