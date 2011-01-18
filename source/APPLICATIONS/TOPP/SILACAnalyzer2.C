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
// $Authors: Lars Nilse, Steffen Sass, Holger Plattfaut $
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
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <OpenMS/DATASTRUCTURES/DataPoint.h>

#ifdef _OPENMP
#ifdef OPENMS_WINDOWSPLATFORM
#include <omp.h>
#endif
#endif

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <fstream>
#include <iomanip>

using namespace OpenMS;
using namespace std;

typedef vector<BinaryTreeNode> Tree;
typedef vector<DataPoint*> Cluster;

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

bool clusterCompare( Cluster v1, Cluster v2)
{
   return v1.size() > v2.size();
}

class TOPPSILACAnalyzer2
: public TOPPBase
{

	private:

		// input and output files
		String in;
		String out;
		String out_clusters;
    String out_filters_mzML;
    String out_filters;
    String in_filters;

		// section "sample"
		String selected_labels;
		Int charge_min;
		Int charge_max;
		Int missed_cleavages;
		Int isotopes_per_peptide_min;
		Int isotopes_per_peptide_max;
		DoubleReal mz_stepwidth;

		// section "algorithm"
		DoubleReal mz_threshold;
		DoubleReal intensity_cutoff;
		DoubleReal rt_threshold;
		DoubleReal rt_scaling;
		DoubleReal model_deviation;

		// section "labels"
    map<String, DoubleReal> label_identifiers;
    vector<vector <String> > SILAClabels; // list of SILAC labels, e.g. selected_labels="[Lys4,Arg6][Lys8,Arg10]" => SILAClabels[0][1]="Arg6"
    vector<vector <DoubleReal> > massShifts; // list of mass shifts

/*		// section "out_clusters"
		Int mass_shift_out_clusters;
		DoubleReal label_out_clusters_mass_shift;
		DoubleReal ms_final_2;	
		Int charge_out_clusters;
		bool out_clusters_flag;
		bool out_clusters_flag_2;
*/
    vector<vector<DataPoint> > data;
		ConsensusMap all_pairs;
		FeatureMap<> all_cluster_points;
		FeatureMap<> subtree_points;
    MSExperiment<Peak1D> filter_exp;


	public:
		TOPPSILACAnalyzer2()
		: TOPPBase("SILACAnalyzer2","Determination of peak ratios in LC-MS data",true)
		{
		}



  //--------------------------------------------------
	// set structure of ini file
  //--------------------------------------------------

	void registerOptionsAndFlags_()
	{
    // create flag for input file (.mzML)
    registerInputFile_("in", "<file>", "", "input file");
		setValidFormats_("in", StringList::create("mzML"));
		// create flag for output file (.consensusXML)
		registerOutputFile_("out", "<file>", "", "output file", false);
		setValidFormats_("out", StringList::create("consensusXML"));
		// create optional flag for additional clusters output file (.featureXML)
		registerOutputFile_("out_clusters", "<file>", "", "Additional output file containing all clusters differed by colours.", false, true);
		setValidFormats_("out_clusters", StringList::create("featureXML"));
    // create optional flag for additional filters output file (.mzML)
    registerOutputFile_("out_filters_mzML", "<file>", "", "Additional output file containing all points that passed the filters as mzML.", false, true);
    setValidFormats_("out_filters_mzML", StringList::create("mzML"));
    // create optional flag for additional output file (.txt) to store filter results
    registerOutputFile_("out_filters", "<file>", "", "Additional output file containing all points that passed the filters as txt. Suitable as input for \"in_filters\" to perform clustering without preceding filtering process.", false, true);
    //setValidFormats_("out_filters", StringList::create("txt"));
    // create optional flag for additional input file (.txt) to load filter results
    registerOutputFile_("in_filters", "<file>", "", "Additional input file containing all points that passed the filters as txt. Use output from \"out_filters\" to perform clustering only.", false, true);
    //setValidFormats_("in_filters", StringList::create("txt"));

		// create flag for additional debug outputs
		registerFlag_("silac_debug","Enables writing of debug information", true);

		// create section "labels" for adjusting masses of labels
		registerSubsection_("labels", "Isotopic labels that can be selected for section \"sample\".");
		// create section "sample" for adjusting sample parameters
		registerSubsection_("sample", "Parameter adjusting for your sample.");
		// create section "algorithm" for adjusting algorithm parameters
		registerSubsection_("algorithm", "Algorithm parameters section.");
		// create section "out_clusters" for adjusting parameters for additional output file
//		registerSubsection_("out_clusters", "Parameters for filtering out one specific label and one specific charge state.");
	}


	// create prameters for sections (set default values and restrictions)
	Param getSubsectionDefaults_(const String& section) const
	{
		Param defaults;



		//--------------------------------------------------
		// section labels
	  //--------------------------------------------------

		if (section == "labels")
		{
			// create labels that can be chosen in section "sample/labels"
			defaults.setValue("Arg6", 6.0202, "Arg6 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Arg6", 0.0);
			defaults.setValue("Arg10", 10.0202, "Arg10 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Arg10", 0.0);
			defaults.setValue("Lys4", 4.0202, "Lys4 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Lys4", 0.0);
			defaults.setValue("Lys6", 6.0202, "Lys6 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Lys6", 0.0);
			defaults.setValue("Lys8", 8.0202, "Lys8 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Lys8", 0.0);
			defaults.setValue("Methyl4", 4.0202, "Methyl4 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Methyl4", 0.0);
			defaults.setValue("Methyl8", 8.0202, "Methyl8 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Methyl8", 0.0);
			defaults.setValue("Methyl12", 12.0202, "Methyl12 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Methyl12", 0.0);
			defaults.setValue("Methyl16", 16.0202, "Methyl16 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Methyl16", 0.0);
			defaults.setValue("Methyl24", 24.0202, "Methyl24 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Methyl24", 0.0);
			defaults.setValue("Methyl32", 32.0202, "Methyl32 mass shift", StringList::create("advanced"));
			defaults.setMinFloat("Methyl32", 0.0);
		}



	  //--------------------------------------------------
		// section sample
	  //--------------------------------------------------

		if (section == "sample")
		{
      defaults.setValue("labels", "[Arg6]", "Specify the labels for your sample. Doublets must be of style [label,label,...]. Triplets must be of style [label,label,...][label,label,..]. Each pair of brackets can contain one or more labels. See section \"labels\" in advanced parameters for allowed labels.");
			defaults.setValue("charge", "2:3", "Specify the charge range for your sample (charge_min:charge_max).");
			defaults.setValue("missed_cleavages", 0 , "Specify the maximum number of missed cleavages.");
			defaults.setValue("isotopes_per_peptide", "3:4", "Specify the range of isotopes per peptide for your sample (isotopes_per_peptide_min:isotopes_per_peptide_max).", StringList::create("advanced"));
//			defaults.setValue("mz_stepwidth", 0.0, "Select 0 for automatic calculation depending on data set (default and recommended). Adjust this parameter by hand only if the deafult calculation leads to bad results.", StringList::create("advanced"));
//			defaults.setMinFloat("mz_stepwidth", 0.0);
		}



	  //--------------------------------------------------
		// section algorithm
	  //--------------------------------------------------

		if (section == "algorithm")
		{
			defaults.setValue("mz_threshold", 0.1, "Specify an upper bound for your m/z range in [Th]. I.e. the range peptides have to elute in to be considered as one peptide.");
			defaults.setMinFloat("mz_threshold", 0.0);
			defaults.setValue("intensity_cutoff", 0.0, "Specify a threshold for intensity. All peaks below that threshold are not considered.");
			defaults.setMinFloat("intensity_cutoff", 0.0);
			defaults.setValue("rt_threshold", 50.0, "Specify an upper bound for your Retention time range in [s]. I.e. the range peptides have to elute in to be considered as one peptide.");
			defaults.setMinFloat("rt_threshold", 0.0);
			defaults.setValue("rt_scaling", 0.002, "Scaling factor for retention times (Cluster height [s] and cluster width [Th] should be of the same order, because the clustering algorithm works better for symmetric clusters.) In the majority of cases simply divide mz_threshold through rt_threshold.");
			defaults.setMinFloat("rt_scaling", 0.0);
			defaults.setValue("maximum_model_deviation", 0.8, "Maximum value of which a predicted SILAC feature may deviate from the averagine model. Value is log calculated. Thus, 0 means total conformity.");
			defaults.setMinFloat("maximum_model_deviation", 0.0);
		}



	  //--------------------------------------------------
		// section out_clusters
	  //--------------------------------------------------
/*
		if (section == "out_clusters")
		{
			defaults.setValue("mass_shift", 0, "Select a mass_shift. Mass shift numbering starts with 1. 0 for no mass shift selection.]", StringList::create("advanced"));
			defaults.setValue("charge", 0, "Select only one charge state from charge range in section \"sample\". 0 for no charge selection.", StringList::create("advanced"));
		}
*/   
	 	return defaults;
	}



  //--------------------------------------------------
	// handle parameters (read in and format given parameters)
  //--------------------------------------------------

	void handleParameters()
	{
    // get input file (.mzML)
		in = getStringOption_("in");
    // get name of output file (.consensusXML)
		out = getStringOption_("out");
    // get name of additional clusters output file (.featureXML)
		out_clusters = getStringOption_("out_clusters");
    // get name of additional filters output file (.mzML)
    out_filters_mzML = getStringOption_("out_filters_mzML");
    // get name of additional filters output file (.txt)
    out_filters = getStringOption_("out_filters");
    // get name of additional filters input file (.txt)
    in_filters = getStringOption_("in_filters");



		//--------------------------------------------------
		// section labels
	  //--------------------------------------------------
				
		// create map of pairs (label as string, mass shift as double)
    label_identifiers.insert(make_pair("Arg6", getParam_().getValue("labels:Arg6")));
    label_identifiers.insert(make_pair("Arg10", getParam_().getValue("labels:Arg10")));
    label_identifiers.insert(make_pair("Lys4", getParam_().getValue("labels:Lys4")));
    label_identifiers.insert(make_pair("Lys6", getParam_().getValue("labels:Lys6")));
    label_identifiers.insert(make_pair("Lys8", getParam_().getValue("labels:Lys8")));
    label_identifiers.insert(make_pair("Methyl4", getParam_().getValue("labels:Methyl4")));
    label_identifiers.insert(make_pair("Methyl8", getParam_().getValue("labels:Methyl8")));
    label_identifiers.insert(make_pair("Methyl12", getParam_().getValue("labels:Methyl12")));
    label_identifiers.insert(make_pair("Methyl16", getParam_().getValue("labels:Methyl16")));
    label_identifiers.insert(make_pair("Methyl24", getParam_().getValue("labels:Methyl24")));
    label_identifiers.insert(make_pair("Methyl32", getParam_().getValue("labels:Methyl32")));

		// create iterators for all labels to get corresponding mass shift
    map<String,DoubleReal>::iterator arg6 = label_identifiers.find("Arg6");
    map<String,DoubleReal>::iterator arg10 = label_identifiers.find("Arg10");
    map<String,DoubleReal>::iterator lys4 = label_identifiers.find("Lys4");
    map<String,DoubleReal>::iterator lys6 = label_identifiers.find("Lys6");
    map<String,DoubleReal>::iterator lys8 = label_identifiers.find("Lys8");
    map<String,DoubleReal>::iterator methyl4 = label_identifiers.find("Methyl4");
    map<String,DoubleReal>::iterator methyl8 = label_identifiers.find("Methyl8");
    map<String,DoubleReal>::iterator methyl12 = label_identifiers.find("Methyl12");
    map<String,DoubleReal>::iterator methyl16 = label_identifiers.find("Methyl16");
    map<String,DoubleReal>::iterator methyl24 = label_identifiers.find("Methyl24");
    map<String,DoubleReal>::iterator methyl32 = label_identifiers.find("Methyl32");

		// create string of all labels from advanced section "labels"
		String labels = "Arg6 Arg10 Lys4 Lys6 Lys8 Methyl4 Methyl8 Methyl12 Methyl16 Methyl24 Methyl32";

				

		//--------------------------------------------------
		// section sample
	  //--------------------------------------------------

		// get selected labels
		selected_labels = getParam_().getValue("sample:labels");
		
		// get selected missed_cleavages
		missed_cleavages = getParam_().getValue("sample:missed_cleavages");

		// get selected charge range
    String charge_string = getParam_().getValue("sample:charge");
		DoubleReal charge_min_temp, charge_max_temp;
		parseRange_(charge_string, charge_min_temp, charge_max_temp);
    charge_min = (Int)charge_min_temp;
    charge_max = (Int)charge_max_temp;

		// get selected isotopes range
    String isotopes_per_peptide_string = getParam_().getValue("sample:isotopes_per_peptide");
		DoubleReal isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
		parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min = (Int)isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max = (Int)isotopes_per_peptide_max_temp;

		// get selected mz_stepwidth
//		mz_stepwidth = getParam_().getValue("sample:mz_stepwidth");



  	//--------------------------------------------------
		// section algorithm
	  //--------------------------------------------------

		// get selected mz_threshold
		mz_threshold = getParam_().getValue("algorithm:mz_threshold");
		// get selected intensity_cutoff
		intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
		// get selected rt_threshold
		rt_threshold = getParam_().getValue("algorithm:rt_threshold");
		// get selected rt_scaling
		rt_scaling = getParam_().getValue("algorithm:rt_scaling");
		// get selected maximum_model_deviation
		model_deviation = getParam_().getValue("algorithm:maximum_model_deviation");



	  //--------------------------------------------------
		// section out clusters
	  //--------------------------------------------------
/*
		// get selected mass_shift
		mass_shift_out_clusters = getParam_().getValue("out_clusters:mass_shift");

		// get selected cahrge state
		charge_out_clusters = getParam_().getValue("out_clusters:charge");
*/
		

	  //--------------------------------------------------
		// calculate all possible mass shifts for labelets from section "sample:labels"
	  //--------------------------------------------------
		
		// split string of SILAC labels (selected_labels) and save in a list (SILAClabels) 
    vector<String> tempList; // temporary list of strings for SILAC labelets, e.g. "Lys6,Arg8"
		boost::split( tempList, selected_labels, boost::is_any_of("[](){}") ); // any bracket allowed to separate labelets
		for (unsigned i = 0; i < tempList.size(); i++)
		{
			if (tempList[i] != "")
			{
        vector<String> tempLabels;
				boost::split( tempLabels, tempList[i], boost::is_any_of(",;: ") ); // various separators allowed to separate labels
				SILAClabels.push_back(tempLabels);
			}
		}

    cout << endl;
		// print SILAC labels
		for (unsigned i = 0; i < SILAClabels.size(); i++)
		{
      cout << "SILAC label " << i+1 << ":   ";
			for (unsigned j = 0; j < SILAClabels[i].size(); j++)
			{
        cout << SILAClabels[i][j] << " ";
			}
      cout << endl;
		}
    cout << endl;

		// check if all selected labels are included in advanced section "labels"
		for (unsigned i = 0; i < SILAClabels.size(); i++)
		{
			for (unsigned j = 0; j < SILAClabels[i].size(); j++)
			{
				int found = labels.find(SILAClabels[i][j]);
				if (found < 0)
        {
					throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,SILAClabels[i][j]);
        }
			}
		}

		// generate list of mass shifts
		for (Int ArgPerPeptide = 0; ArgPerPeptide <= missed_cleavages + 1; ArgPerPeptide++)
		{
			for (Int LysPerPeptide = 0; LysPerPeptide <= missed_cleavages + 1; LysPerPeptide++)
			{
				for (Int MethylPerPeptide = 0; MethylPerPeptide <= missed_cleavages + 1; MethylPerPeptide++)
				{
					if ( ArgPerPeptide + LysPerPeptide + MethylPerPeptide > 0 && ArgPerPeptide + LysPerPeptide + MethylPerPeptide <= missed_cleavages + 1 )
					{
            vector<DoubleReal> massShiftVector;
						for (unsigned i = 0; i < SILAClabels.size(); i++)
						{
							DoubleReal massShift = 0;
							// Considering the case of an amino acid (e.g. LysPerPeptide != 0) for which no label is present (e.g. Lys4There + Lys8There == 0) makes no sense. Therefore each amino acid will have to give its "Go Ahead" before the shift is calculated.
							bool goAhead_Lys = FALSE;
							bool goAhead_Arg = FALSE;
							bool goAhead_Methyl = FALSE;							
							for (unsigned j = 0; j < SILAClabels[i].size(); j++)
							{
								Int Arg6There = 0;	// Is Arg6 in the SILAC label?
								Int Arg10There = 0;
								Int Lys4There = 0;
								Int Lys6There = 0;
								Int Lys8There = 0;
								Int Methyl4There = 0;
								Int Methyl8There = 0;
								Int Methyl12There = 0;
								Int Methyl16There = 0;
								Int Methyl24There = 0;
								Int Methyl32There = 0;

								if ( SILAClabels[i][j].find("Arg6") == 0 ) Arg6There = 1;
								if ( SILAClabels[i][j].find("Arg10") == 0 ) Arg10There = 1;
								if ( SILAClabels[i][j].find("Lys4") == 0 ) Lys4There = 1;
								if ( SILAClabels[i][j].find("Lys6") == 0 ) Lys6There = 1;
								if ( SILAClabels[i][j].find("Lys8") == 0 ) Lys8There = 1;
								if ( SILAClabels[i][j].find("Methyl4") == 0 ) Methyl4There = 1;
								if ( SILAClabels[i][j].find("Methyl8") == 0 ) Methyl8There = 1;
								if ( SILAClabels[i][j].find("Methyl12") == 0 ) Methyl12There = 1;
								if ( SILAClabels[i][j].find("Methyl16") == 0 ) Methyl16There = 1;
								if ( SILAClabels[i][j].find("Methyl24") == 0 ) Methyl24There = 1;
								if ( SILAClabels[i][j].find("Methyl32") == 0 ) Methyl32There = 1;

								goAhead_Arg = goAhead_Arg || !( (ArgPerPeptide != 0 && Arg6There + Arg10There == 0) );
								goAhead_Lys = goAhead_Lys || !( (LysPerPeptide != 0 && Lys4There + Lys6There + Lys8There == 0) ); 
								goAhead_Methyl = goAhead_Methyl || !( (MethylPerPeptide != 0 && Methyl4There + Methyl8There + Methyl12There + Methyl16There + Methyl24There + Methyl32There == 0) );

								massShift = massShift + ArgPerPeptide * ( Arg6There * (arg6->second) + Arg10There * (arg10->second) ) + LysPerPeptide * ( Lys4There * (lys4->second) + Lys6There * (lys6->second) + Lys8There * (lys8->second) ) +  MethylPerPeptide * (Methyl4There * (methyl4->second) + Methyl8There * (methyl8->second) + Methyl12There * (methyl12->second) + Methyl16There * (methyl16->second) + Methyl24There * (methyl24->second) + Methyl32There * (methyl32->second) );
							}

							if (goAhead_Arg && goAhead_Lys && goAhead_Methyl)
							{
								massShiftVector.push_back(massShift);
							}
						}

            if (!massShiftVector.empty())
            {
              massShifts.push_back(massShiftVector);
            }
					}
				}
			}
		}
		
    // sort the mass shift vector
    sort(massShifts.begin(), massShifts.end());
		
		// print mass shifts
		for (unsigned i = 0; i < massShifts.size(); i++)
		{
      cout << "mass shift " << i+1 << ":   ";
			for (unsigned j = 0; j < massShifts[i].size(); j++)
			{
        cout << massShifts[i][j] << " ";
			}
      cout << endl;
		}
    cout << endl;
				


	  //--------------------------------------------------
		// check input from section "out_clusters"
	  //--------------------------------------------------
/*
		out_clusters_flag = false;
		out_clusters_flag_2 = false;

		// check if one charge state is selected
		if (mass_shift_out_clusters > 0 && charge_out_clusters == 0)
		{
			String out_clusters_exception;				
			out_clusters_exception = "Select charge for out_clusters.";
			throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,out_clusters_exception);
		}

		// check if one mass shift is selected
		if (mass_shift_out_clusters == 0 && charge_out_clusters != 0)
		{
			String out_clusters_exception;				
			out_clusters_exception = "Select mass_shift for out_clusters.";
			throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,out_clusters_exception);
		}

		if (mass_shift_out_clusters > 0 && charge_out_clusters != 0)
		{
			// ckeck if selected mass shift is valid (i.e. if selected mass shift is inside possible mass shifts depending on sections "sample:labels" and "sample:missed_cleavages )
			if (mass_shift_out_clusters > (Int)massShifts.size())
			{
				String out_clusters_exception;
				out_clusters_exception = "Selected mass shift is inavlid.";
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,out_clusters_exception);
			} else

			// check if selected charge state is inside charge range from section "sample:charge"
			if (charge_out_clusters < charge_min || charge_out_clusters > charge_max)
			{
				String out_clusters_exception;				
				out_clusters_exception = "Selected charge state \"" + (String)charge_out_clusters + "\" for out_clusters is not in charge range.";
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,out_clusters_exception);
			}
			else
			{
				out_clusters_flag = true;
			}
		}

		if (out_clusters_flag == true)
		{
			// create string containing selected mass shift (necessary to write only corresponding clusters to featureXML output)
			String ms_5;
			String ms_6;
			for(unsigned i = 0; i < massShifts[mass_shift_out_clusters-1].size(); i++)
			{
				ms_5 += "0" + (String)massShifts[mass_shift_out_clusters-1][i];
			}
				for (int i = 0; i < 2; ++i)
				{
					int found_ms = ms_5.find(".");
					if (found_ms > 0)
					{
						ms_6 += ms_5.substr(found_ms - 2, 3);
						ms_5.erase(found_ms, 1);
			
						int found_ms_2 = ms_6.find(".");
						ms_6.erase(found_ms_2, 1);
					}
				}

				ms_6.insert(2, ".");
        istringstream stm;
				stm.str(ms_6);
				stm >> ms_final_2;
		}
*/
		// get output variables
		all_pairs.getFileDescriptions()[0].filename = in;
		all_pairs.getFileDescriptions()[0].label = "light";
		all_pairs.getFileDescriptions()[1].filename = in;
		all_pairs.getFileDescriptions()[1].label = "medium";
		all_pairs.getFileDescriptions()[2].filename = in;
		all_pairs.getFileDescriptions()[2].label = "heavy";
		all_pairs.setExperimentType("silac");
	}



	//--------------------------------------------------
	// filtering
  //--------------------------------------------------

  // exp must only contain level 1 spectra
  DoubleReal estimateMzSpacing(MSExperiment<Peak1D>& exp)
  {
    // estimate m/z step width
    UInt i = 0;
    while (i < exp.size() && exp[i].size() < 5) ++i;
    vector<Real> mz_spacing;
    for (Size j = 1; j < exp[i].size(); ++j)
    {
      mz_spacing.push_back(exp[i][j].getMZ() - exp[i][j-1].getMZ());
    }
    sort(mz_spacing.begin(), mz_spacing.end());
    return mz_spacing[mz_spacing.size() / 2];
  }

  vector<vector<DataPoint> > buildDataStructure(MSExperiment<Peak1D>& exp)
	{
		// extract level 1 spectra
		IntList levels=IntList::create("1");
		exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<MSExperiment<Peak1D>::SpectrumType>(levels, true)), exp.end());
    list<SILACFilter> filters;

    mz_stepwidth = estimateMzSpacing(exp);

		// create filters for all numbers of isotopes per peptide, charge states and mass shifts

		// iterate over all number for isotopes per peptide
    for (int isotopes_per_peptide = isotopes_per_peptide_min; isotopes_per_peptide <= isotopes_per_peptide_max; isotopes_per_peptide++)
		{
      // iterate over all charge states
			for (Int charge = charge_min; charge <= charge_max; charge++)
			{
				// iterate over all mass shifts
				for (unsigned i = 0; i < massShifts.size(); i++)
				{
					// convert vector<DoubleReal> to set<DoubleReal> for SILACFilter
          set<DoubleReal> massShifts_set;
          copy(massShifts[i].begin(), massShifts[i].end(), inserter(massShifts_set, massShifts_set.end()));
					filters.push_back(SILACFilter(massShifts_set, charge, model_deviation, isotopes_per_peptide));
				}
      }
		}

    if (in_filters == "")     // check if option "in_filters" is not specified
    {
      // create filtering
      SILACFiltering filtering(exp, mz_stepwidth, intensity_cutoff);
      filtering.setLogType(log_type_);

      // register filters to the filtering
      for (list<SILACFilter>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
      {
        filtering.addFilter(*filter_it);
      }


      // perform filtering
      filtering.filterDataPoints();

      // retrieve filtered data points
      // vector<vector<DataPoint> > data;
      for (list<SILACFilter>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
      {
 /*        set<DoubleReal> envelope_distances = filter_it->getEnvelopeDistances();

        // check if parameters for out_clusters are specified
        if (out_clusters_flag_2 == false || (abs(*envelope_distances.rbegin() - label_out_clusters_mass_shift) < numeric_limits<DoubleReal>::epsilon() && charge_out_clusters == filter_it->getCharge()))
*/
        data.push_back(filter_it->getElements());
      }
    }

    // delete experiment
		exp.clear(true);


    //--------------------------------------------------
    // store filter results from vector<vector<DataPoint> > data to .txt
    //--------------------------------------------------

    if (out_filters != "" && in_filters == "")     // check if option "out_filters" is specified and "in_filters" is not
    {
      ofstream outfile;
      outfile.open(out_filters.c_str());      // open ofstream to specified output file
      outfile << setprecision(16);      // set precision of outfile to 16 to avoid losing digits

      // vector of DataPoints
      outfile << "<VectorOfDataPoints>" << "\n";

      // iterate over outer DataPoint vector (vector<vector<DataPoint> >)
      for (vector<vector<DataPoint> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
      {
        // DataPoints
        outfile << "<DataPoints>" << "\n";

        // iterate over inner DataPoint vector (vector<DataPoint>)
        for (vector<DataPoint>::iterator it = data_it->begin(); it != data_it->end(); ++it)
        {
          // DataPoint
          outfile << "<DataPoint>" << "\n";

          // feature_id
          outfile << "<feature_id>" << "\n";
          outfile << it->feature_id << "\n";
          outfile << "</feature_id>" << "\n";

          // rt
          outfile << "<rt>" << "\n";
          outfile << it->rt << "\n";
          outfile << "</rt>" << "\n";

          // mz
          outfile << "<mz>" << "\n";
          outfile << it->mz << "\n";
          outfile << "</mz>" << "\n";

          // charge
          outfile << "<charge>" << "\n";
          outfile << it->charge << "\n";
          outfile << "</charge>" << "\n";

          // isotopes_per_peptide
          outfile << "<isotopes_per_peptide>" << "\n";
          outfile << it->isotopes_per_peptide << "\n";
          outfile << "</isotopes_per_peptide>" << "\n";

          // intensities
          outfile << "<intensities>" << "\n";

          // iterate over outer intensities vector (vector<vector<DoubleReal> >)
          for (vector<vector<DoubleReal> >::iterator intensities_it = it->intensities.begin(); intensities_it != it->intensities.end(); ++intensities_it)
          {
            outfile << "<intensities_" << intensities_it - it->intensities.begin() << ">\n";

            // // iterate over inner intensities vector (vector<DoubleReal>)
            for (vector<DoubleReal>::iterator intensity_it = intensities_it->begin(); intensity_it != intensities_it->end(); ++intensity_it)
            {
              outfile << *intensity_it << "\n";
            }
            outfile << "</intensities_" << intensities_it - it->intensities.begin() << ">\n";
          }
          outfile << "</intensities>" << "\n";

          // mass_shifts
          outfile << "<mass_shifts>" << "\n";

          // iterate over mass shifts vector (vector<DoubleReal>)
          for (vector<DoubleReal>::iterator mass_shifts_it = it->mass_shifts.begin(); mass_shifts_it != it->mass_shifts.end(); ++mass_shifts_it)
          {
            outfile << *mass_shifts_it << "\n";
          }
          outfile << "</mass_shifts>" << "\n";

          outfile << "</DataPoint>" << "\n";
        }
        outfile << "</DataPoints>" << "\n";
      }
      outfile << "</VectorOfDataPoints>" << "\n";

      outfile.close();      // close ofstream and store output file
    }


    //--------------------------------------------------
    // load filter results as vector<vector<DataPoint> > data from .txt
    //--------------------------------------------------

    if (in_filters != "")     // check if option "in_filters" is specified
    {
      vector<vector<DataPoint> > data_points_vector_in;
      vector<DataPoint> data_points_in;
      DataPoint data_point_in;

      vector<vector<DoubleReal> > intensities_vector_in;
      vector<DoubleReal> intensities_in;

      vector<DoubleReal> mass_shifts_in;

      ifstream infile;
      infile.open(in_filters.c_str());      // open ifstream from specified input file

      // check if infile can be opened
      if (!infile.is_open())
      {
        cout << "Error: could not open " << in_filters << "..." << endl;
      }
      else
      {
        String temp;

        // name cases and define states
        const int VECTOR_OF_DATA_POINTS_STATE = 0;
        const int DATA_POINTS_STATE = 1;
        const int DATA_POINT_STATE = 2;
        const int FEATURE_ID_STATE = 3;
        const int RT_STATE = 4;
        const int MZ_STATE = 5;
        const int CHARGE_STATE = 6;
        const int ISOTOPES_PER_PEPTIDE_STATE = 7;
        const int VECTOR_OF_INTENSITIES_STATE = 8;
        const int INTENSITIES_STATE = 9;
        const int MASS_SHIFTS_STATE = 10;
        const int END_STATE = 11;

        // set start state
        int state = VECTOR_OF_DATA_POINTS_STATE;

        while(state != END_STATE)
        {
          switch(state)
          {
            // case and state 0: vector of data points
        case VECTOR_OF_DATA_POINTS_STATE:
            if (infile.eof())
            {
              cout << "eof" << endl;
              state = END_STATE;
            } else
            {
              getline(infile, temp);
              if (temp != "")
              {
                state = DATA_POINTS_STATE;
              } else
              {
                cout << "end state" << endl;
                state = END_STATE;
              }
            }
            break;

            // case and state 1: data points
        case DATA_POINTS_STATE:
            getline (infile, temp);
            if (String(temp).hasPrefix("</"))
            {
              state = VECTOR_OF_DATA_POINTS_STATE;
            } else
            {
              state = DATA_POINT_STATE;
            }
            break;

            // case and state 2: data point
        case DATA_POINT_STATE:
            getline(infile, temp);
            if (String(temp).hasPrefix("</"))
            {
              data_points_vector_in.push_back(data_points_in);
              // cout << "size of data_points_vector_in: " << data_points_vector_in.size() << endl;
              data_points_in.clear();
              // cout << "size of data_points_in.clear(): " << data_points_in.size() << endl;
              state = DATA_POINTS_STATE;
            } else
            {
              data_point_in.intensities.clear();
              data_point_in.mass_shifts.clear();
              state = FEATURE_ID_STATE;
            }
            break;

            // case and state 3: feature_id
        case FEATURE_ID_STATE:
            getline(infile,temp);
            getline(infile, temp);
            data_point_in.feature_id = String(temp).toInt();
            getline (infile, temp);
            state = RT_STATE;
            break;

            // case and state 4: rt
        case RT_STATE:
            getline(infile, temp);
            getline(infile, temp);
            data_point_in.rt = String(temp).toDouble();
            getline (infile, temp);
            state = MZ_STATE;
            break;

            // case and state 5: mz
        case MZ_STATE:
            getline(infile, temp);
            getline(infile, temp);
            data_point_in.mz = String(temp).toDouble();
            getline(infile, temp);
            state = CHARGE_STATE;
            break;

            // case and state 6: charge
        case CHARGE_STATE:
            getline(infile, temp);
            getline(infile, temp);
            data_point_in.charge = String(temp).toInt();
            getline(infile, temp );
            state = ISOTOPES_PER_PEPTIDE_STATE;
            break;

            //case and state  7: isotopes_per_peptide
        case ISOTOPES_PER_PEPTIDE_STATE:
            getline(infile, temp);
            getline(infile, temp);
            data_point_in.isotopes_per_peptide = String(temp).toInt();
            getline(infile, temp);
            state = VECTOR_OF_INTENSITIES_STATE;
            break;

            // case and state 8: vector of intensities
        case VECTOR_OF_INTENSITIES_STATE:
            getline( infile, temp );
            if (!String(temp).hasPrefix("</"))
            {
              state = INTENSITIES_STATE;
            } else if (String(temp).hasPrefix("</"))
            {
              data_point_in.intensities.insert(data_point_in.intensities.end(), intensities_vector_in.begin(), intensities_vector_in.end());
              intensities_vector_in.clear();
              state = MASS_SHIFTS_STATE;
            }
            break;

            // case and state 9: intesities
        case INTENSITIES_STATE:
            getline( infile, temp );
            do
            {
              if(!String(temp).hasPrefix("<"))
              {
                DoubleReal intensity_in = String(temp).toDouble();
                intensities_in.push_back(intensity_in);
                // cout << "size of intensities_in: " << intensities_in.size() << endl;
              }
              getline( infile, temp );
            } while (!String(temp).hasPrefix("</"));
            intensities_vector_in.push_back(intensities_in);
            // cout << "size of intensities_vector_in: " << intensities_vector_in.size() << endl;
            intensities_in.clear();
            state = VECTOR_OF_INTENSITIES_STATE;
            break;

            // case and state 10: mass_shifts
        case MASS_SHIFTS_STATE:
            getline( infile, temp );
            do
            {
              if(!String(temp).hasPrefix("<"))
              {
                DoubleReal mass_shift_in = String(temp).toDouble();
                mass_shifts_in.push_back(mass_shift_in);
                // cout << "size of mass_shifts_in: " << mass_shifts_in.size() << endl;
              }
              getline( infile, temp );
            } while (!String(temp).hasPrefix("</"));
            data_point_in.mass_shifts.insert(data_point_in.mass_shifts.begin(), mass_shifts_in.begin(), mass_shifts_in.end());
            mass_shifts_in.clear();
            getline( infile, temp );      // temp steht auf </DataPoint>
            data_points_in.push_back(data_point_in);
            // cout << "size of data_points_in: " << data_points_in.size() << endl;
            state = DATA_POINT_STATE;
            break;

         default:
            break;
          }
        }
      }
      infile.close();


      //--------------------------------------------------
      // store read in filter results from .txt to .txt to enable comparison of original and read in filter results
      //--------------------------------------------------

      ofstream outfil;
      string out_filters_check = "check_" + in_filters;
      outfil.open(out_filters_check.c_str());     // open ofstream to specified output file
      outfil << setprecision(16);      // set precision of outfil to 16 to avoid losing digits

      outfil << "<VectorOfDataPoints>" << "\n";
      for (vector<vector<DataPoint> >::iterator data_it = data_points_vector_in.begin(); data_it != data_points_vector_in.end(); ++data_it)
      {
        outfil << "<DataPoints>" << "\n";

        for (vector<DataPoint>::iterator it = data_it->begin(); it != data_it->end(); ++it)
        {
          // DataPoint
          outfil << "<DataPoint>" << "\n";

          // feature_id
          outfil << "<feature_id>" << "\n";
          outfil << it->feature_id << "\n";
          outfil << "</feature_id>" << "\n";

          // rt
          outfil << "<rt>" << "\n";
          outfil << it->rt << "\n";
          outfil << "</rt>" << "\n";

          // mz
          outfil << "<mz>" << "\n";
          outfil << it->mz << "\n";
          outfil << "</mz>" << "\n";

          // charge
          outfil << "<charge>" << "\n";
          outfil << it->charge << "\n";
          outfil << "</charge>" << "\n";

          // isotopes_per_peptide
          outfil << "<isotopes_per_peptide>" << "\n";
          outfil << it->isotopes_per_peptide << "\n";
          outfil << "</isotopes_per_peptide>" << "\n";

          // intensities
          outfil << "<intensities>" << "\n";
          for (vector<vector<DoubleReal> >::iterator intensities_it = it->intensities.begin(); intensities_it != it->intensities.end(); ++intensities_it)
          {
            outfil << "<intensities_" << intensities_it - it->intensities.begin() << ">\n";
            for (vector<DoubleReal>::iterator intensity_it = intensities_it->begin(); intensity_it != intensities_it->end(); ++intensity_it)
            {
              outfil << *intensity_it << "\n";
            }
            outfil << "</intensities_" << intensities_it - it->intensities.begin() << ">\n";
          }
          outfil << "</intensities>" << "\n";

          // mass_shifts
          outfil << "<mass_shifts>" << "\n";
          for (vector<DoubleReal>::iterator mass_shifts_it = it->mass_shifts.begin(); mass_shifts_it != it->mass_shifts.end(); ++mass_shifts_it)
          {
            outfil << *mass_shifts_it << "\n";
          }
          outfil << "</mass_shifts>" << "\n";

          outfil << "</DataPoint>" << "\n";
        }
        outfil << "</DataPoints>" << "\n";
      }
      outfil << "</VectorOfDataPoints>" << "\n";

      outfil.close();     // close ofstream

      // set loaded filter results from .txt as input for clustering
      data = data_points_vector_in;
    }

    return data;
	}


	ExitCodes main_(int , const char**)
	{		
		handleParameters();

		//--------------------------------------------------
		// loading input from .mzML
		//--------------------------------------------------

		MzMLFile file;
		MSExperiment<Peak1D> exp;

		file.setLogType(log_type_);
    file.load(in, exp);

		// set size of input map
		exp.updateRanges();

		all_pairs.getFileDescriptions()[0].size = exp.getSize();
		all_pairs.getFileDescriptions()[1].size = exp.getSize();
		all_pairs.getFileDescriptions()[2].size = exp.getSize();


		//--------------------------------------------------
		// build SILACData structure
		//--------------------------------------------------

    vector<vector<DataPoint> > data=buildDataStructure(exp);


		//--------------------------------------------------
		// clustering
		//--------------------------------------------------

    vector<vector<Real> > silhouettes;
    vector<Cluster> clusters;
    vector<Tree> subtrees;

    for (vector<vector<DataPoint> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
		{

cout << "size of vector<vector<DataPoint>>: " << data.size() << std::endl;

      // check if there are at least two points for clustering
      if (data_it->size() >= 2)
      {

cout << "size of vector<DataPoint>: " << data_it->size() << endl;

        // hierarchical clustering
        CentroidLinkage method(rt_scaling);
        HashClustering c(*data_it, rt_threshold, mz_threshold, method);
        c.setLogType(log_type_);
        c.performClustering();
        vector<Tree> current_subtrees;
        c.getSubtrees(current_subtrees);
        subtrees.insert(subtrees.end(), current_subtrees.begin(), current_subtrees.end());
        vector<Cluster> current_clusters;
        c.createClusters(current_clusters);
        const vector<vector<Real> >& current_silhouettes = c.getSilhouetteValues();
        silhouettes.insert(silhouettes.end(), current_silhouettes.begin(), current_silhouettes.end());

        // QT clustering
/*        DoubleReal isotope_distance = 1.000495 / (DoubleReal)data_it->front().charge;
        QTClustering c(*data_it, rt_threshold, mz_threshold, isotope_distance);
        c.setLogType(log_type_);
        vector<Cluster> current_clusters =  c.performClustering();
*/
        clusters.insert(clusters.end(), current_clusters.begin(), current_clusters.end());
      }
    }

    sort(clusters.begin(), clusters.end(), clusterCompare);



		//--------------------------------------------------
		// subtree output (for debug)
		//--------------------------------------------------

		if (getFlag_("silac_debug"))
		{
      vector<String> colors;
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

			Size subtree_number = 1;
      for (vector<vector<BinaryTreeNode> >::iterator subtree_it = subtrees.begin(); subtree_it != subtrees.end(); ++subtree_it)
			{
        set<DataPoint*> leafs;
        for (vector<BinaryTreeNode>::iterator tree_it = subtree_it->begin(); tree_it != subtree_it->end(); ++tree_it)
				{
					leafs.insert(tree_it->data1);
					leafs.insert(tree_it->data2);
				}
        for (set<DataPoint*>::iterator leafs_it = leafs.begin();leafs_it != leafs.end(); ++leafs_it)
				{
					Feature tree_point;
					tree_point.setRT((*leafs_it)->rt);
					tree_point.setMZ((*leafs_it)->mz);
					tree_point.setIntensity((*leafs_it)->intensities[0][0]);
					tree_point.setCharge((*leafs_it)->charge);
					tree_point.setMetaValue("subtree", subtree_number);
					tree_point.setMetaValue("color", colors[subtree_number%colors.size()]);
					subtree_points.push_back(tree_point);
				}
				++subtree_number;
			}

			// required, as somehow the order of features on some datasets between Win & Linux is different and thus the TOPPtest might fail
			subtree_points.sortByPosition();
		}



		//--------------------------------------------------
		// create a map of SILAC type names
		//--------------------------------------------------

    map<Size,String> silac_types;
    silac_types.insert(make_pair(0, "Singlet"));
    silac_types.insert(make_pair(1, "Doublet"));
    silac_types.insert(make_pair(2, "Triplet"));
    silac_types.insert(make_pair(3, "Quadruplet"));



		//--------------------------------------------------------------
		// determine file name for debug output
		//--------------------------------------------------------------
		String debug_trunk = in;
		if (in.has('.'))
		{
			debug_trunk = in.substr(0, in.find_first_of('.'));
		}



		//--------------------------------------------------
    // mzML output (out_filters_mzML)
		//--------------------------------------------------

    if (out_filters_mzML != "")
    {
      // build new experiment to store only points that passed the filters as .mzML
      for (vector<vector<DataPoint> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
      {
        for (vector<DataPoint>::iterator it = data_it->begin(); it != data_it->end(); ++it)
        {
          MSExperiment<Peak1D>::SpectrumType spectrum;
          spectrum.setMSLevel(1);
          spectrum.setRT(it->rt);
          MSExperiment<Peak1D>::PeakType peak;
          peak.setIntensity(it->intensities[0][0]);
          peak.setMZ(it->mz);
          spectrum.push_back(peak);
          filter_exp.push_back(spectrum);
        }
      }
    }



		//--------------------------------------------------
		// consensusXML output
		//--------------------------------------------------

		if (out != "")
		{
			// write ratios of all cluster to additional *.dat
			String debug_dat = debug_trunk + ".dat";
      ofstream stream_ratios(debug_dat.c_str());
			Size id = 0;

			// iterate over all clusters
      for (vector<Cluster>::iterator cluster_it = clusters.begin(); cluster_it != clusters.end(); ++cluster_it)
			{
				DoubleReal rt = 0.0;
				DoubleReal mz = 0.0;
				DoubleReal total_intensity = 0.0;
				Size mass_shifts_size = (*(cluster_it->begin()))->mass_shifts.size();
				// create a vector with the maximum intensity of each isotope peak
        vector<DoubleReal> max_intensities(mass_shifts_size,0.0);
				// specify SILAC type
				String silac_type = silac_types[mass_shifts_size];
				Int charge = (*(cluster_it->begin()))->charge;

				// add mass shifts as value for Quality
				String mass_shift = "";
				String mass_shift_2 = "";
				DoubleReal mass_shift_final;

        for (vector<DoubleReal>::iterator shift_it = (*(cluster_it->begin()))->mass_shifts.begin(); shift_it != (*(cluster_it->begin()))->mass_shifts.end(); ++shift_it)
				{							
					DoubleReal mass_shift_current = *shift_it * charge;			// convert mass shift from Th to Da
					mass_shift += "0" + (String)mass_shift_current;			// combine mass shifts as string
				}

				// format string of mass shifts
				for (int i = 0; i < 2; ++i)
				{
					int found = mass_shift.find(".");
					if (found > 0)
					{
						mass_shift_2 += mass_shift.substr(found - 2, 3);
						mass_shift.erase(found, 1);
							
						found = mass_shift_2.find(".");
						mass_shift_2.erase(found, 1);
					}
        }
						
				mass_shift_2.insert(2, ".");			// insert "." between the two mass shifts

				// convert mass shits from string to DoubleReal (value for Quality has to be of type DoubleReal)
        mass_shift_final = String(mass_shift_2).toDouble();

				// intensity vector used for linear regression
        vector<vector<DoubleReal> > intensities(mass_shifts_size);
				// iterate over the cluster elements
				for (Cluster::iterator el_it = cluster_it->begin(); el_it != cluster_it->end(); ++el_it)
				{
					// extract all SILAC pattern intensities of the current element
          vector<vector<DoubleReal> >& element_intensities = (*el_it)->intensities;
					// add monoisotopic intensity of the light peak to the total intensity
					total_intensity += element_intensities[0][0];
					// add monoisotopic RT position of the light peak to the total RT value, weighted by the intensity
					rt += element_intensities[0][0] * (*el_it)->rt;
					// iterate over all mass shifts of the element
					for (Size i = 0; i < mass_shifts_size; ++i)
					{
						// create a temporary vector with all isotope pattern intensities of the current mass shift
            vector<DoubleReal>& current_intensities=element_intensities[i];
            intensities[i].insert(intensities[i].end(),current_intensities.begin(),current_intensities.end());

						// find maximum intensity
            vector<DoubleReal>::iterator max_position = max_element(current_intensities.begin(), current_intensities.end());
						if (*max_position > max_intensities[i])
						{
							max_intensities[i] = *max_position;
							mz = (*el_it)->mz;
						}
					}
				}

				// average retention time
				rt /= total_intensity;

				stream_ratios << id << "\t" << cluster_it->size() << "\t" << rt << "\t" << mz;

/*				for (Size k = 1;	k < mass_shifts_size; ++k)
				{
					// perform linear regression for each mass shift != 0
					Math::LinearRegression linear_reg;
					linear_reg.computeRegressionNoIntercept(0.95, intensities[0].begin(), intensities[0].end(), intensities[k].begin());
					stream_ratios  << "\t" << linear_reg.getSlope();

					// create consensus feature for each mass shift !=0
					ConsensusFeature consensus_feature;
					consensus_feature.setRT(rt);
					consensus_feature.setMZ(mz);
					consensus_feature.setIntensity(linear_reg.getSlope());
					consensus_feature.setCharge(charge);
					consensus_feature.setQuality(linear_reg.getRSquared());

					// insert feature handle for each mass shift
					for (Size l = 0; l < mass_shifts_size; ++l)
					{
						FeatureHandle handle;
						handle.setRT(rt);
						handle.setMZ(mz+(*(cluster_it->begin()))->mass_shifts[l]);
						handle.setIntensity(max_intensities[l]);
						handle.setCharge(charge);
						handle.setMapIndex(l);
						handle.setUniqueId(id);
						consensus_feature.insert(handle);
					}
					all_pairs.push_back(consensus_feature);
				}
*/

					// create consensus feature for each mass shift !=0
					ConsensusFeature consensus_feature;
					consensus_feature.setRT(rt);
					consensus_feature.setMZ(mz);
					consensus_feature.setIntensity(max_intensities[0]);			// set intensity of light peptiide to intensity of consensus
					consensus_feature.setCharge(charge);
          consensus_feature.setQuality(mass_shift_final);			// set mass shifts as value for Quality (format: Da.Da)

					// insert feature handle for each mass shift
					for (Size l = 0; l < mass_shifts_size; ++l)
					{

						// perform linear regression for each mass shift != 0

						Math::LinearRegression linear_reg;
						linear_reg.computeRegressionNoIntercept(0.95, intensities[0].begin(), intensities[0].end(), intensities[l].begin());
						stream_ratios  << "\t" << linear_reg.getSlope();

						FeatureHandle handle;
						handle.setRT(rt);
						handle.setMZ(mz+(*(cluster_it->begin()))->mass_shifts[l]);
						handle.setIntensity(linear_reg.getSlope());			// set ratios as values for intensities: it(map="0") = l/l, it(map="1") = m/l, it(map="2") = h/l
						handle.setCharge(charge);
						handle.setMapIndex(l);
						handle.setUniqueId(id);
						consensus_feature.insert(handle);
					}
					all_pairs.push_back(consensus_feature);

        stream_ratios << endl;
				++id;
			}
			stream_ratios.close();
		}



		//--------------------------------------------------
		// featureXML output (out_clusters)
		//--------------------------------------------------

		if (out_clusters != "")
		{
      vector<String> colors;
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

      for (vector<vector<DataPoint> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
			{
        for (vector<DataPoint>::iterator it = data_it->begin(); it != data_it->end(); ++it)
        {
          // visualize the light variant
          Feature cluster_point;
          cluster_point.setRT(it->rt);
          cluster_point.setMZ(it->mz);
          cluster_point.setIntensity(it->intensities[0][0]);
          cluster_point.setCharge(it->charge);
          //cluster_point.setOverallQuality(it->quality);
          cluster_point.setQuality(0, it->quality);
          cluster_point.setMetaValue("SILAC type", silac_types[it->mass_shifts.size() - 1]);
						
          // add mass shifts as meta value and as value for OverallQuality
          Int charge = it->charge;
          String mass_shift_meta_value = "";
          String mass_shift = "";
          String mass_shift_2 = "";
          DoubleReal mass_shift_final;

          for (vector<DoubleReal>::iterator shift_it = it->mass_shifts.begin() + 1; shift_it != it->mass_shifts.end(); ++shift_it)
          {
            mass_shift_meta_value += ((String)*shift_it) + " ";			// mass shifts as meta value
							
            DoubleReal mass_shift_current = *shift_it * charge;			// convert mass shift from Th to Da
            mass_shift += "0" + (String)mass_shift_current;			// combine mass shifts as string
          }
						
          // format string of mass shifts
          for (int i = 0; i < 2; ++i)
          {
            int found = mass_shift.find(".");
            if (found > 0)
            {
              mass_shift_2 += mass_shift.substr(found - 2, 3);
              mass_shift.erase(found, 1);
								
              found = mass_shift_2.find(".");
              mass_shift_2.erase(found, 1);
            }
          }
						
          mass_shift_2.insert(2, ".");			// insert "." between the two mass shifts

          // convert mass shits from string to DoubleReal (value for OverallQuality has to be of type DoubleReal)
          mass_shift_final = String(mass_shift_2).toDouble();
						
          cluster_point.setOverallQuality(mass_shift_final);			// set mass shifts as value for OverallQuality (format: Da.Da)
						
          if (mass_shift_meta_value != "" && it->mass_shifts.size() - 1 == 1)
          {
            cluster_point.setMetaValue("Mass shift (l/h)", mass_shift_meta_value);
          }
          else if (mass_shift_meta_value != "" && it->mass_shifts.size() - 1 == 2)
          {
            cluster_point.setMetaValue("Mass shift (l/m l/h)", mass_shift_meta_value);
          }
          else if (mass_shift_meta_value != "" )
          {
            cluster_point.setMetaValue("Mass shift", mass_shift_meta_value);
          }

          cluster_point.setMetaValue("Cluster id", it->cluster_id);
          cluster_point.setMetaValue("color", colors[it->cluster_id%colors.size()]);
          cluster_point.setMetaValue("isotopes per peptide", it->isotopes_per_peptide);
						
          if (getFlag_("silac_debug"))
          {
            cluster_point.setMetaValue("Cluster size", it->cluster_size);
            cluster_point.setMetaValue("feature_id", it->feature_id);
          }

          // if no parameters for out_cluster are sprcified write all clusters
//					if (out_clusters_flag == false)
              all_cluster_points.push_back(cluster_point);

          // if parameters for out_clusters are specified write only coresponding clusters
//					else if (it->charge == charge_out_clusters && ms_final == ms_final_2)
//						all_cluster_points.push_back(cluster_point);
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
			// names of dat files
			String debug_clusters_dat = debug_trunk + "_cluster_sizes.dat";

			// write all cluster data points to *_clusters.dat
      ofstream stream_clusters(debug_clusters_dat.c_str());
      for(vector<Cluster>::iterator cluster_it = clusters.begin(); cluster_it != clusters.end(); ++cluster_it)
			{
        stream_clusters << cluster_it->size() << endl;
			}
		}



		//--------------------------------------------------------------
		// write output
		//--------------------------------------------------------------

		// consensusXML
		if (out != "")
		{
			// assign unique ids
			all_pairs.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			// annotate output with data processing info
			addDataProcessing_(all_pairs, getProcessingInfo_(DataProcessing::QUANTITATION));

			ConsensusXMLFile c_file;
			c_file.store(out, all_pairs);
		}

    // mzML
    if (out_filters_mzML != "")
		{
      MzMLFile m_file;
      m_file.store(out_filters_mzML, filter_exp);
		}

		// featureXML
		if (out_clusters != "")
		{
			// assign unique ids
			all_cluster_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			FeatureXMLFile f_file;
			f_file.store(out_clusters, all_cluster_points);
		}

		// debug output
		if (getFlag_("silac_debug"))
		{
			// assign unique ids
			subtree_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);

			FeatureXMLFile t_file;
			t_file.store(debug_trunk+ "_subtrees.featureXML", subtree_points);
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

