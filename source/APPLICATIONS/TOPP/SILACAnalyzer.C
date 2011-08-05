// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Lars Nilse, Steffen Sass, Holger Plattfaut, Bastian Blank $
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
#include <OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace OpenMS;
using namespace std;

typedef vector<SILACPattern*> Cluster;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SILACAnalyzer SILACAnalyzer

  @brief Identifies peptide pairs in LC-MS data and determines their relative abundance.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ SILACAnalyzer \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
    </tr>
  </table>
</CENTER>

  SILACAnalyzer is a tool for the fully automated analysis of quantitative proteomics data. It identifies pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we first explain the algorithm and then discuss the tuning of its parameters.

  <b>Algorithm</b>

  The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).

  We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
  - all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold,
  - the intensity profiles in neighbourhoods around all six m/z positions show a good correlation and
  - the relative intensity ratios within a peptide agree up to a factor with the ratios of a theoretic averagine model.

  Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster corresponds to the mono-isotopic mass trace of the lightest peptide of a SILAC pattern. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

  @image html SILACAnalyzer_algorithm.png

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SILACAnalyzer.cli

  <b>Parameter Tuning</b>

  SILACAnalyzer can detect SILAC patterns of any number of peptides, i.e. doublets (pairs), triplets, quadruplets et cetera.

  <i>input:</i>
  - in [*.mzML] - LC-MS dataset to be analyzed
  - ini [*.ini] - file containing all parameters (see discussion below)

  <i>standard output:</i>
  - out [*.consensusXML] - contains the list of identified peptides (retention time and m/z of the lightest peptide, ratios)

  <i>optional output:</i>
  - out_clusters [*.consensusXML] - contains the complete set of data points passing the filters, see Fig. (e)

  The results of an analysis can easily visualized within TOPPView. Simply load *.consensusXML and *.featureXML as layers over the original *.mzML.

  Parameters in section <i>algorithm:</i>
  - <i>allow_missing_peaks</i> - Low intensity peaks might be missing from the isotopic pattern of some of the peptides. Specify if such peptides should be included in the analysis.
  - <i>mz_threshold</i> - Upper bound for the width [Th] of an isotopic peak.
  - <i>rt_threshold</i> - Upper bound for the retention time [s] over which a characteristic peptide elutes.
  - <i>intensity_cutoff</i> - Lower bound for the intensity of isotopic peaks in a SILAC pattern.
  - <i>intensity_correlation</i> - Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.
  - <i>model_deviation</i> - Upper bound on the factor by which the ratios of observed isotopic peaks are allowed to differ from the ratios of the theoretic averagine model, i.e. ( theoretic_ratio / model_deviation ) < observed_ratio < ( theoretic_ratio * model_deviation ).

  Parameters in section <i>sample:</i>
  - <i>labels</i> - Labels used for labelling the sample. [...] specifies the labels for a single sample. For example, [Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. For permitted labels see section <i>labels</i>.
  - <i>charge</i> - Range of charge states in the sample, i.e. min charge : max charge.
  - <i>missed_cleavages</i> - Maximum number of missed cleavages.
  - <i>peaks_per_peptide</i> - Range of peaks per peptide in the sample, i.e. min peaks per peptide : max peaks per peptide.

 Parameters in section <i>labels:</i>
 This section contains a list of all isotopic labels currently available for analysis of SILAC data with SILACAnalyzer.
 
 <b>References:</b>
  @n L. Nilse, M. Sturm, D. Trudgian, M. Salek, P. Sims, K. Carroll, S. Hubbard,  <a href="http://www.springerlink.com/content/u40057754100v71t">SILACAnalyzer - a tool for differential quantitation of stable isotope derived data</a>, in F. Masulli, L. Peterson, and R. Tagliaferri (Eds.): CIBB 2009, LNBI 6160, pp. 4555, 2010.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

bool clusterCompare( Cluster v1, Cluster v2)
{
   return v1.size() > v2.size();
}

class TOPPSILACAnalyzer
: public TOPPBase
{
  private:

    // input and output files
    String in;
    String out;
    String out_clusters;    

    String out_filters;
    String in_filters;
    String out_debug;

    // section "sample"
    String selected_labels;
    Int charge_min;
    Int charge_max;
    Int missed_cleavages;
    Int isotopes_per_peptide_min;
    Int isotopes_per_peptide_max;

    // section "algorithm"
    DoubleReal mz_threshold;
    DoubleReal rt_threshold;
    DoubleReal intensity_cutoff;
    DoubleReal intensity_correlation;
    DoubleReal model_deviation;
    bool allow_missing_peaks;

    // section "labels"
    map<String, DoubleReal> label_identifiers;
    vector<vector <String> > SILAClabels;     // list of SILAC labels, e.g. selected_labels="[Lys4,Arg6][Lys8,Arg10]" => SILAClabels[0][1]="Arg6"
    vector<vector <DoubleReal> > massShifts;      // list of mass shifts

    typedef HierarchicalClustering<SILACPattern *> Clustering;

    vector<vector<SILACPattern> > data;
    vector<Clustering *> cluster_data;
    FeatureMap<> subtree_points;
    MSExperiment<Peak1D> filter_exp;


  public:
    TOPPSILACAnalyzer()
    : TOPPBase("SILACAnalyzer","Determination of peak ratios in LC-MS data",true)
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
    registerOutputFile_("out_clusters", "<file>", "", "Optional output file containing data points passing all filters, hence belonging to a SILAC pattern. Points of the same colour correspond to the mono-isotopic peak of the lightest peptide in a pattern.", false, true);
    setValidFormats_("out_clusters", StringList::create("consensusXML"));

    // create optional flag for additional output file (.consensusXML) to store filter results
    registerOutputFile_("out_filters", "<file>", "", "Additional output file containing all points that passed the filters as txt. Suitable as input for \"in_filters\" to perform clustering without preceding filtering process.", false, true);
    setValidFormats_("out_filters", StringList::create("consensusXML"));
    // create optional flag for additional input file (.consensusXML) to load filter results
    registerOutputFile_("in_filters", "<file>", "", "Additional input file containing all points that passed the filters as txt. Use output from \"out_filters\" to perform clustering only.", false, true);
    setValidFormats_("in_filters", StringList::create("consensusXML"));
    registerStringOption_("out_debug", "<filebase>", "", "Filename base for debug output.", false, true);

    // create section "labels" for adjusting masses of labels
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'sample\'.");
    // create section "sample" for adjusting sample parameters
    registerSubsection_("sample", "Parameter describing the sample and its labels.");
    // create section "algorithm" for adjusting algorithm parameters
    registerSubsection_("algorithm", "Parameters for the algorithm.");

    // create flag for missing peaks
    registerFlag_("algorithm:allow_missing_peaks", "Low intensity peaks might be missing from the isotopic pattern of some of the peptides. Should such peptides be included in the analysis?", true);
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
      defaults.setValue("Arg6", 6.0201290268, "Arg6 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Arg6", 0.0);
      defaults.setValue("Arg10", 10.008268600, "Arg10 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Arg10", 0.0);
      defaults.setValue("Lys4", 4.0251069836, "Lys4 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Lys4", 0.0);
      defaults.setValue("Lys6", 6.0201290268, "Lys6 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Lys6", 0.0);
      defaults.setValue("Lys8", 8.0141988132, "Lys8 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Lys8", 0.0);
      defaults.setValue("dICPL4", 4.025107, "mass difference between isotope-coded protein labels ICPL 4 and ICPL 0", StringList::create("advanced"));
      defaults.setMinFloat("dICPL4", 0.0);
      defaults.setValue("dICPL6", 6.020129, "mass difference between isotope-coded protein labels ICPL 6 and ICPL 0", StringList::create("advanced"));
      defaults.setMinFloat("dICPL6", 0.0);
      defaults.setValue("dICPL10", 10.045236, "mass difference between isotope-coded protein labels ICPL 10 and ICPL 0", StringList::create("advanced"));
      defaults.setMinFloat("dICPL10", 0.0);
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
      defaults.setValue("labels", "[Arg6]", "Labels used for labelling the sample. [...] specifies the labels for a single sample. For example, [Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. For permitted labels see \'advanced parameters\', section \'labels\'.");
      defaults.setValue("charge", "2:3", "Range of charge states in the sample, i.e. min charge : max charge.");
      defaults.setValue("missed_cleavages", 0 , "Maximum number of missed cleavages.");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("peaks_per_peptide", "3:4", "Range of peaks per peptide in the sample, i.e. min peaks per peptide : max peaks per peptide.", StringList::create("advanced"));
    }


    //--------------------------------------------------
    // section algorithm
    //--------------------------------------------------

    if (section == "algorithm")
    {
      defaults.setValue("mz_threshold", 0.1, "Upper bound for the width [Th] of an isotopic peak.");
      defaults.setMinFloat("mz_threshold", 0.0);
      defaults.setValue("rt_threshold", 50.0, "Upper bound for the retention time [s] over which a characteristic peptide elutes. ");
      defaults.setMinFloat("rt_threshold", 0.0);
      defaults.setValue("intensity_cutoff", 10000.0, "Lower bound for the intensity of isotopic peaks in a SILAC pattern.");
      defaults.setMinFloat("intensity_cutoff", 10.0);
      defaults.setValue("intensity_correlation", 0.9, "Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.");
      defaults.setMinFloat("intensity_correlation", 0.0);
      defaults.setMaxFloat("intensity_correlation", 1.0);
      defaults.setValue("model_deviation", 6.0, "Upper bound on the factor by which the ratios of observed isotopic peaks are allowed to differ from the ratios of the theoretic averagine model, i.e. ( theoretic_ratio / model_deviation ) < observed_ratio < ( theoretic_ratio * model_deviation ).");
      defaults.setMinFloat("model_deviation", 1.0);
    }

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
    // get name of additional clusters output file (.consensusXML)
    out_clusters = getStringOption_("out_clusters");

    // get name of additional filters output file (.consensusXML)
    out_filters = getStringOption_("out_filters");
    // get name of additional filters input file (.consensusXML)
    in_filters = getStringOption_("in_filters");
    out_debug = getStringOption_("out_debug");


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
    label_identifiers.insert(make_pair("dICPL4", getParam_().getValue("labels:dICPL4")));
    label_identifiers.insert(make_pair("dICPL6", getParam_().getValue("labels:dICPL6")));
    label_identifiers.insert(make_pair("dICPL10", getParam_().getValue("labels:dICPL10")));

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
    map<String,DoubleReal>::iterator dicpl4 = label_identifiers.find("dICPL4");
    map<String,DoubleReal>::iterator dicpl6 = label_identifiers.find("dICPL6");
    map<String,DoubleReal>::iterator dicpl10 = label_identifiers.find("dICPL10");

    // create string of all labels from advanced section "labels"
    String labels = "Arg6 Arg10 Lys4 Lys6 Lys8 Methyl4 Methyl8 Methyl12 Methyl16 Methyl24 Methyl32 dICPL4 dICPL6 dICPL10";


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

    // check if charge_min is smaller than charge max, if not swap
    if (charge_min > charge_max)
      swap(charge_min, charge_max);

    // get selected peaks range
    String isotopes_per_peptide_string = getParam_().getValue("sample:peaks_per_peptide");
    DoubleReal isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
    parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min = (Int)isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max = (Int)isotopes_per_peptide_max_temp;

    //check if isotopes_per_peptide_min is smaller than isotopes_per_peptide_max, if not swap
    if (isotopes_per_peptide_min > isotopes_per_peptide_max)
      swap(isotopes_per_peptide_min, isotopes_per_peptide_max);


    //--------------------------------------------------
    // section algorithm
    //--------------------------------------------------

    // get selected mz_threshold
    mz_threshold = getParam_().getValue("algorithm:mz_threshold");
    // get selected rt_threshold
    rt_threshold = getParam_().getValue("algorithm:rt_threshold");
    // get selected intensity_cutoff
    intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
    // get selected intensity_correlation
    intensity_correlation = getParam_().getValue("algorithm:intensity_correlation");
    // get selected model_deviation
    model_deviation = getParam_().getValue("algorithm:model_deviation");
    // get flag for missing peaks
    allow_missing_peaks = getFlag_("algorithm:allow_missing_peaks");


    //--------------------------------------------------
    // calculate all possible mass shifts for labelets from section "sample:labels" (concernig missed_cleavage)
    //--------------------------------------------------

    // split string of SILAC labels (selected_labels) and save in a list (SILAClabels)
    vector<String> tempList; // temporary list of strings for SILAC labelets, e.g. "Lys6,Arg8"
    boost::split( tempList, selected_labels, boost::is_any_of("[](){}") ); // any bracket allowed to separate labelets
    for (UInt i = 0; i < tempList.size(); i++)
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
    for (UInt i = 0; i < SILAClabels.size(); i++)
    {
      cout << "SILAC label " << i + 1 << ":   ";
      for (UInt j = 0; j < SILAClabels[i].size(); j++)
      {
        cout << SILAClabels[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;

    // check if all selected labels are included in advanced section "labels"
    for (UInt i = 0; i < SILAClabels.size(); i++)
    {
      for (UInt j = 0; j < SILAClabels[i].size(); ++j)
      {
        Int found = (Int) labels.find(SILAClabels[i][j]);

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
          for (Int dICPLPerPeptide = 0; dICPLPerPeptide <= missed_cleavages + 1; dICPLPerPeptide++)
          {
            if ( ArgPerPeptide + LysPerPeptide + MethylPerPeptide + dICPLPerPeptide > 0 && ArgPerPeptide + LysPerPeptide + MethylPerPeptide + dICPLPerPeptide <= missed_cleavages + 1 )
            {
              vector<DoubleReal> massShiftVector;
              for (UInt i = 0; i < SILAClabels.size(); i++)
              {
                DoubleReal massShift = 0;
                // Considering the case of an amino acid (e.g. LysPerPeptide != 0) for which no label is present (e.g. Lys4There + Lys8There == 0) makes no sense. Therefore each amino acid will have to give its "Go Ahead" before the shift is calculated.
                bool goAhead_Lys = false;
                bool goAhead_Arg = false;
                bool goAhead_Methyl = false;
                bool goAhead_dICPL = false;

                for (UInt j = 0; j < SILAClabels[i].size(); j++)
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
                  Int dICPL4There = 0;
                  Int dICPL6There = 0;
                  Int dICPL10There = 0;

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
                  if ( SILAClabels[i][j].find("dICPL4") == 0 ) dICPL4There = 1;
                  if ( SILAClabels[i][j].find("dICPL6") == 0 ) dICPL6There = 1;
                  if ( SILAClabels[i][j].find("dICPL10") == 0 ) dICPL10There = 1;

                  goAhead_Arg = goAhead_Arg || !( (ArgPerPeptide != 0 && Arg6There + Arg10There == 0) );
                  goAhead_Lys = goAhead_Lys || !( (LysPerPeptide != 0 && Lys4There + Lys6There + Lys8There == 0) );
                  goAhead_Methyl = goAhead_Methyl || !( (MethylPerPeptide != 0 && Methyl4There + Methyl8There + Methyl12There + Methyl16There + Methyl24There + Methyl32There == 0) );
                  goAhead_dICPL = goAhead_dICPL || !( (dICPLPerPeptide != 0 && dICPL4There + dICPL6There + dICPL10There == 0) );

                  massShift = massShift + ArgPerPeptide * ( Arg6There * (arg6->second) + Arg10There * (arg10->second) ) + LysPerPeptide * ( Lys4There * (lys4->second) + Lys6There * (lys6->second) + Lys8There * (lys8->second) ) +  MethylPerPeptide * (Methyl4There * (methyl4->second) + Methyl8There * (methyl8->second) + Methyl12There * (methyl12->second) + Methyl16There * (methyl16->second) + Methyl24There * (methyl24->second) + Methyl32There * (methyl32->second) ) + dICPLPerPeptide * ( dICPL4There * (dicpl4->second) + dICPL6There * (dicpl6->second) + dICPL10There * (dicpl10->second) );
                }

                if (goAhead_Arg && goAhead_Lys && goAhead_Methyl && goAhead_dICPL)
                  massShiftVector.push_back(massShift);
              }

              if (!massShiftVector.empty())
                massShifts.push_back(massShiftVector);
            }
          }
        }
      }
    }

    // create zero-mass-shift to search for peptides if no label is specified
    if (massShifts.size() == 0)
    {
      vector<DoubleReal> mass_shift_vector_peptide(1, 0.0);
      massShifts.push_back(mass_shift_vector_peptide);
    }

    // sort the mass shift vector
    sort(massShifts.begin(), massShifts.end());

    // print mass shifts
    for (UInt i = 0; i < massShifts.size(); i++)
    {
      cout << "mass shift " << i + 1 << ":   ";
      for (UInt j = 0; j < massShifts[i].size(); j++)
      {
        cout << massShifts[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }


  //--------------------------------------------------
  // filtering
  //--------------------------------------------------

  void filterData(MSExperiment<Peak1D>& exp)
  {
    // extract level 1 spectra
    IntList levels=IntList::create("1");
    exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<MSExperiment<Peak1D>::SpectrumType>(levels, true)), exp.end());
    list<SILACFilter> filters;

    // create filters for all numbers of isotopes per peptide, charge states and mass shifts
    // iterate over all number for peaks per peptide (from max to min)
    for (Int isotopes_per_peptide = isotopes_per_peptide_max; isotopes_per_peptide >= isotopes_per_peptide_min; isotopes_per_peptide--)
    {
      // iterate over all charge states (from max to min)
      for (Int charge = charge_max; charge >= charge_min; charge--)
      {
        // iterate over all mass shifts
        for (UInt i = 0; i < massShifts.size(); i++)
        {
          // convert vector<DoubleReal> to set<DoubleReal> for SILACFilter
          vector<DoubleReal> massShifts_set = massShifts[i];

          //copy(massShifts[i].begin(), massShifts[i].end(), inserter(massShifts_set, massShifts_set.end()));
          filters.push_back(SILACFilter(massShifts_set, charge, model_deviation, isotopes_per_peptide));
        }
      }
    }

    if (in_filters == "")     // check if option "in_filters" is not specified
    {
      // create filtering
      SILACFiltering filtering(exp, intensity_cutoff, intensity_correlation, allow_missing_peaks);
      filtering.setLogType(log_type_);

      // register filters to the filtering
      for (list<SILACFilter>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
      {
        filtering.addFilter(*filter_it);
      }

      // perform filtering
      filtering.filterDataPoints();

      // retrieve filtered data points
      for (list<SILACFilter>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
      {
        data.push_back(filter_it->getElements());
      }
    }

    // delete experiment
    exp.clear(true);


    //--------------------------------------------------
    // combine DataPoints to improve the clustering
    //--------------------------------------------------

    // DataPoints that originate from filters with same charge state and mass shift(s)
    // and whose filters only differ in number of isotopes per peptide are combined
    // to get one cluster for peptides whose elution profile varies in number of isotopes per peptide

    // perform combination only if the user specified a peaks_per_peptide range > 1
    if (isotopes_per_peptide_min != isotopes_per_peptide_max)
    {
      // erase empty filter results from "data"
      vector<vector<SILACPattern> > data_temp;

      for (vector<vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
      {
        if (data_it->size() != 0)
        {
          data_temp.push_back(*data_it);     // keep DataPoint if it is not empty
        }
      }

      data.swap(data_temp);     // data = data_temp
      data_temp.clear();      // clear "data_temp"

      if (data.size() >= 2)
      {
        Int temp = 0;
        // combine corresponding DataPoints
        vector<vector<SILACPattern> >::iterator data_it_1 = data.begin();      // first iterator over "data" to get first DataPoint for combining
        vector<vector<SILACPattern> >::iterator data_it_2 = data_it_1 + 1;     // second iterator over "data" to get second DataPoint for combining
        vector<vector<SILACPattern> >::iterator data_it_end = data.end() - 1;      // pointer to second last elemnt of "data"
        vector<SILACPattern>::iterator it_1;     // first inner iterator over elements of first DataPoint
        vector<SILACPattern>::iterator it_2;     // second inner iterator over elements of second DataPoint

        while (data_it_1 < data_it_end)      // check for combining as long as first DataPoint is not second last elment of "data"
        {          
          while (data_it_1->size() == 0 && data_it_1 < data_it_end)
          {
            ++data_it_1;      // get next first DataPoint
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          if (data_it_1 == data_it_end && data_it_2 == data.end())     // if first iterator points to last element of "data" and second iterator points to end of "data"
          {            
            break;      // stop combining
          }

          while (data_it_2 < data.end() && data_it_2->size() == 0)      // as long as current second DataPoint is empty and second iterator does not point to end of "data"
          {
            ++data_it_2;      // get next second DataPoint
          }

          if (data_it_2 == data.end())      // if second iterator points to end of "data"
          {            
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          it_1 = data_it_1->begin();      // set first inner iterator to first element of first DataPoint
          it_2 = data_it_2->begin();      // set second inner iterator to first element of second DataPoint

          // check if DataPoints are not empty
          if (data_it_1->size() != 0 && data_it_2->size() != 0)
          {
            // check if DataPoints have the same charge state and mass shifts
            if (it_1->charge != it_2->charge || it_1->mass_shifts != it_2->mass_shifts)
            {              
              if (data_it_2 < data_it_end)     // if DataPpoints differ and second DataPoint is not second last element of "data"
              {
                temp++;
                ++data_it_2;      // get next second DataPoint
                if (temp > 50000)
                {                  
                  ++data_it_1;
                  temp = 0;
                }
              }

              else if (data_it_2 == data_it_end && data_it_1 < data.end() - 2)     // if DataPpoints differ and second DataPoint is second last element of "data" and first DataPoint is not third last element of "data"
              {
                ++data_it_1;      // get next first DataPoint
                data_it_2 = data_it_1 + 1;      // reset second iterator
              }

              else
              {
                ++data_it_1;      // get next first DataPoint
              }
            }

            else
            {              
              // perform combining
              (*data_it_1).insert(data_it_1->end(), data_it_2->begin(), data_it_2->end());      // append second DataPoint to first DataPoint
              (*data_it_2).clear();     // clear second Datapoint to keep iterators valid and to keep size of "data"

              if (data_it_2 < data_it_end)     // if second DataPoint is not second last element of "data"
              {
                ++data_it_2;      // get next second DataPoint
              }
              else
              {
                data_it_2 = data_it_1 + 1;      // reset second iterator
              }
            }
          }
          else
          {
            ++data_it_1;      // get next first DataPoint
          }
        }

        // erase empty DataPoints from "data"
        vector<vector<SILACPattern> > data_temp;

        for (vector<vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
        {
          if (data_it->size() != 0)
          {
            data_temp.push_back(*data_it);     // keep DataPoint if it is not empty
          }
        }

        data.swap(data_temp);     // data = data_temp
        data_temp.clear();      // clear "data_temp"
      }
    }

#if 0
    //--------------------------------------------------
    // remove isolated DataPoints from filter results (i.e. DataPoints with only few immediate neighbours)
    //--------------------------------------------------

    {
      Int immediate_neighbour_threshold = 5;      // maximum number of DataPoints within neighbourhood
      DoubleReal rt_neighbourhood = 10;     // size of neighbourhood for RT (+ and -)
      DoubleReal mz_neighbourhood = 0.02;     // size of neighbourhood for m/z (+ and -)

      vector<vector<DataPoint> > data_2;

      for (vector<vector<DataPoint> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
      {
        vector<DataPoint> single_layer;

        for (vector<DataPoint>::iterator it = data_it->begin(); it != data_it->end(); ++it)
        {
          Int immediate_neighbours = 0;

          for (vector<DataPoint>::iterator it_2 = data_it->begin(); it_2 != data_it->end(); ++it_2)
          {
            DoubleReal distance_mz = abs(it->mz - it_2->mz);
            DoubleReal distance_rt = abs(it->rt - it_2->rt);

            if (distance_rt < rt_neighbourhood && distance_mz < mz_neighbourhood)
            {
              ++immediate_neighbours;
            }
          }

          if (immediate_neighbours > immediate_neighbour_threshold)
          {
            single_layer.push_back(*it);
          }
        }

        data_2.push_back(single_layer);
      }

      data.swap(data_2);      // data = data_2
    }
#endif

    //--------------------------------------------------
    // store filter results from vector<vector<DataPoint> > data to .featureXML
    //--------------------------------------------------

    if (out_filters != "" && in_filters == "")     // check if option "out_filters" is specified and "in_filters" is not
    {
      ConsensusMap map;
      for (std::vector<std::vector<SILACPattern> >::const_iterator it = data.begin(); it != data.end(); ++it)
        generateFilterConsensusByPattern(map, *it);
      writeConsensus(out_filters, map);
    }

    //--------------------------------------------------
    // load filter results as vector<vector<DataPoint> > data from .featureXML
    //--------------------------------------------------

    if (in_filters != "")     // check if option "in_filters" is specified
    {
      ConsensusMap map;
      readConsensus(in_filters, map);
      readFilterConsensusByPattern(map);
    }
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


    //--------------------------------------------------
    // filter input data
    //--------------------------------------------------

    filterData(exp);


    //--------------------------------------------------
    // clustering
    //--------------------------------------------------

    vector<Cluster> clusters;

    clusterData();

/*    // store silhouettes for all subtrees as .csv
    Size silhouette_number = 1;

    ofstream silhouettesFile;
    silhouettesFile.open ("silhouettes.csv");

    for (std::vector<std::vector<Real> >::iterator asw_vector_it = silhouettes.begin(); asw_vector_it != silhouettes.end(); ++asw_vector_it)
    {
      silhouettesFile << silhouette_number;

      for (std::vector<Real>::iterator asw_it = asw_vector_it->begin(); asw_it != asw_vector_it->end(); ++asw_it)
      {
        silhouettesFile << "," << *asw_it;
      }

      silhouettesFile << endl;
      silhouette_number++;
    }

    silhouettesFile.close();

    // store all subtrees as .featureXML
    std::vector<String> colors;
    // 15 HTML colors
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

    for (std::vector<std::vector<SILACTreeNode> >::iterator subtree_it = subtrees.begin(); subtree_it != subtrees.end(); ++subtree_it)
    {
      std::set<DataPoint*> leafs;

      for (std::vector<SILACTreeNode>::iterator tree_it = subtree_it->begin(); tree_it != subtree_it->end(); ++tree_it)
      {
        leafs.insert(tree_it->data1);
        leafs.insert(tree_it->data2);
      }

      for (std::set<DataPoint*>::iterator leafs_it = leafs.begin(); leafs_it != leafs.end(); ++leafs_it)
      {
        //visualize the light variant
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

    //subtree_points.sortByPosition();
*/

#if 0
    //--------------------------------------------------
    // consensusXML output
    //--------------------------------------------------

    if (out != "")
    {
      Size id = 0;

      // iterate over all clusters
      for (vector<Cluster>::iterator cluster_it = clusters.begin(); cluster_it != clusters.end(); ++cluster_it)
      {
        DoubleReal rt = 0.0;
        DoubleReal mz = 0.0;
        DoubleReal total_intensity = 0.0;
        Size mass_shifts_size = (*(cluster_it->begin()))->mass_shifts.size();

        // create a vector with the maximum intensity of each isotope peak
        vector<DoubleReal> max_intensities(mass_shifts_size, 0.0);
        Int charge = (*(cluster_it->begin()))->charge;

        DoubleReal mass_shift_final = 0;

        // mass shifts as value for Quality in [Da] (i.e. 6008)
        for (vector<DoubleReal>::iterator shift_it = (*(cluster_it->begin()))->mass_shifts.begin(); shift_it != (*(cluster_it->begin()))->mass_shifts.end(); ++shift_it)
        {

          // mass shifts for doublets
          if ((*(cluster_it->begin()))->mass_shifts.size() == 2)
          {
            mass_shift_final = floor(*shift_it * charge);     // mass shift as value for Quality
          }

          // mass shifts for triplets, quadruplets, ...
          if ((*(cluster_it->begin()))->mass_shifts.size() > 2)
          {
            if (shift_it == (*(cluster_it->begin()))->mass_shifts.begin() + 1)
            {
              mass_shift_final = floor(*shift_it * charge) * 1000;      // mass shift as value for Quality
            }

            if ((shift_it > (*(cluster_it->begin()))->mass_shifts.begin() + 1) && (shift_it < (*(cluster_it->begin()))->mass_shifts.end() - 1))
            {
              mass_shift_final = (mass_shift_final + floor(*shift_it * charge)) * 1000;     // mass shift as value for Quality
            }

            if (shift_it == (*(cluster_it->begin()))->mass_shifts.end() - 1)
            {
              mass_shift_final += floor(*shift_it * charge);      // mass shift as value for Quality
            }
          }
        }

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

         // create consensus feature for each mass shift !=0
        ConsensusFeature consensus_feature;
        consensus_feature.setRT(rt);
        consensus_feature.setMZ(mz);
        consensus_feature.setIntensity(max_intensities[0]);			// set intensity of light peptiide to intensity of consensus
        consensus_feature.setCharge(charge);
        consensus_feature.setQuality(mass_shift_final);			// set mass shifts as value for Quality in [Da]

        // insert feature handle for each mass shift
        for (Size l = 0; l < mass_shifts_size; ++l)
        {
          // perform linear regression for each mass shift != 0

          Math::LinearRegression linear_reg;
          linear_reg.computeRegressionNoIntercept(0.95, intensities[0].begin(), intensities[0].end(), intensities[l].begin());

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
        ++id;
      }

      // set type of experiment
      all_pairs.setExperimentType("silac");

      // set mapList entries
      for (UInt i = 0; i <= massShifts[0].size(); i++)
      {
        all_pairs.getFileDescriptions()[i].filename = in;
        all_pairs.getFileDescriptions()[i].size = id;
      }      
    }
#endif

    //--------------------------------------------------------------
    // write output
    //--------------------------------------------------------------

    if (out != "")
    {
      ConsensusMap map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        generateClusterConsensusByCluster(map, **it);
      }
      writeConsensus(out, map);
    }

    if (out_clusters != "")
    {
      ConsensusMap map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        generateClusterConsensusByPattern(map, **it);
      }
      writeConsensus(out_clusters, map);
    }

    return EXECUTION_OK;
  }

  void clusterData();

private:
  void generateClusterConsensusByCluster(ConsensusMap &, const Clustering &) const;
  void generateClusterConsensusByPattern(ConsensusMap &, const Clustering &) const;
  void generateFilterConsensusByPattern(ConsensusMap &, const std::vector<SILACPattern> &) const;
  ConsensusFeature generateSingleConsensusByPattern(const SILACPattern &) const;
  void readFilterConsensusByPattern(ConsensusMap &);

  static const String &selectColor(UInt nr);
 
  void readConsensus(const String &filename, ConsensusMap &in) const
  {
    ConsensusXMLFile c_file;
    c_file.load(filename, in);
  }

  void writeConsensus(const String &filename, ConsensusMap &out) const
  {
    out.sortByPosition();
    out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    out.setExperimentType("silac");

    ConsensusXMLFile c_file;
    c_file.store(filename, out);
  }
};

void TOPPSILACAnalyzer::clusterData()
{
  typedef Clustering::PointCoordinate PointCoordinate;

  ProgressLogger progresslogger;
  progresslogger.setLogType(log_type_);

  progresslogger.startProgress(0, data.size(), "clustering data");

  UInt nr = 0;

  for (vector<vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it, ++nr)
  {
    const PointCoordinate max_delta = {{rt_threshold, mz_threshold}};
    Clustering *clustering = new Clustering(max_delta);

    for (vector<SILACPattern>::iterator it = data_it->begin(); it != data_it->end(); ++it)
    {
      const PointCoordinate key = {{it->rt, it->mz}};
      SILACPattern &p = *it;
      clustering->insertPoint(key, &p);
    }

    clustering->cluster();

    cluster_data.push_back(clustering);
  }

  progresslogger.endProgress();
}

void TOPPSILACAnalyzer::generateClusterConsensusByCluster(ConsensusMap &out, const Clustering &clustering) const
{
  for (Clustering::Grid::const_grid_iterator cell_it = clustering.grid.grid_begin(); cell_it != clustering.grid.grid_end(); ++cell_it)
  {
    std::ostringstream o;
    o << cell_it->first[0] << ':' << cell_it->first[1];
    std::string cell_id = o.str();

    for (Clustering::Grid::const_cell_iterator cluster_it = cell_it->second.begin(); cluster_it != cell_it->second.end(); ++cluster_it)
    {
      ConsensusFeature cluster;
      cluster.setMetaValue("Cell ID", cell_id);

      DoubleReal total_rt = 0;
      DoubleReal total_mz = 0;
      DoubleReal total_intensity = 0;
      std::vector<DoubleReal> total_intensities(10, 0);

      Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
      UInt pattern_id = 0;

      SILACPattern *pattern_max = &*pattern_it->second;

      for (; pattern_it != cluster_it->second.end();
             ++pattern_it, ++pattern_id)
      {
        SILACPattern *pattern = &*pattern_it->second;

        if (pattern_max->intensities[0][0] < pattern->intensities[0][0])
        {
          pattern_max = pattern;
        }

        for (UInt i = 0; i < pattern->mass_shifts.size(); ++i)
        {
          std::vector<DoubleReal> &current_intensities = pattern->intensities[i];
          // Find maximum intensity over all peaks
          std::vector<DoubleReal>::const_iterator max_intensity = max_element(current_intensities.begin(), current_intensities.end());
          if (max_intensity != current_intensities.end() && *max_intensity > 0)
            total_intensities[i] += *max_intensity;
        }

        total_intensity += pattern->intensities[0][0];
        // add monoisotopic RT position of the light peak to the total RT value, weighted by the intensity
        total_rt += pattern->intensities[0][0] * pattern->rt;
        total_mz += pattern->intensities[0][0] * pattern->mz;
      }

      // Average RT of monoisotopic peak
      total_rt /= total_intensity;
      cluster.setRT(total_rt);
      // Average MZ of monoisotopic peak
      total_mz /= total_intensity;
      cluster.setMZ(total_mz);
      // XXX: Max intensity or sum?
      cluster.setIntensity(total_intensities[0]);

      UInt charge = pattern_max->charge;
      cluster.setCharge(charge);

      {
        std::ostringstream mass_shifts_out;
        mass_shifts_out << std::fixed << std::setprecision(4);

        for (UInt i = 0; i < pattern_max->mass_shifts.size(); ++i)
        {
          // Product feature for each mass shift
          FeatureHandle point;
          point.setRT(total_rt);
          point.setMZ(total_mz + pattern_max->mass_shifts[i]);
          point.setIntensity(total_intensities[i]);
          point.setUniqueId(i);

          cluster.insert(point);

          // Product mass shifts string
          mass_shifts_out << pattern_max->mass_shifts[i] * charge << ';';
        }

        // Remove the last delimiter
        std::string mass_shifts_outs = mass_shifts_out.str(); mass_shifts_outs.erase(mass_shifts_outs.end() - 1);
        cluster.setQuality(std::floor(pattern_max->mass_shifts.at(1) * charge));
        cluster.setMetaValue("Mass shifts [Da]", mass_shifts_outs);
      }

      out.push_back(cluster);
    }
  }
}

void TOPPSILACAnalyzer::generateClusterConsensusByPattern(ConsensusMap &out, const Clustering &clustering) const
{
  UInt cluster_id = 0;

  for (Clustering::Grid::const_grid_iterator cell_it = clustering.grid.grid_begin(); cell_it != clustering.grid.grid_end(); ++cell_it)
  {
    std::ostringstream o;
    o << cell_it->first[0] << ':' << cell_it->first[1];
    std::string cell_id = o.str();

    for (Clustering::Grid::const_cell_iterator cluster_it = cell_it->second.begin(); cluster_it != cell_it->second.end(); ++cluster_it, ++cluster_id)
    {
      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin(); pattern_it != cluster_it->second.end(); ++pattern_it)
      {
        ConsensusFeature consensus = generateSingleConsensusByPattern(*pattern_it->second);

        consensus.setMetaValue("color", selectColor(cluster_id));
        consensus.setMetaValue("Cluster ID", cluster_id);
        consensus.setMetaValue("Cell ID", cell_id);

        out.push_back(consensus);
      }
    }
  }
}

void TOPPSILACAnalyzer::generateFilterConsensusByPattern(ConsensusMap &out, const std::vector<SILACPattern> &pattern) const
{
  for (std::vector<SILACPattern>::const_iterator pattern_it = pattern.begin(); pattern_it != pattern.end(); ++pattern_it)
  {
    out.push_back(generateSingleConsensusByPattern(*pattern_it));
  }
}

ConsensusFeature TOPPSILACAnalyzer::generateSingleConsensusByPattern(const SILACPattern &pattern) const
{
  // XXX: get from experiment
  Int charge = pattern.charge;

  ConsensusFeature consensus;
  consensus.setRT(pattern.rt);
  consensus.setMZ(pattern.mz);
  consensus.setIntensity(pattern.intensities[0][0]);
  consensus.setCharge(charge);

  consensus.setMetaValue("Peaks per peptide", pattern.isotopes_per_peptide);

  // Output mass shifts
  {
    std::ostringstream out;
    out << std::fixed << std::setprecision(4);
    for (vector<DoubleReal>::const_iterator shift_it = pattern.mass_shifts.begin() + 1; shift_it != pattern.mass_shifts.end(); ++shift_it)
    {
      out << *shift_it * charge << ';';
    }
    // Remove the last delimiter
    std::string outs = out.str(); outs.erase(outs.end() - 1);
    consensus.setQuality(std::floor(pattern.mass_shifts.at(1) * charge));
    consensus.setMetaValue("Mass shifts [Da]", outs);
  }

  // Output all intensities per peptide as list
  {
    std::ostringstream out;
    for (vector<vector<DoubleReal> >::const_iterator inten_it = pattern.intensities.begin(); inten_it != pattern.intensities.end(); ++inten_it)
    {
      std::ostringstream out2;
      out2 << std::fixed << std::setprecision(4);
      for (vector<DoubleReal>::const_iterator inten2_it = inten_it->begin(); inten2_it != inten_it->end(); ++inten2_it)
      {
        out2 << *inten2_it << ',';
      }
      // Remove the last delimiter
      std::string out2s = out2.str(); out2s.erase(out2s.end() - 1);
      out << out2s << ';';
    }
    // Remove the last delimiter
    std::string outs = out.str(); outs.erase(outs.end() - 1);
    consensus.setMetaValue("Intensities", outs);
  }

  UInt point_id = 0;
  for (std::vector<DataPoint>::const_iterator point_it = pattern.points.begin();
       point_it != pattern.points.end();
       ++point_it, ++point_id)
  {
    FeatureHandle point;
    point.setRT(point_it->rt);
    point.setMZ(point_it->mz);
    point.setUniqueId(point_id);

    consensus.insert(point);
  }

  return consensus;
}

void TOPPSILACAnalyzer::readFilterConsensusByPattern(ConsensusMap &in)
{
  std::map<std::pair<Int, Int>, std::vector<SILACPattern> > layers;

  for (ConsensusMap::const_iterator pattern_it = in.begin(); pattern_it != in.end(); ++pattern_it)
  {
    SILACPattern pattern;
    pattern.rt = pattern_it->getRT();
    pattern.mz = pattern_it->getMZ();
    pattern.charge = pattern_it->getCharge();
    pattern.quality = pattern_it->getQuality();

    pattern.isotopes_per_peptide = pattern_it->getMetaValue("Peaks per peptide");

    StringList text = StringList::create(pattern_it->getMetaValue("Mass shifts [Da]"), ';');
    pattern.mass_shifts.push_back(0);
    for (StringList::const_iterator text_it = text.begin(); text_it != text.end(); ++text_it)
    {
      pattern.mass_shifts.push_back(text_it->toDouble() / pattern.charge);
    }

    text = StringList::create(pattern_it->getMetaValue("Intensities"), ';');
    for (StringList::const_iterator text_it = text.begin(); text_it != text.end(); ++text_it)
    {
      StringList text2 = StringList::create(*text_it, ',');
      vector<DoubleReal> inten;
      for (StringList::const_iterator text2_it = text2.begin(); text2_it != text2.end(); ++text2_it)
      {
        inten.push_back(text2_it->toDouble());
      }
      pattern.intensities.push_back(inten);
    }

    for (ConsensusFeature::const_iterator point_it = pattern_it->begin(); point_it != pattern_it->end(); ++point_it)
    {
      DataPoint point;
      point.rt = point_it->getRT();
      point.mz = point_it->getMZ();

      pattern.points.push_back(point);
    }

    layers[std::make_pair(Int(pattern.mass_shifts.at(1)), pattern.charge)].push_back(pattern);
  }

  for (std::map<std::pair<Int, Int>, std::vector<SILACPattern> >::iterator it = layers.begin(); it != layers.end(); ++it)
  {
    data.push_back(it->second);
  }
}

const String &TOPPSILACAnalyzer::selectColor(UInt nr)
{
  // 15 HTML colors
  const static String colors[] = {
    "#00FFFF", "#000000", "#0000FF", "#FF00FF", "#008000",
    "#808080", "#00FF00", "#800000", "#000080", "#808000",
    "#800080", "#FF0000", "#C0C0C0", "#008080", "#FFFF00",
  };
  const Int colors_len = 15;

  return colors[nr % colors_len];
}

//@endcond

int main(int argc, const char** argv )
{
  TOPPSILACAnalyzer tool;
  return tool.main(argc, argv);
}
