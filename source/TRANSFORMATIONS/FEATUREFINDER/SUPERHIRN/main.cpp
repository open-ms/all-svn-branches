
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <vector>
#include <sstream>
#include <stdlib.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/IsotopicDist.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Deisotoper.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Process_Data.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms2_info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2_feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/featureLCprofile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS.h>;
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>;
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1_feature_merger.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PEAK_DETEC_mzXML_reader.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PeakDetectController.h>

using namespace std;

void initParams() {
  Process_Data::CENTROID_DATA_MODUS = 1; // data is centroided already 
  
   
  FT_PEAK_DETEC_mzXML_reader::PEAK_EXTRACTION_SCAN_LEVELS.push_back(1);
  FT_PEAK_DETEC_mzXML_reader::FRAGMENT_MASS_SCAN_LEVELS.push_back(2);
  FT_PEAK_DETEC_mzXML_reader::MS1_base_inter_scan_distance = 0.1;  
  Process_Data::MS1_TR_RESOLUTION = 0.01;
  
  float thresh = 1000;
  LCMSCData::intensity_min_threshold = thresh;
  Process_Data::INTENSITY_THRESHOLD = thresh;
  Process_Data::max_inter_scan_retention_time_distance = 0.1;
  Process_Data::min_nb_cluster_members = 4;
  
  //def->search_tag("MS1 feature CHRG range min", &INT);
  Deisotoper::sfMinCharge = 1;
  //def->search_tag("MS1 feature CHRG range max", &INT);
  Deisotoper::sfMaxCharge = 5;
  //def->search_tag("Detectable isotope factor",&DB);
  IsotopicDist::sfDetectableIsoFact = 0.05;
  // def->search_tag("IntensityCV",&DB);
  IsotopicDist::sfIntensityCV = 0.9;
  
  MS1_feature_merger::MS1_FEATURE_CLUSTERING = 1;
  MS1_feature_merger::MS1_PEAK_AREA_TR_RESOLUTION = 0.01;
  MS1_feature_merger::INITIAL_TR_TOLERANCE = 5.0;
  MS1_feature_merger::MS1_FEATURE_MERGING_TR_TOLERANCE = 1.0;
  MS1_feature_merger::PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION = 25;
  MS1_feature_merger::PPM_TOLERANCE_FOR_MZ_CLUSTERING = 10;
  
  
  // ----------------------------------------------------------------------
  // aus void FT_PEAK_DETECT_initializer::init_all(){
  // ----------------------------------------------------------------------
  
  //def->search_tag("Centroid window width",&INT);
  CentroidPeak::sfCentroidWindowWidth = 5;
  //def->search_tag("Absolute isotope mass precision",&DB);
  CentroidData::sfMassTolDa = 0.01;
  //def->search_tag("Relative isotope mass precision",&DB);
  CentroidData::sfMassTolPpm = 10;
  //def->search_tag("Minimal peak height",&DB);
  CentroidData::sfMinIntensity = 0.0;
  //def->search_tag("Min. Centroid MS Signal Intensity",&DB);
  CentroidData::sfIntensityFloor = 50; //1.0; // in config its 50, but in CentroidData it's 1;
  
  //def->search_tag("Report mono peaks",&INT);
  FT_PEAK_DETEC_mzXML_reader::sfReportMonoPeaks = 0;
  //def->search_tag("Report scan number",&INT);
  FT_PEAK_DETEC_mzXML_reader::sfReportScanNumber = 0;
  
  
  
  // ----------------------------------------------------------------------
  // aus initializer
  // ----------------------------------------------------------------------
  
  //def->search_tag("start elution window", &TMP);
  FT_PEAK_DETEC_mzXML_reader::TR_MIN=0;
  
  //def->search_tag("end elution window", &TMP);
  FT_PEAK_DETEC_mzXML_reader::TR_MAX=180;
  
  // feature parameters:
  //def->search_tag("MS1 retention time tolerance", &TMP);
  feature::TR_TOL = 0.5;
//  def->search_tag("MS1 m/z tolerance", &TMP);
  feature::PPM_MZ_TOL = 0;
  consensIsotopePattern::FT_MZ_TOLERANCE = feature::PPM_MZ_TOL;
  
  // MS2_M2_matcher parameters:
  //def->search_tag("MS2 mass matching modus", &TMP_B);
  ms2_info::THEO_MATCH_MODUS = 1;
  
  
  /*
   
  //def->search_tag("Peptide Prophet Threshold", &TMP);
  MS2_MS1_matcher::PEP_PROB_CUT_OFF = 0.9;
  //def->search_tag("MS2 SCAN tolerance", &TMP_I);
  MS2_MS1_matcher::SCAN_TOL = 200;
  //def->search_tag("INCLUSIONS LIST MS2 SCAN tolerance", &TMP_I);
  MS2_MS1_matcher::POST_SCAN_TOL = 100000;
  
   */
   
   
  // MS2 matching PPM parameters:
  //def->search_tag("MS2 PPM m/z tolerance", &TMP);
  ms2_info::MS2_MZ_PPM_TOLERANCE = 30;
  
  // MS2 retention time tolerance:
  //def->search_tag("MS2 retention time tolerance", &TMP);
  //if( TMP > 0 ){
  //  ms2_info::MS2_TR_TOL = TMP; 
  //}
  //else{
    ms2_info::MS2_TR_TOL = feature::TR_TOL;
  //}
  
  // peptide_DELTA_group:
  //def->search_tag("MS1 retention time tolerance", &TMP);
  // ------ peptide_DELTA_group::TR_TOL = 0.5;
  //def->search_tag("Delta pair TR tolerance", &TMP);  
  //peptide_DELTA_group::MOD_TR_TOL = TMP;
  // if smaller than -1, then use the general TR tolerance:
  //if( peptide_DELTA_group::MOD_TR_TOL < peptide_DELTA_group::TR_TOL ){
    // ------ peptide_DELTA_group::MOD_TR_TOL = peptide_DELTA_group::TR_TOL;
  //}
  //def->search_tag("Delta pair TR tolerance", &TMP);  
  
  
  // protein group:
  //def->search_tag("Peptide Proteotype Mode", &TMP_B);
  // ------ protein_group::PEPTIDE_PROTEOTYPE_MODE = 1;
  
  
//  def->search_tag("Peptide Prophet Threshold", &TMP);
  double d = 0.9;
  // ------ peptide_DELTA_group::PEPTIDE_PROBABILITY_THRESHOLD = d;
  feature::PEPTIDE_PROBABILITY_THRESHOLD = d;
  // ------ interact_parser::PEPTIDE_PROBABILITY_THRESHOLD = d;
  LC_MS::PEP_PROPHET_THERSHOLD = d;
  
  // FLO: TODO - what is this about?
  // slide window over TR
  //map<int, vector<double> > TMP2;
  //def->search_tag("Delta M/z list", &TMP2);
  //peptide_DELTA_group::LC_MS_Modification_Masses = TMP2; 
  //MS1_feature_ratiolizer::LC_MS_Modification_Masses = TMP2;
  //peptide_ratio_analyzer::LC_MS_Modification_Masses = TMP2;
  
  
  
  
  /* ------------------------- classes I do not have...
  
  
  // DELTA Pair MATCHING
  //def->search_tag("MS2 Delta clustering filter", &TMP_B);
  DELTA_grouper::MS_2_SELECTION = 0;
  
  
  // STATIC MODIFICATIONS:
  //def->search_tag("static glyco modifictation", &TMP_B);
  interact_parser::STATIC_GLYCOSYLATION_MOD = 0;
  //def->search_tag("static C-term. modification", &TMP);
  interact_parser::STATIC_C_TERM_MODIFICATION = 0;
  
  // FLO: TODO What kind of paramters are this...
  //map<double,double> M_TMP;
  //def->search_tag("INTERACT AA MOD transform table", &M_TMP);
  //interact_parser::MOD_MASS_TRANSFORM_TABLE = M_TMP;
  
  //def->search_tag("MS2 mass type check", &TMP_B);
  interact_parser::MS2_MASS_COMPARE = 0;
  
  //////////////////
  // log modus or not:
  bool LOG_MODE_INTENSITY = false;
  consens_profile_builder::LOG_INTENSITY_TRANSFORMATION = LOG_MODE_INTENSITY;
  profile_scorer::LOG_INTENSITY_TRANSFORMATION = LOG_MODE_INTENSITY;
  
  
  /////////////////////////
  // LC-MS correlation:
  
  //def->search_tag("pairwise correlation analysis", &TMP_B);
  LC_MS_correlation::VIEW = 0;
  
  // Tr tolerance used to compare
  // spec A to B after smoothing and
  // also used in the merging
  //def->search_tag("MS1 retention time tolerance", &TMP);
  LC_MS_correlation::tr_tol = 0.5;
  
  // M/z tolerance used to compare
  // spec A to B after smoothing and
  // also used in the merging
  //def->search_tag("MS1 m/z tolerance", &TMP);
  LC_MS_correlation::mz_tol = 10;
  
  // numbers of bins to divide the peak areas
  // used afterwards to compare peak areas of 2 peaks
  //def->search_tag("intensity bin size", &TMP_I);
  LC_MS_correlation::intensity_bin_size = 2000;
  
  // minimal alignment error = TR tolerance / 2:
//  def->search_tag("MS1 retention time tolerance", &TMP);
//  TMP /= 2.0;
  LC_MS_correlation::min_align_error = 0.25;
  
  // maximal alignment error:
  //def->search_tag("maximal smoothing error", &TMP);
  LC_MS_correlation::max_align_error = 3.0;
  
  // tolerance of bins, how far 2 bins can be apart and still
  // be a match
  //def->search_tag("intensity bin tolerance", &TMP);
  LC_MS_correlation::intensity_bin_tolerance = 2;
  
  // represents the worst score possible, this one will be used to 
  // normaize the observed scores between 0(bad) and 1(good) [ 0 ... 1]
  //def->search_tag("minimal LC/MS score", &TMP);
  LC_MS_correlation::min_LC_MS_score = 0.1;
  
  
  
  //////////////////////////////////////////////////////////
  // LC_MS alignmnt parameters:
  
  // retention time tolerance:
  //def->search_tag("retention time window", &TMP);
  LC_MS_aligner::Tr_window = 5.0;
  
  // mass tolerance:
  //def->search_tag("mass / charge window", &TMP);
  LC_MS_aligner::Mz_window = 20;  
  
  // percentage of outside delta error points
  //def->search_tag("perc. outside error delta points", &TMP);
  LC_MS_aligner::ERROR_DELTA_POINTS = 0.75;
  
  // Tr tolerance used to compare
  // spec A to B after smoothing and
  // also used in the merging
  //def->search_tag("MS1 retention time tolerance", &TMP);
  LC_MS_aligner::Tr_tolerance = 0.5;
  
  // M/z tolerance used to compare
  // spec A to B after smoothing and
  // also used in the merging
  //def->search_tag("MS1 m/z tolerance", &TMP);
  LC_MS_aligner::Mz_tolerance = 10;
  
  // defines if also identifications should be checked
  // in finding common lc/ms peaks between runs:
  //def->search_tag("MS2 info alignment weight", &TMP_B);
  LC_MS_aligner::ID_ALIGNMENT_WEIGHT = 1;
  
  // plotting:
  //def->search_tag("pairwise alignment plotting", &TMP_B);
  LC_MS_aligner::print_out = 0;
  
  
  /////////////////////////////////////////////////
  // these are the important parameters for the
  // for the lowess fitter druing the LC/MS alignment:
  
  // minimal alignment error = TR tolerance / 2:
  //def->search_tag("MS1 retention time tolerance", &TMP);
  regressor::min_error = 0.5;
  
  // number of boostrape cycles:
  //def->search_tag("maximal smoothing error", &TMP);
  regressor::max_error = 3.0;
  
  // max stripes for the plot smoother:
  //def->search_tag("max. nb. stripes", &TMP_I);
  plot_smoother::max_nb_stripes = 1;  
  
  // used to copmute the alignment error, use a tr window to
  // calculate the standard deviations of raw data to predicted
  // delta shift
  //def->search_tag("smoothing error TR window", &TMP);
  regressor::tr_error_smooth = 1.0;
  
  
  
  */
  
  
  /////////////////////////////////////////////////
  // parameter if the LC elution profile will be stored 
  // in the XML:
  //def->search_tag("Create monoisotopic LC profile", &TMP_B);
  // ------- LCMSDataImporter::CREATE_FEATURE_ELUTION_PROFILES = 1;
  FT_PeakDetectController::CREATE_FEATURE_ELUTION_PROFILES = 1;
  // ------- LC_MS_XML_writer::STORE_FEATURE_ELUTION_PROFILES = 1;
  
  /////////////////////////////////////////////////
  // what and how XML data is stored in the mastermap:
  // ms2 information of a feature:
  // only the best ms2 info / feature stored:
  
  /*def->search_tag("store only best MS2 per feature", &TMP_B);
  LC_MS_XML_writer::STORE_BEST_MS2_SCAN_PER_FEATURE = TMP_B;
  // only the best ms2 info / aligned feature stored:
  def->search_tag("store only best MS2 per ALIGNED feature", &TMP_B);
  LC_MS_XML_writer::STORE_BEST_MS2_SCAN_PER_FEATURE = TMP_B;
  // how many alternative protein names to store:
  def->search_tag("nb. max. alternative protein names", &TMP_I);
  LC_MS_XML_writer::MAXIMAL_NB_ALTERNATIVE_PROTEIN_NAMES = TMP_I;
  // if to store ms2 traces
  def->search_tag("MS2 fragment mass tracing", &TMP_B);
  LC_MS_XML_writer::STORE_MS2_FRAGMENT_TRACE_DATA = TMP_B;
   */
  
  
  
  /////////////////////////////////////////////////
  // Parameters for the peak merging:
  //def->search_tag("Activation of MS1 feature merging post processing", &TMP_B);
  MS1_feature_merger::MS1_FEATURE_CLUSTERING = 1;
  //def->search_tag("MS1 LC retention time resolution", &TMP);
  MS1_feature_merger::MS1_PEAK_AREA_TR_RESOLUTION = 0.01;
  
  //def->search_tag("Initial Apex Tr tolerance", &TMP);
  MS1_feature_merger::INITIAL_TR_TOLERANCE = 5.0;
  //def->search_tag("MS1 feature Tr merging tolerance", &TMP);
  MS1_feature_merger::MS1_FEATURE_MERGING_TR_TOLERANCE = 1.0;
  //def->search_tag("Percentage of intensity variation between LC border peaks", &TMP);
  MS1_feature_merger::PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION = 25;
  //def->search_tag("PPM value for the m/z clustering of merging candidates", &TMP);
  MS1_feature_merger::PPM_TOLERANCE_FOR_MZ_CLUSTERING = 10;
  
  
  /////////////////////////////////////////////////
  // what information is extracted from the LC/MS or mastermap:
  // TR min:
  
//  def->search_tag("start elution window", &TMP);
  LC_MS_XML_reader::TR_MIN = 0;
  // TR max:
//  def->search_tag("end elution window", &TMP);
  LC_MS_XML_reader::TR_MAX = 180;
  // mz min.
  //def->search_tag("MS1 feature mz range min", &TMP);
  LC_MS_XML_reader::FEATURE_MZ_MIN = 0;
  // mz max.
  //def->search_tag("MS1 feature mz range max", &TMP );
  LC_MS_XML_reader::FEATURE_MZ_MAX = 2000;
  // signal to noise min.
  //def->search_tag("MS1 feature signal to noise threshold", &TMP );
  //LC_MS_XML_reader::SIGNAL_TO_NOISE_THERSHOLD = TMP;
  // intensity min.
//  def->search_tag("MS1 feature intensity cutoff", &TMP );
//  LC_MS_XML_reader::PEAK_INTENSITY_THRESHOLD = TMP; // ich glaube 10000
  // charge state min.
  //def->search_tag("MS1 feature CHRG range min", &TMP_I );
  LC_MS_XML_reader::FEATURE_CHRG_MIN = 1;
  // charge state max.
  //def->search_tag("MS1 feature CHRG range max", &TMP_I );
  LC_MS_XML_reader::FEATURE_CHRG_MAX = 5;
  //  Create monoisotopic LC profile:	to create and store the original profile of the detected
  //					monosiotopic pecursors in the XML (!!! increases the
  //					XML file size!!! (on[1]/off[0])
  //def->search_tag("Create monoisotopic LC profile", &TMP_B);
  //LC_MS_XML_reader::EXTRACT_MONO_ISOTOPE_PROFILE = TMP_B;
  
  
  
  /////////////////////////////////////////////////
  // what and how data is stored during superhirn processing:
  // ms2 information of a feature:
  // only the best ms2 info / feature stored:
  //def->search_tag("progress all low probability ms2 info in MS1 feature", &TMP_B);
  feature::STORE_ALL_LOW_PROBABILITY_MS2_SCANS = 0;
  
  
  // XML Data format to use during SuperHirn processing:
  //LC_MS_XML_reader::DATA_STORAGE_XML_FORMAT_TYPE = def->search_tag("SuperHirn Data Storage XML Output Format");
}

typedef map<double, RawData*> Map;
typedef vector<Map> Vec;

int main (int argc, char * const argv[]) {
  
  Map* datamap = new Map();  
  Vec* datavec = new Vec();
  
  ifstream file;
  string fileName ="/Users/zellerf/nulogs/ffsh-prec.txt";
  file.open(fileName.c_str());
  string line;
  
  while(getline(file, line))
  {
    if (line.substr(0, 8) == "Spectrum") {
      string scan = line.substr(8);
    }
    
    stringstream linestream (stringstream::in | stringstream::out);
    stringstream linestream2 (stringstream::in | stringstream::out);
    stringstream ss (stringstream::in | stringstream::out);
    double rt = 0;
    double val = 0;
    string valWithComma;
    vector<double>* vmzvals = new vector<double>();
    vector<double>* vintvals = new vector<double>();
    
    // RT value
    getline(file, line);
    if (line.substr(0, 2) != "RT") {
      cout << "IS NOT RT!\n";
      return 0;
    }
    ss << line.substr(5);
    ss >> rt;
    
    // mz values
    getline(file, line); 
    linestream << line;
    while (linestream >> valWithComma) {
      string tmp = valWithComma.substr(0, valWithComma.length()-1);
      
      // did not manage to stream it, shame on me
      val = atof(tmp.c_str());
      vmzvals->push_back(val);
    }
    
    // int values
    getline(file, line); 
    linestream2 << line;
    while (linestream2 >> valWithComma) {
      string tmp = valWithComma.substr(0, valWithComma.length()-1);
      
      // did not manage to stream it, shame on me
      val = atof(tmp.c_str());
      vintvals->push_back(val);
      
      if (val < 50.0) { //  < CentroidData::sfIntensityFloor
        cout << val << " < CentroidData::sfIntensityFloor\n";
      }
    }
    
    RawData* data = new RawData(*vmzvals, *vintvals);
    int t = vmzvals->size();
    
    
    vector<double> masses,intens;
    data->get(masses,intens);
    int t2 = masses.size();
    
    if (vmzvals->size() != vintvals->size()) {
      cout << vmzvals->size() << " != " << vintvals->size() << "!";
      return 0;
    }
    
    // Correcting the data (but it did not mather)
    if (rt == 1.0) {
      rt = -1.0;
    }
    
    Map m;
    m[rt/60.0] = data;
    
    datavec->push_back(m);
  }
  
  file.close();
  
  initParams();
  
  IsotopicDist::init(); // wichtig. mit test. with vim test
  
  FT_PeakDetectController controller;
  controller.start_scan_parsing_of_mzXML_file(*datavec);
}



