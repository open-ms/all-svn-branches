/*
 *  simple_math2.h
 *  CleanUpSH2
 *
 *  Created by zellerf on 11/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _USE_simple_math2
#define _USE_simple_math2

class simple_math2 {
  
private:
  
  bool HIGH_CHECK;
  bool LOW_CHECK;

public:
  
  simple_math2();
  
  static const double T_TEST_001[12];
  static const double T_TEST_002[12];
  static const double T_TEST_01[12];
  static const double T_TEST_02[12];
  static const double T_TEST_05[12];
  static std::string ALPHA_VALUE;
  
  // this structure provides the function to compare
  // in the sorting algorithm:
  struct VECTOR_OPERATOR{
    // provide the compare function for sort:
    bool operator()(const std::pair<double,  void*> A,const std::pair<double,  void*> B) const{
      // check if they have same mass
      if(A.first == B.first){
        return false;
      }
      else{
        return A.first > B.first;
      }
    }
  };
  
  // if they fall into the m/z tolerance window
  static bool compareMassValuesAtPPMLevel( double, double , double );
  // get the masse error at the PPM value     
  static double getMassErrorAtPPMLevel( double , double );
  
  void ITERATIVE_OUTLIER_DETECTION_BY_DIXON(std::vector<double>* IN);
  
  void ITERATIVE_OUTLIER_DETECTION_BY_DIXON(std::vector< std::pair<double, void*> >* IN);
  
  void OUTLIER_DETECTION_BY_DIXON(std::vector<double>* IN);
  
  void OUTLIER_DETECTION_BY_DIXON(std::vector< std::pair<double, double> >* IN);
  
  void OUTLIER_DETECTION_BY_DIXON(std::vector< std::pair<double, void*> >* IN);
  
  bool check_T_TEST( double IN , int SAMPLE_NB);
    
};


#endif
