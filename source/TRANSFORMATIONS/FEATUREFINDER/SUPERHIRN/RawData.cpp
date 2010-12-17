/*
 *  RawData.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *
 */


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

#include <iostream>
#include <iomanip>
using namespace std;

// Constructor & destructor
///////////////////////////////////////////////////////////////////////////////


// More general constructor
RawData::RawData(
                 vector<double>& pMassValues, // mass sample values 
                 vector<double>& pIntensValues // intensity sample values 
                 )
{
  fProfileMasses = pMassValues;
  fProfileIntens = pIntensValues;
  LOW_INTENSITY_MS_SIGNAL_THRESHOLD = 1.0; 
}

// Destructor
RawData::~RawData()
{
}

// Operators

// Writes data to out stream using the << operator
ostream& operator<<(
                    ostream& pOut, // output stream 
                    RawData& pRawData) // 
{
  vector<double> m,h;
  vector<double>::iterator mi,hi;
  
  pRawData.get(m,h);
  for (mi=m.begin(),hi=h.begin();mi!=m.end();++mi,++hi) {
    pOut << fixed << setprecision(4) << *mi << " " << fixed << setprecision(2) << *hi << endl;
	}
  
  return pOut;
}


// Public methods

// Retrieve raw data as mass and intensity vectors
void RawData::get(
                  vector<double> &pProfileMasses, // Mass sample values in profile mode
                  vector<double> &pProfileIntens  // Intensity sample values in profile mode
                  ){
  
  pProfileMasses = fProfileMasses;
  pProfileIntens = fProfileIntens;

}

// Set raw data as mass and intensity vectors
void RawData::set(	vector<double> &pProfileMasses, // Mass sample values in profile mode
                   vector<double> &pProfileIntens  // Intensity sample values in profile mode
                   ){
  
	fProfileMasses = pProfileMasses;
	fProfileIntens = pProfileIntens;

}

