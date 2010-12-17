/*
 *  LC_MS_XML_reader.cpp
 *  CleanUpSH3
 *
 *  Created by zellerf on 12/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>

// These values are overwritten by config
double LC_MS_XML_reader::TR_MIN = 0;
double LC_MS_XML_reader::TR_MAX = 0; // 180
double LC_MS_XML_reader::FEATURE_MZ_MIN = 0; // 200
double LC_MS_XML_reader::FEATURE_MZ_MAX = 0; //1800;
int LC_MS_XML_reader::FEATURE_CHRG_MIN = 0; //1;
int LC_MS_XML_reader::FEATURE_CHRG_MAX = 0;

