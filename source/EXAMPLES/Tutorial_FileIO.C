#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  
  MzXMLFile mzxml;
  MzMLFile mzml; 
  MzDataFile mzdata;
  MascotGenericFile mascot;

  // temporary data storage
  MSExperiment<Peak1D> map;

  // convert MzXML to MzML
  mzxml.load("C:/OpenMS/SWP2012/source/EXAMPLES/data/Tutorial_FileIO.mzXML",map);
  mzml.store("C:/OpenMS/SWP2012/source/EXAMPLES/output/Tutorial_FileIO.mzML",map);
  /*
  // convert MzXML to MzData
  mzxml.load("C:/OpenMS/SWP2012/source/EXAMPLES/data/Tutorial_FileIO.mzXML",map);
  mzdata.store("C:/OpenMS/SWP2012/source/EXAMPLES/output/Tutorial_FileIO.mzData",map);

  // convert MzXML to Mascot Generic 
  mzxml.load("C:/OpenMS/SWP2012/source/EXAMPLES/data/Tutorial_FileIO.mzXML",map);
  mascot.store("C:/OpenMS/SWP2012/source/EXAMPLES/output/Tutorial_FileIO.MGF",map);
  */
  return 0;
} //end of main
