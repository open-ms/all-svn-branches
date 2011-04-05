// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DTAFILE_H
#define OPENMS_FORMAT_DTAFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
 	/**
 		@brief File adapter for DTA files.
 		
  	The first line contains the singly protonated peptide mass (MH+) and the peptide charge state separated by a space.
  	Subsequent lines contain space separated pairs of fragment ion m/z and intensity values.
  	
  	From precusor mass and charge state the mass-charge-ratio is calculated and stored in the spectrum as precursor mass.
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI DTAFile
  {
    public:
      /// Default constructor
      DTAFile();
			/// Destructor
      virtual ~DTAFile();
      
      /**
 				@brief Loads a DTA file to a spectrum.
 				
 				The content of the file is stored in @p spectrum.
 				@p spectrum has to be a MSSpectrum or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      template <typename SpectrumType>
      void load(const String& filename, SpectrumType& spectrum)
      {
				std::ifstream is(filename.c_str());
				if (!is)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
		
				//  Delete old spectrum
				spectrum.clear(true);
				
				//temporary variables
				String line;
				std::vector<String> strings(2);
				typename SpectrumType::PeakType p;
				char delimiter;
			
				// line number counter
				Size line_number = 1;
							
				//read first line and store precursor m/z and charge
				getline(is,line,'\n');
				line.trim();
		
				//test which delimiter is used in the line
				if (line.has('\t'))
				{
					delimiter = '\t';
				}
				else
				{
					delimiter = ' ';
				}
				
				line.split(delimiter,strings);
				if (strings.size()!=2)
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line (" + String(line_number) + "): \"")+line+"\" (got  " + String(strings.size()) + ", expected 2 entries)" ,filename);
				}
				Precursor precursor;
        DoubleReal mh_mass;
        Int charge;
				try
				{
          // by convention the first line holds: singly protonated peptide mass, charge state
				  mh_mass = strings[0].toDouble();
					charge = strings[1].toInt();
				}
				catch(...)
				{
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line (" + String(line_number) + "): \"")+line+"\": not a float number." ,filename);
				}
        if (charge != 0)
				{
					precursor.setMZ( (mh_mass - Constants::PROTON_MASS_U) / charge + Constants::PROTON_MASS_U);
				}
				else
				{
					precursor.setMZ( mh_mass );
				}
				precursor.setCharge(charge);
				spectrum.getPrecursors().push_back(precursor);
				
		    while (getline(is,line,'\n'))
		    {
					++line_number;
		    	line.trim();
					if ( line.empty() ) continue;
					
					//test which delimiter is used in the line
					if (line.has('\t'))
					{
						delimiter = '\t';
					}
					else
					{
						delimiter = ' ';
					}
		
					line.split(delimiter,strings);
					if (strings.size()!=2)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line (" + String(line_number) + "): \"")+line+"\" (got  " + String(strings.size()) + ", expected 2 entries)" ,filename);
					}				
					try 
					{
						//fill peak
						p.setPosition((typename SpectrumType::PeakType::PositionType)strings[0].toDouble());
						p.setIntensity((typename SpectrumType::PeakType::IntensityType)strings[1].toDouble());
					} 
					catch ( Exception::BaseException & /*e*/ )
					{
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line (" + String(line_number) + "): \"")+line+"\": not a float number." ,filename);
					}
					spectrum.push_back(p);
				}
			
				spectrum.setName(File::basename(filename));	
				is.close();  	
      }

      /**
      	@brief Stores a spectrum in a DTA file.
      	
      	The content of @p spectrum is stored in a file.
      	@p spectrum has to be a MSSpectrum or have the same interface.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
      */
      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum) const
      {
				std::ofstream os(filename.c_str());
				if (!os)
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
				os.precision(writtenDigits<DoubleReal>());
				
				//write precursor information
				Precursor precursor;
				if (spectrum.getPrecursors().size()>0)
				{
					precursor = spectrum.getPrecursors()[0];
				}
				if (spectrum.getPrecursors().size()>1)
				{
					std::cerr << "Warning: The spectrum written to the DTA file '" << filename << "' has more than one precursor. The first precursor is used!" << "\n";
				}
				//unknown charge
				if (precursor.getCharge()==0)
				{
					os << precursor.getMZ();
				}
				//known charge
				else
				{
					os << ((precursor.getMZ() - 1.0) * precursor.getCharge() +1.0);
				}
				//charge
				os << " " << precursor.getCharge() << "\n";

				// Iterate over all peaks of the spectrum and
				// write one line for each peak of the spectrum.
				typename SpectrumType::ConstIterator it(spectrum.begin());
				for (; it != spectrum.end(); ++it)
				{
					// Write m/z and intensity.
					os << it->getPosition() << " " << it->getIntensity() << "\n";
				}
		
				// Done.
				os.close();		
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_DTAFILE_H

