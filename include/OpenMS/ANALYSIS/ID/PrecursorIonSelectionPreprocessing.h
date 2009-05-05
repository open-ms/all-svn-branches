// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_ANALYSIS_ID_PRECURSORIONSELECTIONPREPROCESSING_H
#define OPENMS_ANALYSIS_ID_PRECURSORIONSELECTIONPREPROCESSING_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/AASequence.h>


#include <algorithm>
#include <cmath>
#include <set>
#include <fstream>
namespace OpenMS
{
	
  /**
		 @brief This class implements the database preprocessing needing for precursor ion selection.

		 

		 @htmlinclude OpenMS_PrecursorIonSelectionPreprocessing.parameters

	*/
  class OPENMS_DLLAPI PrecursorIonSelectionPreprocessing : public DefaultParamHandler
  {
  public:
      PrecursorIonSelectionPreprocessing();
      PrecursorIonSelectionPreprocessing(const PrecursorIonSelectionPreprocessing& source);
      ~PrecursorIonSelectionPreprocessing();
		
		PrecursorIonSelectionPreprocessing& operator = (const PrecursorIonSelectionPreprocessing& source);
		
		//const std::set<AASequence>& getSequences() const;

		const std::map<String,std::vector<DoubleReal> >& getProtMasses() const;


		const std::vector<DoubleReal> & getMasses(String acc);
	

		/**
		 *	@brief Calculates tryptric peptide masses of a given database and stores masses and peptide sequences
		 *	
		 *	@param db_path Path to database file (fasta)
		 *	@param save Flag if preprocessing should be stored.
		 *
		 *	@throws Exception::FileNotFound if file with preprocessing or db can't be found
		 *  @throws Exception::UnableToCreateFile if preprocessing file can't be written
		 */
		void dbPreprocessing(String db_path,bool save=true);

		
		/**
		 *	@brief Loads tryptric peptide masses of a given database.
		 *
		 *	@throws Exception::FileNotFound if file with preprocessing can't be found
		 *	@throws Exception::InvalidParameter if precursor_mass_tolerance_unit is ppm and
		 *  file containing bin borders can't be found
		 */
		void loadPreprocessing();

		/// get the weighted frequency of a mass
		DoubleReal getWeight(DoubleReal mass);
  protected:
		/// saves the preprocessed db
		void savePreprocessedDB_(String db_path,String path);
		/// loads the preprocessed db
		void loadPreprocessedDB_(String path);

		/// all tryptic masses of the distinct peptides in the database
		std::vector<DoubleReal> masses_;
    /// the sequences of the tryptic peptides
    std::set<AASequence> sequences_;
    /// stores masses of tryptic peptides for proteins, key is the accession number
    std::map<String,std::vector<DoubleReal> > prot_masses_;
    /// the masses of the bins used for preprocessing (only used if bins are not equidistant, i.e. with ppm)
    std::vector<DoubleReal> bin_masses_;
		/// counter for the bins
		std::vector<UInt> counter_;
		/// maximal relative frequency of a mass
    UInt f_max_;

  };
}
    
#endif //#ifndef OPENMS_ANALYSIS_ID_PRECURSORIONSELECTIONPREPROCESSING_H
