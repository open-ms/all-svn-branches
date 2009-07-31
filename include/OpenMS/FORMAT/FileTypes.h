// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FILETYPES_H
#define OPENMS_FORMAT_FILETYPES_H

namespace OpenMS
{
	struct FileTypes
	{
	/**
		@brief Centralizes the file types recognized by FileHandler. FileType seperate from FileHandler to avoid circular inclusions by DocumentIdentifier, ExperimentalSettings and FileHandler and respective fileclasses (e.g. DTA2DFile). See also: FileHandler::nameToType, FileHandler::typeToName and FileHandler::NameOfTypes .

				@ingroup FileIO
	*/
		//if you change here, do not forget to change FileHandler::NameOfTypes[]
		enum Type
		{
			UNKNOWN,        		///< Unknown file extension
			DTA,            		///< DTA file (.dta)
			DTA2D,          		///< DTA2D file (.dta2d)
			MZDATA,         		///< MzData file (.MzData)
			MZXML,          		///< MzXML file (.MzXML)
			FEATUREXML,     		///< %OpenMS feature file (.featureXML)
			ANDIMS,         		///< ANDI\\MS file (.cdf)
			IDXML,  						///< %OpenMS identification format (.idXML)
			CONSENSUSXML,  			///< %OpenMS consensus map format (.consensusXML)
			MGF,								///< Mascot Generic Format (.mgf)
			INI,          			///< %OpenMS parameters file (.ini)
			TRANSFORMATIONXML,  ///< Tranformation description file (.trafoXML)
			MZML,								///< MzML file (.mzML)
			MS2,								///< MS2 file (.ms2)
			PEPXML,							///< TPP pepXML file (.pepXML)
			MZIDENTML,				  ///< mzIdentML (HUPO PSI AnalysisXML format) (.mzIdentML)
			GELML,							///< GelML (HUPO PSI format) (.GelML)
			TRAML,							///< TraML (HUPO PSI format) for transitions (.TraML)
			MSP,								///< NIST spectra library file format (.msp)
			AUTOEXECUTE,     		///< AutoExecute file (.txt)
			XMASS,           		///< XMass Analysis file (fid)			
			SIZE_OF_TYPE    		///< No file type. Simply stores the number of types
		};
	}; //struct
} //namespace OpenMS

#endif //OPENMS_FORMAT_FILETYPES_H
