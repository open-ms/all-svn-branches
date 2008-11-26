// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

#include <fstream>

// STL includes
#include <iostream>


namespace OpenMS
{
	namespace Internal
	{
		/**
		@brief XML Handler for a consensusXML.
		*/
		class ConsensusXMLHandler
			: public XMLHandler
		{
		 public:

			/// Constructor
			ConsensusXMLHandler(ConsensusMap& consensus_map , const String& filename, const String& version)
				: XMLHandler(filename, version),
					consensus_map_(&consensus_map),
					act_cons_element_(),
					last_meta_(0)
			{
			}

			/// Destructor
			virtual ~ConsensusXMLHandler()
			{
			}

			// Docu in base class
			virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

			// Docu in base class
			virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

			// Docu in base class
			virtual void characters(const XMLCh* const chars, const unsigned int length);

			///Writes the contents to a stream
			void writeTo(std::ostream& os);

			///Sets the options
			void setOptions(const PeakFileOptions& options);

		 protected:
			/// Options that can be set
			PeakFileOptions options_;
			///@name Temporary variables for parsing
			//@{
			ConsensusMap* consensus_map_;
			ConsensusFeature act_cons_element_;
			DPosition<2> pos_;
			DoubleReal it_;
			UInt last_map_;

			/// Pointer to last read object as a MetaInfoInterface, or null.
			MetaInfoInterface* last_meta_;
			/// Temporary protein ProteinIdentification
			ProteinIdentification prot_id_;
			/// Temporary peptide ProteinIdentification
			PeptideIdentification pep_id_;
			/// Temporary protein hit
			ProteinHit prot_hit_;
			/// Temporary peptide hit
			PeptideHit pep_hit_;
			/// Map from protein id to accession
			std::map<String,String> proteinid_to_accession_;


			//@}
		};

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H
