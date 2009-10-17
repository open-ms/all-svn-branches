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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_COMPRESSEDINPUTSTREAM_H
#define OPENMS_FORMAT_COMPRESSEDINPUTSTREAM_H

#include <xercesc/util/BinInputStream.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

namespace OpenMS
{
	class CompressedInputStream
		:	public xercesc::BinInputStream
	{
		public:
			//based on LocalInputSource
			CompressedInputStream(const   XMLCh* const    fileName, xercesc::MemoryManager* const manager = xercesc::XMLPlatformUtils::fgMemoryManager);

   		CompressedInputStream(const   char* const     fileName, xercesc::MemoryManager* const  manager = xercesc::XMLPlatformUtils::fgMemoryManager);	 
   		~CompressedInputStream();
   		
   		 bool getIsOpen() const;
    	// -----------------------------------------------------------------------
    	//  Implementation of the input stream interface
    	// -----------------------------------------------------------------------
    	virtual XMLFilePos curPos() const;

   	 	virtual XMLSize_t readBytes(XMLByte* const  toFill, const XMLSize_t maxToRead);

	    virtual const XMLCh* getContentType() const;

   		
    private:
    	Bzip2Ifstream* 	bzip2_;
    	XMLSize_t       fCurIndex;
	};
	
	inline XMLFilePos CompressedInputStream::curPos() const
	{
    return fCurIndex;
	}
	
	inline bool CompressedInputStream::getIsOpen() const
	{
			return bzip2_->isOpen();
	}
} // namespace OpenMS

#endif // OPENMS_FORMAT_COMPRESSEDINPUTSTREAM_H