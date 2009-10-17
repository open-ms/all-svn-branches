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


//#include <xercesc/util/CompressedInputStream.h>
#include <OpenMS/FORMAT/CompressedInputStream.h>
using namespace xercesc;

namespace OpenMS
{
	CompressedInputStream::CompressedInputStream(const   XMLCh* const    fileName, xercesc::MemoryManager* const manager)
	:bzip2_(new Bzip2Ifstream("/Users/david/Studium/OpenMS/FileCompress2/source/TEST/data/MzMLFile_6_uncompressed.mzML.bz2"))
	{
	}

	CompressedInputStream::CompressedInputStream(const   char* const     fileName, xercesc::MemoryManager* const  manager)
	:bzip2_(new Bzip2Ifstream(static_cast<const char*>(fileName)))
  {
  }
	
	CompressedInputStream::~CompressedInputStream()
	{
		delete bzip2_;
	}
	XMLSize_t CompressedInputStream::readBytes(XMLByte* const  toFill, const XMLSize_t  maxToRead)
	{
    //  Figure out whether we can really read. 
   if(bzip2_->streamEnd())
   {
   	return 0;
   }
   
   unsigned char* fill_it = static_cast<unsigned char*>(toFill);
   XMLSize_t actualRead = (XMLSize_t) bzip2_->read((char*)fill_it, static_cast<const size_t>(maxToRead));
    return actualRead;
	}

	const XMLCh* CompressedInputStream::getContentType() const
	{
    return 0;
	}	
	
} // namespace OpenMS