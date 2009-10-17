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
#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/CompressedInputStream.h>

#include <xercesc/internal/MemoryManagerImpl.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUniDefs.hpp>


using namespace xercesc;
namespace OpenMS
{
	CompressedInputSource::CompressedInputSource( const XMLCh* const basePath, const XMLCh* const relativePath, xercesc::MemoryManager* const manager)
    : xercesc::InputSource(manager)
	{
	    //
	    //  If the relative part is really relative, then weave it together
  	  //  with the base path. If not, just take the relative path as the
   	 //  entire path.
   	 //
   	 if(xercesc::XMLPlatformUtils::isRelative(relativePath, manager))
   	 {
   	 	XMLCh* tmpBuf = xercesc::XMLPlatformUtils::weavePaths(basePath, relativePath, manager);
   	  setSystemId(tmpBuf);
      manager->deallocate(tmpBuf); //delete [] tmpBuf;
    	}
    	else
    	{
    	  XMLCh* tmpBuf = xercesc::XMLString::replicate(relativePath, manager);
        xercesc::XMLPlatformUtils::removeDotSlash(tmpBuf, manager);
        setSystemId(tmpBuf);
        manager->deallocate(tmpBuf);//delete [] tmpBuf;
    	}

	}

	CompressedInputSource::CompressedInputSource(const XMLCh* const filePath, MemoryManager* const manager)
    : xercesc::InputSource(manager)
	{

    	//
    	//  If the path is relative, then complete it acording to the current
    	//  working directory rules of the current platform. Else, just take
    	//  it as is.
    	//
    	if (xercesc::XMLPlatformUtils::isRelative(filePath, manager))
    	{
        XMLCh* curDir = xercesc::XMLPlatformUtils::getCurrentDirectory(manager);

        XMLSize_t curDirLen = XMLString::stringLen(curDir);
        XMLSize_t filePathLen = XMLString::stringLen(filePath);
        XMLCh* fullDir = (XMLCh*) manager->allocate
        (
            (curDirLen + filePathLen + 2) * sizeof(XMLCh)
        );//new XMLCh [ curDirLen + filePathLen + 2];

        XMLString::copyString(fullDir, curDir);
        fullDir[curDirLen] = chForwardSlash;
        XMLString::copyString(&fullDir[curDirLen+1], filePath);
        
        XMLPlatformUtils::removeDotSlash(fullDir, manager);
        XMLPlatformUtils::removeDotDotSlash(fullDir, manager);

        setSystemId(fullDir);

        manager->deallocate(curDir);//delete [] curDir;
        manager->deallocate(fullDir);//delete [] fullDir;
    	}
     	else
    	{
        XMLCh* tmpBuf = XMLString::replicate(filePath, manager);
        XMLPlatformUtils::removeDotSlash(tmpBuf, manager);
        setSystemId(tmpBuf);
      	manager->deallocate(tmpBuf);//delete [] tmpBuf;
  	  }

	}

	CompressedInputSource::~CompressedInputSource()
	{
	}
	

BinInputStream* CompressedInputSource::makeStream() const
{
    CompressedInputStream* retStrm = new CompressedInputStream(getSystemId(), getMemoryManager());
    if (!retStrm->getIsOpen())
    {
       delete retStrm;
        return 0;
    }
    return retStrm;
}	
	
	
} // namespace OpenMS