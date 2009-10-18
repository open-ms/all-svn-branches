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
#include <iostream>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <cstdlib>
using namespace std;
namespace OpenMS
{
	Bzip2Ifstream::Bzip2Ifstream(const char * filename) : n_buffer(0),stream_at_end(false)
	{
		file = fopen( filename, "rb" ); //read binary: always open in binary mode because windows and mac open in text mode
		
		//aborting, ahhh!
		if( !file ) 
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		bzip2file = BZ2_bzReadOpen ( &bzerror, file, 0, 0, NULL, 0 );
		if ( bzerror != BZ_OK ) 
		{
	  	BZ2_bzReadClose ( &bzerror, bzip2file );
	  	throw Exception::ConversionError(__FILE__,__LINE__,__PRETTY_FUNCTION__,"bzip2 compression failed: ");
		}
	}
	
	Bzip2Ifstream::Bzip2Ifstream()
		: file(NULL),bzip2file(NULL),n_buffer(0),bzerror(0),stream_at_end(true)
	{
	}
	
	Bzip2Ifstream::~Bzip2Ifstream()
	{
		BZ2_bzReadClose(&bzerror,bzip2file);
		fclose(file);
	}
	
	size_t Bzip2Ifstream::read(char* s, size_t n)
	{
		if(bzip2file != NULL)
		{
			bzerror = BZ_OK;
		//while is just needed if the whole file should be read at once
		//while ( bzerror == BZ_OK && /* arbitrary other conditions */) 
		//{
  			n_buffer = BZ2_bzRead ( &bzerror, bzip2file, s, n/* size of buf */ );		
	  		if ( bzerror == BZ_OK ) 
	  		{
    			return n_buffer;
    			/* do something with buf[0 .. nBuf-1] */
  			}
		//}
			if ( bzerror != BZ_STREAM_END ) 
			{
   			BZ2_bzReadClose ( &bzerror, bzip2file );
   			throw Exception::ConversionError(__FILE__,__LINE__,__PRETTY_FUNCTION__,"bzip2 compression failed: ");
			} 
			else 
			{
   			BZ2_bzReadClose ( &bzerror, bzip2file);
   			stream_at_end = true;
   			bzip2file = NULL;
   			return n_buffer;
			}
		}
		else
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"no file for decompression initialized");
		}
	}
	
	void Bzip2Ifstream::open(const char* filename)
	{
		if(file != NULL)
		{
			fclose(file);
		}
		if(bzip2file != NULL)
		{
			BZ2_bzReadClose(&bzerror,bzip2file);
		}
		file = fopen( filename, "rb" ); //read binary: always open in binary mode because windows and mac open in text mode
		
		//aborting, ahhh!
		if( !file ) 
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		bzip2file = BZ2_bzReadOpen ( &bzerror, file, 0, 0, NULL, 0 );
		if ( bzerror != BZ_OK ) 
		{
	  	BZ2_bzReadClose ( &bzerror, bzip2file );
	  	throw Exception::ConversionError(__FILE__,__LINE__,__PRETTY_FUNCTION__,"bzip2 compression failed: ");
		}
	}
	
	void Bzip2Ifstream::close()
	{
		if(file != NULL)
		{
			fclose(file);
		}
		if(bzip2file != NULL)
		{
			BZ2_bzReadClose(&bzerror,bzip2file);
		}
		file = NULL;
		bzip2file = NULL;
	}	

} //namespace OpenMS