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
#include <cstdlib>
using namespace std;
namespace OpenMS
{
	Bzip2Ifstream::Bzip2Ifstream(const char * filename) : stream_at_end(false)
	{
		f = fopen( filename, "rb" ); //read binary: always open in binary mode because windows and mac open in text mode
		
		//aborting, ahhh!
		if( !f ) 
		{
  		cout<<"FEHLER in Datei!"<<endl;
  		exit(1);
  		/* handle error */
		}
		
		b = BZ2_bzReadOpen ( &bzerror, f, 0, 0, NULL, 0 );
		if ( bzerror != BZ_OK ) 
		{
	  	BZ2_bzReadClose ( &bzerror, b );
	  	cout<<"FEHLER!in ReadOpen"<<endl;
	  	exit(1);
	  	/* handle error */
		}
	}
	
	size_t Bzip2Ifstream::read(char* s, size_t n)
	{
		bzerror = BZ_OK;
	//while is just needed if the whole file should be read at once
	//	while ( bzerror == BZ_OK && /* arbitrary other conditions */) 
	//	{
  		nBuf = BZ2_bzRead ( &bzerror, b, s, n/* size of buf */ );		
	  	if ( bzerror == BZ_OK ) 
	  	{
    		return nBuf;
    		/* do something with buf[0 .. nBuf-1] */
  		}
		//}
		if ( bzerror != BZ_STREAM_END ) 
		{
   		BZ2_bzReadClose ( &bzerror, b );
   		cout<<"Fehler beim decompressen;"<<endl;
   		exit(1);
   		/* handle error */
		} 
		else 
		{
   		BZ2_bzReadClose ( &bzerror, b );
   		stream_at_end = true;
   		return nBuf;
		}
	}
	

} //namespace OpenMS