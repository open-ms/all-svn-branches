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

#ifndef OPENMS_FORMAT_BZIP2_IFSTREAM_H
#define	OPENMS_FORMAT_BZIP2_IFSTREAM_H

#include "/opt/local/var/macports/software/bzip2/1.0.5_2/opt/local/include/bzlib.h"
#include <istream>

namespace OpenMS
{
	class Bzip2Ifstream 
	{
		public: 
			Bzip2Ifstream(const char * filename);
			Bzip2Ifstream();
			~Bzip2Ifstream();
			
			//operator>>();
			size_t read(char* s, size_t n);
			bool streamEnd() const;
			bool isOpen() const;
			void open(const char* filename);
			void close();
		protected:
			FILE*   file;
			BZFILE* bzip2file;
			size_t     n_buffer;
			int     bzerror;
			bool stream_at_end;
			
			//not implemented
			Bzip2Ifstream(const Bzip2Ifstream& bzip2);
			Bzip2Ifstream& operator=(const Bzip2Ifstream& bzip2);
	};
	
	inline bool Bzip2Ifstream::isOpen() const
	{
		return (file);
	}
	
		inline bool Bzip2Ifstream::streamEnd() const
	{
		return stream_at_end;
	}

} //namespace OpenMS
#endif //OPENMS_FORMAT_BZIP2_IFSTREAM_H