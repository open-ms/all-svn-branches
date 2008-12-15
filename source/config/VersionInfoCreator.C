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
// $Maintainer: Chris Bielow
// -------------------------------------------------------------


// read in untouched VersionInfo.C and insert a Macro-Value which will correspond to the current SVN version

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int main( int argc, char** argv )
{
	if (argc!=4)
	{
		cout << "Usage: " << argv[0] << " <infile> <outfile> <svnversionfile>\n";
		return 0;
	}

	string svn_version = "<unknown>";
  ifstream infile_svn (argv[3]);
  if (infile_svn.is_open())
  {
		if (! infile_svn.eof() ) getline (infile_svn,svn_version);
	}
	else
  {
  	cerr << "Unable to open file " << argv[3] << "\n";
  	return 0;
  }


	char * buffer;
	std::streamoff length;
	ifstream infile;
  infile.open (argv[1], ios::binary );
	if (infile.is_open())
  {
	  // get length of file:
	  infile.seekg (0, ios::end);
	  length = infile.tellg();
	  infile.seekg (0, ios::beg);

	  // allocate memory:
	  buffer = new char [length];

	  // read data as a block:
	  infile.read (buffer,length);
	  infile.close();

		ofstream outfile (argv[2],ofstream::binary);
		if (outfile.is_open())
  	{
  		string s_infile = (string) buffer;
  		// replace "//@TAG@" with data and write result
			string to_replace = "//@TAG@";
			size_t pos1 = s_infile.find(to_replace);
			if (pos1==string::npos)
			{
				cerr << "The search string " << to_replace << " was not found in " << argv[1] << "! Who altered it?! Exiting...\n";
				return 0;
			}

			cout << "found tag at pos: " << pos1 << "\n";
			string replace_by= "#define OPENMS_SVNREVISION \"" + svn_version + "\"";
			cout << "Replacing by " << replace_by << "\n";
			s_infile.replace(pos1,to_replace.length(), replace_by);

			// write to outfile
		  outfile.write (s_infile.c_str(), s_infile.length());
			outfile.close();
		}
		else
		{
			cerr << "Unable to open file " << argv[2] << "\n";
	  	return 0;
		}
		delete[] buffer;

	  infile.close();

	}
  else
  {
  	cerr << "Unable to open file " << argv[1] << "\n";
  	return 0;
  }

	return 0;
}