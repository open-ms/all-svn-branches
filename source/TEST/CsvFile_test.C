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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/CsvFile.h>
//#include <OpenMS/FORMAT/FileHandler.h>
//#include <OpenMS/KERNEL/MSExperiment.h>
//#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>


using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(CsvFile, "$Id:$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CsvFile* ptr = 0;
START_SECTION((CsvFile()))
	ptr = new CsvFile;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~CsvFile()))
	delete ptr;
END_SECTION

START_SECTION((void CsvFile(const String& filename, char itemsepartor,bool itemenclose, Int first_n)))
	CsvFile file(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"),',');

	TEST_STRING_EQUAL(file[0],"a1,b1,c1")
	TEST_STRING_EQUAL(file[1],"a2,b2,c2")
	TEST_STRING_EQUAL(file[2],"a3,b3,c3")
	
	CsvFile file2(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),',',true);

	TEST_STRING_EQUAL(file2[0],"\"a1\",\"b1\",\"c1\"")
	TEST_STRING_EQUAL(file2[1],"\"a2\",\"b2\",\"c2\"")
	TEST_STRING_EQUAL(file2[2],"\"a3\",\"b3\",\"c3\"")
END_SECTION


START_SECTION((void fload(const String& filename, char itemsepartor,bool itemenclose, Int first_n)))
	CsvFile file;
	file.fload(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"),',');
	TEST_STRING_EQUAL(file[0],"a1,b1,c1")
	TEST_STRING_EQUAL(file[1],"a2,b2,c2")
	TEST_STRING_EQUAL(file[2],"a3,b3,c3")
	
	file.fload(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),',',true);
	TEST_STRING_EQUAL(file[0],"\"a1\",\"b1\",\"c1\"")
	TEST_STRING_EQUAL(file[1],"\"a2\",\"b2\",\"c2\"")
	TEST_STRING_EQUAL(file[2],"\"a3\",\"b3\",\"c3\"")
END_SECTION

START_SECTION((bool getRow(UInt row, StringList &list)))
	CsvFile file;
	file.fload(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"),',');
	StringList list;
	TEST_EXCEPTION(Exception::InvalidIterator, file.getRow(5, list))
	file.getRow(0,list);
	TEST_EQUAL(list,StringList::create("a1,b1,c1"))
	file.getRow(1,list);
	TEST_EQUAL(list,StringList::create("a2,b2,c2"))
	file.getRow(2,list);
	TEST_EQUAL(list,StringList::create("a3,b3,c3"))
	
		file.fload(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),',',true);
	TEST_EXCEPTION(Exception::InvalidIterator, file.getRow(5, list))
	file.getRow(0,list);
	TEST_EQUAL(list,StringList::create("a1,b1,c1"))
	file.getRow(1,list);
	TEST_EQUAL(list,StringList::create("a2,b2,c2"))
	file.getRow(2,list);
	TEST_EQUAL(list,StringList::create("a3,b3,c3"))
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
