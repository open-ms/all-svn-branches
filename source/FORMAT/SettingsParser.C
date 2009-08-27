// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Import settings
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SettingsParser.h>

namespace OpenMS
{

	SettingsParser::SettingsParser()
	{
	}

	int SettingsParser::findIndex_(const String& value, const std::string* table, const Size tableSize) const
	{
	  int index = tableSize;
	  for(; index!=0 && table[index]!=value; --index)
	  {
	  }  
	  return index;
	}
		
	void SettingsParser::stringToIntVector_(const String& str, std::vector<Int>& list) const
	{
	  std::vector<String> strings;
	  str.split(',', strings);
	  list.clear();
	  list.resize(strings.size());
	  
	  for(Size index=0; index<strings.size(); ++index)
	    list[index] = strings[index].toInt();
	}
	
  void SettingsParser::addMetaInfo_(MetaInfoInterface& metaInfo, const String& name, const String& value) const
  {
    std::vector<String> strings;
    name.split('.', strings);
    if(strings.size() > 1)
      metaInfo.setMetaValue(strings[strings.size()-1], value);
  }    	
    
} // namespace OpenMS

