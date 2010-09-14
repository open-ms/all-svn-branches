// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOPPIOINFO_H
#define OPENMS_APPLICATIONS_TOPPIOINFO_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

namespace OpenMS
{
  /// Stores the information for input/output files/lists
  struct TOPPIOInfo
  {
    ///Standard constructor
    TOPPIOInfo()
      :	type(IOT_FILE),
        param_name(),
        valid_types()
    {
    }

    ///Copy constructor
    TOPPIOInfo(const TOPPIOInfo& rhs)
      :	type(rhs.type),
        param_name(rhs.param_name),
        valid_types(rhs.valid_types)
    {
    }

    ///The type
    enum IOType
    {
      IOT_FILE,
      IOT_LIST
    };

    ///Comparison operator
    bool operator< (const TOPPIOInfo& rhs) const
    {
      if (type != rhs.type)
      {
        return type == IOT_FILE;
      }
      else
      {
        return param_name.compare(rhs.param_name) < 0;
      }
    }

    ///Assignment operator
    TOPPIOInfo& operator= (const TOPPIOInfo& rhs)
    {
      type = rhs.type;
      param_name = rhs.param_name;
      valid_types = rhs.valid_types;

      return *this;
    }

    ///The type of the parameter
    IOType type;
    ///The name of the parameter
    String param_name;
    ///The valid file types for this parameter
    StringList valid_types;
  };
}
#endif // OPENMS_APPLICATIONS_TOPPIOINFO_H
