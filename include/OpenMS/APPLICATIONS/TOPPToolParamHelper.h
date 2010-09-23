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

#ifndef OPENMS_APPLICATIONS_TOPPTOOLPARAMHELPER_H
#define OPENMS_APPLICATIONS_TOPPTOOLPARAMHELPER_H

// OpenMS
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/APPLICATIONS/TOPPIOInfo.h>

// Qt
#include <QtCore/QVector>

namespace OpenMS
{
  /**
    @brief Parameter helper class used in TOPPView and TOPPAS

   Provides simplified methods for common parameter functions.
  */
  class TOPPToolParamHelper
  {
    public:
      /// Refreshes the parameters of this tool, returns true if it has been a change
      static bool refreshParameters(Param&, const String& tool_name, const String& tool_type);

      /// Extracts the types (=different algorithms) of a given TOPP Tool. Empty if no types are defined (=only one algorithm).
      StringList getToolTypes(String tool_name);

      /// Initializes the parameters with standard values from -write_ini, uses the parameters from the old_ini_file if given, returns if parameters have changed (if old_ini_file was given)
      static bool initParam(Param& tool_param, String tool_name, String tool_type, bool show_messagebox_on_error = true, const String& old_ini_file = "");

      /// Fills @p io_infos with the required input file/list parameters.
      static void getInputParameters(const Param&, QVector<TOPPIOInfo>& io_infos);

      /// Fills @p io_infos with the required output file/list parameters.
      static void getOutputParameters(const Param&, QVector<TOPPIOInfo>& io_infos);

      /// Writes @p param to the @p ini_file
      static void writeParam(const Param& param, const String& tool_name, const QString& ini_file);

      /// Extracts input file parameters and tests whether extension matches
      static bool toolAcceptsFileExtension(const Param& target_tool_params, String extension);

      /// Extracts output file parameters and tests against supported_outfile_extensions
      static bool toolDeliversFileExtension(const Param& topp_tool_param, StringList supported_outfile_extensions);

    protected:
      /// Fills @p io_infos with the required input/output file/list parameters. If @p input_params is true, input params are returned, otherwise output params.
      static void getParameters_(const Param&, QVector<TOPPIOInfo>& io_infos, bool input_params);

  };
}
#endif // OPENMS_APPLICATIONS_TOPPTOOLPARAMHELPER_H
