// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DWIDGETPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DWIDGETPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

class QComboBox;

namespace OpenMS
{
	class Spectrum2DWidget;

	namespace Internal
	{

		///Preferences dialog page of a Spectrum2DWidget (internal use only)	
		class Spectrum2DWidgetPDP: public PreferencesDialogPage
		{
			Q_OBJECT
			
			public:
				/// Constructor
				Spectrum2DWidgetPDP( Spectrum2DWidget* manager, QWidget* parent = 0, const char* name = "Spectrum2DWidgetPDP", Qt::WFlags f = 0);
				///  Destructor
				virtual ~Spectrum2DWidgetPDP();
				// Docu in base class
				virtual void load();
				// Docu in base class
				virtual void save();
			protected:
				PreferencesDialogPage* canvas_;
			  QComboBox* axis_mapping_;
		};
	
	} //namespace Internal
	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUM2DWIDGETPDP_H

