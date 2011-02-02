// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DPREFDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DPREFDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_Spectrum2DPrefDialog.h>

namespace OpenMS 
{
	namespace Internal
	{
		///Preferences dialog for Spectrum2DWidget
		class OPENMS_DLLAPI Spectrum2DPrefDialog
			: public QDialog,
	  		public Ui::Spectrum2DPrefDialogTemplate
		{
			Q_OBJECT
			
			public:
				///Constructor
				Spectrum2DPrefDialog(QWidget * parent);
		};
	}
}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM2DPREFDIALOG_H

