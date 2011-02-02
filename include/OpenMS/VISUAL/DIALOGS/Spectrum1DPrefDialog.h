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


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM1DPREFDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM1DPREFDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_Spectrum1DPrefDialog.h>

namespace OpenMS 
{
	namespace Internal
	{
		///Preferences dialog for Spectrum1DWidget
		class OPENMS_DLLAPI Spectrum1DPrefDialog
			: public QDialog,
	  		public Ui::Spectrum1DPrefDialogTemplate
		{
			Q_OBJECT
			
			public:
				///Constructor
				Spectrum1DPrefDialog(QWidget * parent);
		};
	}
}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM1DPREFDIALOG_H

