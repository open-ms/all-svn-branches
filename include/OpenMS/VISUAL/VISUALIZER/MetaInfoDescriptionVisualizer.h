// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: stefan_heess $
// --------------------------------------------------------------------------

 
#ifndef OPENMS_VISUAL_VISUALIZER_METAINFODESCRIPTIONVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_METAINFODESCRIPTIONVISUALIZER_H


//OpenMS
#include <OpenMS/config.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

//QT
#include <q3textedit.h>
#include <qpushbutton.h>
#include <QLabel>
#include <iostream>
#include <qwidget.h>
#include <qvalidator.h>

class QLabel;
class QLineEdit;


namespace OpenMS {
/**
@brief Class that displays all meta information for MetaInfoDescription objects

This class provides all functionality to view the meta information of an object of type MetaInfoDescription.
*/
	class MetaInfoDescriptionVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	   /// Default constructor
		MetaInfoDescriptionVisualizer(bool editable= FALSE, QWidget *parent =0, const char *name = 0);
		/// Loads the meta data from the object to the viewer.
		void load(MetaInfoDescription &h);

	private slots:
	  /// Saves the changes made to the meta data into the object.
		void store();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject();

	private:  
		
		/** @name Some edit fields
		*/
		//@{
   	QLineEdit *metainfodescription_name_;
		QTextEdit *metainfodescription_comment_;
		//@}
		
		/** @name Some buttons
		*/
		//@{
		QPushButton *savebutton_;
		QPushButton *cancelbutton_;
		//@}
		
		/// Pointer to current object	to keep track of the actual object.
		MetaInfoDescription *ptr_;
		/// Copy of current object for restoring the original values
		MetaInfoDescription tempMetaInfoDescription_;
		
		
	};


}
#endif
