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

 
#ifndef OPENMS_VISUAL_VISUALIZER_TAGGINGVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_TAGGINGVISUALIZER_H

//OpenMS
//#include <OpenMS/config.h>
#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

//QT
//#include <QtGui/qtextedit.h>
//#include <QtGui/qpushbutton.h>
//#include <QtGui/QLabel>
//#include <iostream>
//#include <QtGui/qwidget.h>
//#include <QtGui/qvalidator.h>

class QLabel;
class QLineEdit;
class QComboBox;
class QDoubleValidator;

namespace OpenMS {
/**
@brief Class that displays all meta information of tagging objects.

This class provides all functionality to view the meta information of an object of type Tagging.
*/
	class TaggingVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		TaggingVisualizer(bool editable= FALSE, QWidget *parent =0);
		/// Loads the meta data from the object to the viewer.
		void load(Tagging &t);

	private slots:
	  /// Saves the changes made to the meta data into the object.
		void store();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject();

	private:  
		
	  /// A validator to check the input for the mass
		QDoubleValidator *massvali_;
		
		/// A Validator to check the input for the shift
		QDoubleValidator *shiftvali_;
			
		/// Sets the fields and comboboxes with current values
		void updateTag_();
		
		/** @name Edit fields and buttons
   */
    //@{
		QLineEdit *treatmenttype_;
		QTextEdit *treatmentcomment_;
		QLineEdit *modificationname_;
		QLineEdit *modificationmass_;
		QComboBox *modificationspecificity_;
		QLineEdit *modificationAA_;
		QLineEdit *taggingmass_shift_;
		QComboBox *taggingvariant_;
		
		
		QPushButton *savebutton_;
		QPushButton *cancelbutton_;
		//@}
		
		/// Pointer to current object	 to keep track of the actual object.
		Tagging *ptr_;
		/// Copy of current object for restoring the original values.
		Tagging temptag_;
		
		
	};


}
#endif
