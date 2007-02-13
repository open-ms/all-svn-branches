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

 
#ifndef OPENMS_VISUAL_VISUALIZER_PEPTIDEHITVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_PEPTIDEHITVISUALIZER_H


//OpenMS
#include <OpenMS/config.h>
#include <OpenMS/METADATA/PeptideHit.h>
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
@brief Class that displays all meta information for PeptideHit objects

This class provides all functionality to view the meta information of an object of type PeptideHit.
*/
	class PeptideHitVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	   /// Default constructor
		PeptideHitVisualizer(bool editable= FALSE, QWidget *parent =0, const char *name = 0);
		/// Loads the meta data from the object to the viewer.
		void load(PeptideHit &h);

	private slots:
	  /// Saves the changes made to the meta data into the object.
		void store();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject();

	private:  
		
		/** @name Edit fields and buttons
   */
    //@{
		QLineEdit *peptidehit_score_;
		QLineEdit *peptidehit_charge_;
		QLineEdit *peptidehit_score_type_;
		QLineEdit *peptidehit_rank_;
		Q3TextEdit *peptidehit_sequence_;
		
		QPushButton *savebutton_;
		QPushButton *cancelbutton_;
		//@}
		
		/// Pointer to current object	to keep track of the actual object.
		PeptideHit *ptr_;
		/// Copy of current object for restoring the original values
		PeptideHit tempPeptideHit_;
		
		
	};


}
#endif
