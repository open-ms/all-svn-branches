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

#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

#include <iostream>

#include <QtGui/QLabel>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>

using namespace std;

namespace OpenMS
{

	PreferencesDialogPage::PreferencesDialogPage(PreferencesManager* manager,QWidget* parent)
	: QWidget(parent),
	manager_(manager),
	help_()
	{
		
	}
	
	PreferencesDialogPage::~PreferencesDialogPage()
	{
		
	}
	
	const string& PreferencesDialogPage::getHelpText() const
	{
		return help_;
	}
	
	void PreferencesDialogPage::setHelpText(const std::string& text)
	{
		help_ = text;	
	}
	
	
	void PreferencesDialogPage::load()
	{
		
	}
	
	void PreferencesDialogPage::save()
	{
		
	}

	QGroupBox* PreferencesDialogPage::addBox(QGridLayout* grid, int row, int column, const QString& label)
	{
		QGroupBox* box = new QGroupBox(label);
		box->setLayout(new QGridLayout(box));
		grid->addWidget(box, row, column);
		return box;
	}
	
	void PreferencesDialogPage::addWidget(QLayout* grid, int row, const QString& label, QWidget* widget) const
	{
		QGridLayout* grid2 = qobject_cast<QGridLayout*>(grid);
		if (grid2!=0)
		{
			QLabel* lab = new QLabel(label);
			grid2->addWidget(lab,row,0);
			grid2->addWidget(widget,row,1);
		}
		else
		{
			cout << __PRETTY_FUNCTION__ << ": Warning, could not cast grid to QGridlayout!" << endl;
		}
	}

	void PreferencesDialogPage::addLayout(QLayout* grid, int row, const QString& label, QLayout* layout) const
	{
		QGridLayout* grid2 = qobject_cast<QGridLayout*>(grid);
		if (grid2!=0)
		{
			QLabel* lab = new QLabel(label);
			grid2->addWidget(lab,row,0);
			grid2->addLayout(layout,row,1);
		}
		else
		{
			cout << __PRETTY_FUNCTION__ << ": Warning, could not cast grid to QGridlayout!" << endl;
		}
	}

	void PreferencesDialogPage::finish(QLayout* grid, int stretch) const
	{
		QGridLayout* grid2 = qobject_cast<QGridLayout*>(grid);
		if (grid2!=0)
		{
			grid2->setColumnStretch(grid2->columnCount(),stretch);	
			grid2->setRowStretch(grid2->rowCount(),stretch);
		}
		else
		{
			cout << __PRETTY_FUNCTION__ << ": Warning, could not cast grid to QGridlayout!" << endl;
		}
	}

} //namespace OpenMS
