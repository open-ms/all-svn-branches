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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/Spectrum3DCanvasPDP.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>


// Qt
#include <QtGui/QLayout>
#include <QtGui/QRadioButton>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QSpinBox>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		Spectrum3DCanvasPDP::Spectrum3DCanvasPDP(Spectrum3DCanvas* manager, QWidget* parent,  Qt::WFlags f)
				: PreferencesDialogPage(manager,parent,f)
		{
			help_ = "This is the preferences dialog of 3D spectrum!"
								"<br>";
			QGridLayout* grid;
			QLabel * label;
			grid = new QGridLayout(this);
			grid->setMargin(6);
			grid->setSpacing(4);	

			QGroupBox* box = new QGroupBox("Dot coloring",this);

			QGroupBox* coloring_group = new QGroupBox("Color Mode:",box);
			dot_mode_black_ = new QRadioButton("Black",coloring_group);
			dot_mode_gradient_ = new QRadioButton("Gradient",coloring_group);
			dot_gradient_ = new MultiGradientSelector(coloring_group);

			QGroupBox* interpolation_box = new QGroupBox("Interpolation steps",box);
			label = new QLabel("Interpolation steps: ",interpolation_box);
			dot_interpolation_steps_ = new QSpinBox(interpolation_box);
			dot_interpolation_steps_->setMinimum(10);
			dot_interpolation_steps_->setMaximum(1000);
			dot_interpolation_steps_->setSingleStep(1);	
			dot_interpolation_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);

			QGroupBox* shading_group = new QGroupBox("Shade Mode:",box);
			//shading_group->setFrameStyle(QFrame::NoFrame);
			shade_mode_flat_ = new QRadioButton("Flat",shading_group);
			shade_mode_smooth_ = new QRadioButton("Smooth",shading_group);	
			grid->addWidget(box,0,0,1,3);

			box = new QGroupBox("Line Width",this);
			label = new QLabel("Line Width: ",box);
			dot_line_width_ = new QSpinBox(box);
			dot_line_width_->setMinimum(1);
			dot_line_width_->setMaximum(10);
			dot_line_width_->setSingleStep(1);	
			dot_line_width_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,0,1);	

			box = new QGroupBox("Colors",this);
			label = new QLabel("Background color: ",box);
			background_color_ = new ColorSelector(box);
			label = new QLabel("Axes Color: ",box);
			axes_color_ = new ColorSelector(box);
			grid->addWidget(box,1,1);	

			box = new QGroupBox("Data",this);
			label = new QLabel("Reduction Mode:  ",box);
			data_reduction_ = new QComboBox( box);
			data_reduction_->insertItem(0,"Reduction OFF");
			data_reduction_->insertItem(1,"MaxReduction");
			data_reduction_->insertItem(2,"SumReduction");
			label = new QLabel("Displayed Peaks : ",box);
			reduction_diplay_peaks_ = new QSpinBox(box);
			reduction_diplay_peaks_->setMinimum(10000);
			reduction_diplay_peaks_->setMaximum(50000);
			reduction_diplay_peaks_->setSingleStep(5000);	
			grid->addWidget(box,2,1);
		
			load();
		}

		Spectrum3DCanvasPDP::~Spectrum3DCanvasPDP()
		{
			
		}
				
		void Spectrum3DCanvasPDP::load()
		{
			Spectrum3DCanvas* man = static_cast<Spectrum3DCanvas*>(manager_);
			
			if (man->getDotMode()==Spectrum3DCanvas::DOT_GRADIENT)
			{
				dot_mode_gradient_->setChecked(true);
			}
			else
			{
				if (man->getDotMode()==Spectrum3DCanvas::DOT_BLACK)
				{
					dot_mode_black_->setChecked(true);
				}
			}
			if (man->getShadeMode()==Spectrum3DCanvas::SHADE_FLAT)
			{
				shade_mode_flat_->setChecked(true);
			}
			else
			{
				if (man->getShadeMode()==Spectrum3DCanvas::SHADE_SMOOTH)
				{
					shade_mode_smooth_->setChecked(true);
				}
			}
			data_reduction_->setCurrentIndex(data_reduction_->findText(man->getPrefAsString("Preferences:3D:Reduction:Mode").c_str()));	
			reduction_diplay_peaks_->setValue(UnsignedInt(man->getPrefAsInt("Preferences:3D:DisplayedPeaks")));
		
			background_color_->setColor(QColor(man->getPrefAsString("Preferences:3D:BackgroundColor").c_str()));
			dot_gradient_->gradient().fromString(man->getPref("Preferences:3D:Dot:Gradient"));
			dot_interpolation_steps_->setValue(UnsignedInt(man->getPref("Preferences:3D:Dot:InterpolationSteps")));
			dot_line_width_->setValue(UnsignedInt(man->getPref("Preferences:3D:Dot:LineWidth")));
			axes_color_->setColor(QColor(man->getPrefAsString("Preferences:3D:AxesColor").c_str()));
		}
		
		void Spectrum3DCanvasPDP::save()
		{
			Spectrum3DCanvas* man = static_cast<Spectrum3DCanvas*>(manager_);
			if(dot_mode_gradient_->isChecked())
			{
				man->setPref("Preferences:3D:Dot:Gradient",dot_gradient_->gradient().toString());
				man->setDotGradient(dot_gradient_->gradient().toString());
		
				man->setPref("Preferences:3D:Dot:InterpolationSteps",dot_interpolation_steps_->value());
		
				man->setPref("Preferences:3D:Dot:Mode",Spectrum3DCanvas::DOT_GRADIENT);
		
				if(shade_mode_flat_ -> isChecked())
				{
					man->setPref("Preferences:3D:Shade:Mode",Spectrum3DCanvas::SHADE_FLAT);
				}
				else if (shade_mode_smooth_->isChecked())
				{
					man->setPref("Preferences:3D:Shade:Mode",Spectrum3DCanvas::SHADE_SMOOTH);
				}
			}	
			else
			{ 
				if(dot_mode_black_->isChecked())
				{
					man->setPref("Preferences:3D:Dot:Mode",Spectrum3DCanvas::DOT_BLACK);
				}
			} 

			if(data_reduction_->currentText().toAscii().data()=="MaxReduction")
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_MAX);
			}
			else if(data_reduction_->currentText().toAscii().data()=="Reduction OFF")
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_OFF);
			}
			else if(data_reduction_->currentText().toAscii().data()=="SumReduction")
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_SUM);
			}
			
			man->setPref("Preferences:3D:Reduction:Mode", data_reduction_->currentText().toAscii().data());
			man->setPref("Preferences:3D:DisplayedPeaks",	reduction_diplay_peaks_->value());
			
			man->setPref("Preferences:3D:BackgroundColor",background_color_->getColor().name().toAscii().data());
			man->setPref("Preferences:3D:AxesColor",axes_color_->getColor().name().toAscii().data());
			man->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());
		
			man->setDataMode();
			man->repaintAll();	
		

		}
	} // namespace Internal
} //namespace

