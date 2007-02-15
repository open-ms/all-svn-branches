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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_PEAKPICKINGDIALOG_H
#define OPENMS_VISUAL_DIALOGS_PEAKPICKINGDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/PeakPickingDialogTemplate.h>

namespace OpenMS
{

  /**
  	@brief GoTo dialog used to zoom to a m/z and retention time range.
  	
  	@ingroup Dialogs
  */
  class PeakPickingDialog
  	: public QDialog,
  		public Ui::PeakPickingDialogTemplate
  {
    Q_OBJECT

  public:
    PeakPickingDialog(QWidget* parent = 0, Qt::WFlags fl = 0);
    ~PeakPickingDialog();
    void setPeakHeight(float h);
    float getPeakHeight();
    void setPeakHeightMs2(float h2);
    float getPeakHeightMs2();
    void setSignalToNoise(float sn);
    float getSignalToNoise();
    void setFwhm(float fwhm);
    float getFwhm();
    void setOptimization(bool opt);
    bool getOptimization();
  protected slots:
    virtual void resetButton_clicked();
    virtual void startButton_clicked();
  };

}
#endif // OPENMS_VISUAL_DIALOGS_PeakPickingDialog_H

