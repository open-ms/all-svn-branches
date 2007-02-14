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

#ifndef OPENMS_APPLICATIONS_TOPPVIEWBASE_H
#define OPENMS_APPLICATIONS_TOPPVIEWBASE_H

//OpenGL
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/VISUAL/DIALOGS/OpenDialog.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

//STL
#include <map>

//QT
#include <q3mainwindow.h>
#include <qworkspace.h>
#include <QCloseEvent>
#include <QLabel>
#include <Q3PopupMenu>
#include <Q3ActionGroup>

class Q3Action;
class QComboBox;
class QToolButton;
class QMenuBar;
class QLabel;
class QRadioButton;
class Q3ActionGroup;
class Q3PopupMenu;
class Q3ListView;
class Q3ListViewItem;

namespace OpenMS
{
  class MultiGradientSelector;
  class Spectrum1DWindow;
  class Spectrum2DWindow;
  class Spectrum3DWindow;
  class SpectrumWindow;

  /**
  	@brief MDI window of TOPPView tool
  	
  	@todo Add preferences for layers (Marc)
  	@todo Reimplement PreferenceDialogPages (Marc)
  	@todo Use right mouse button for navigation in data (Marc)
		@todo: Use log-scaled intensities intensity-cutoff plot (Marc)
		
  */
  class TOPPViewBase : public Q3MainWindow, public PreferencesManager
  {
      Q_OBJECT

    public:
      /// Access is possible only through this method as TOPPViewBase is a singleton
      static TOPPViewBase* instance();

      /**
      	@brief Opens and displays a spectrum form a file
      	
      	@param filename The file to open
      	@param as_new_window If the data is displayed in the current window or in a new window
      	@param maps_as_2d If maps are displayed 2D or 3D
      	@param maximize If the window the new file was added to should be displayed maximized
      	@param use_mower If a mower should be used to suppress noise in the data
      	@param force_type File type to force
      */
      void addSpectrum(const String& filename, bool as_new_window=true, bool maps_as_2d=true, bool maximize=false, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER, FileHandler::Type force_type=FileHandler::UNKNOWN);
      /**
      	@brief Opens and displays a spectrum form the database
      	
      	@param db_id The id in the database
      	@param as_new_window If the data is displayed in the current window or in a new window
      	@param maps_as_2d If maps are displayed 2D or 3D
      	@param maximize If the window the new file was added to should be displayed maximized
      	@param use_mower If a mower should be used to suppress noise in the data
      */
      void addDBSpectrum(UnsignedInt db_id, bool as_new_window=true, bool maps_as_2d=true, bool maximize=false, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);

      /// maximizes the size of the active window
      void maximizeActiveSpectrum();
      /// opens all the files that are inside the handed over iterator range
      template <class StringListIterator>
      void loadFiles(const StringListIterator& begin, const StringListIterator& end)
      {
        //use mower?
        OpenDialog::Mower mow = OpenDialog::NO_MOWER;
        if ( getPrefAsString("Preferences:MapIntensityCutoff")=="Noise Estimator")
        {
          mow = OpenDialog::NOISE_ESTIMATOR;
        }
				
				bool last_was_plus = false;
        for (StringListIterator it=begin; it!=end; ++it)
        {
          if (*it=="+")
          {
          	last_was_plus = true;
          	continue;
        	}
        	else if (!last_was_plus)
        	{
        		addSpectrum(*it,true,getPrefAsString("Preferences:DefaultMapView")=="2D",true,mow);
        	}
        	else
        	{
        		last_was_plus = false;
        		addSpectrum(*it,false,getPrefAsString("Preferences:DefaultMapView")=="2D",true,mow);
        	}
        }
        maximizeActiveSpectrum();
        for (int i=0; i<tab_bar_->count(); ++i)
        {
        	if (tab_bar_->tabData(i).toInt() == qlonglong(ws_->activeWindow()))
        	{
        		tab_bar_->setCurrentIndex(i);
        		break;
        	}
        }
      }
      /// returns selected peaks of the active spectrum framed by \c layer_index_.begin() and the last peak BEFORE \c layer_index_.end();
      std::vector<MSExperiment<>::SpectrumType::Iterator> getActiveSpectrumSelectedPeaks();

      /**
      	@brief Loads the preferences from the filename given.
      	
      	If the filename is empty, the application name + ".ini" is used as filename
      */
      void loadPreferences(std::string filename="");
      /// stores the preferences (used when this window is closed)
      void savePreferences();

    public slots:
      /// shows the dialog for opening spectra from file or the database
      void openSpectrumDialog();
      /// shows the goto dialog
      void gotoDialog();
      /// shows the preferences dialog
      void preferencesDialog();
      /// lets the user edit the preferences of a layer
      void layerPreferencesDialog();
      /// Lets the user change the intensity distribution of a layer
      void layerIntensityDistribution();
			/// Changes the axis visibility of the current window
			void changeAxisVisibility();
      /// lets the user edit the meta data of a layer
      void editMetadata();
      /// saves the contents of the active window
      void saveImage();
      /// saves the content of active window to an image
      void print();
      /// closes the active window
      void closeFile();
      /// saves the current view of the current layer
      void saveLayer();
      /// updates the toolbar, when the active window changes
      void updateToolbar(QWidget* widget);
      /// adapts the layer bar to the active window
      void updateLayerbar();
      /// brings the tab corresponding to the active window in front
      void updateTabBar(QWidget* w);
      /// tile the open windows vertically
      void tileVertical();
      /// tile the open windows horizontally
      void tileHorizontal();
      /// Links or unlinks two spectra (for zooming)
      void linkActiveTo(const QString&);
      /**
      	@brief Shows a status message in the status bar.
      	
      	If time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
      	Otherwise the message is displayed for @p time ms.
      */
      void showStatusMessage(std::string msg, OpenMS::UnsignedInt time);
      /// shows m/z, intensity and rt in the status bar
      void showCursorStatus(double mz, double intensity, double rt);

      /// Shows a list all selected peaks
      void showPeaklistActiveSpectrum();
      /// Picks peaks in the active spectrum
      void pickActiveSpectrum();
      /// Filter baseline
      void baselineFilteringActiveSpectrum();
      /// Smooth data in the active spectrum
      void smoothActiveSpectrum();
      /// Finds features in the active spectrum
      void findFeaturesActiveSpectrum();
			
    protected slots:
    	/// slot for layer manager selection change
    	void layerSelectionChange();
    	/// slot for layer manager context menu
    	void layerContextMenu(Q3ListViewItem* item, const QPoint& pos, int col);
    	/// signal for layer manager visibility change (check box)
    	void layerVisibilityChange(Q3ListViewItem* item, const QPoint& pnt, int col);
    	
      void closeFileByTab(int);
      void focusByTab(int);
      void removeWidgetFromBar(QObject*);
      void openRecentFile(int i);

      /** @name Toolbar slots
      */
      //@{
      void setActionMode(Q3Action*);
      void setDrawMode1D(Q3Action*);
      void setIntensityMode(Q3Action* a);
      void showGridLines(bool);
      void showPoints(bool);
      void showSurface(bool);
      void showContours(bool);
      void resetZoom();
      //@}

      ///use this event to do the cleanup
      virtual void closeEvent(QCloseEvent * e);
      /// Call whenever a window is closed
      virtual void windowClosed();

    protected:
      ///singleton instance
      static TOPPViewBase* instance_;
      ///not accessable as this class is a singleton
      TOPPViewBase(QWidget* parent=0, Qt::WindowFlags f=0);
      ///not accessable as this class is a singleton
      ~TOPPViewBase();

      /// Adds a tab for the window in the tabbar
      void addTab_(SpectrumWindow*, const String&);
      /// connect the slots/signals for status messages and mode changes (paint or mouse mode)
      void connectWindowSignals_(SpectrumWindow* sw);
      ///returns a pointer to the active SpectrumWindow (0 if none is active)
      SpectrumWindow*  activeWindow_() const;
      ///returns a pointer to the active SpectrumCanvas (0 if none is active)
      SpectrumCanvas*  activeCanvas_() const;
      ///returns a pointer to the active Spectrum1DWindow (0 the active window is no Spectrum1DWindow or there is no active window)
      Spectrum1DWindow* active1DWindow_() const;
      ///returns a pointer to the active Spectrum2DWindow (0 the active window is no Spectrum2DWindow or there is no active window)
      Spectrum2DWindow* active2DWindow_() const;
      ///returns a pointer to the active Spectrum3DWindow (0 the active window is no Spectrum2DWindow or there is no active window)
      Spectrum3DWindow* active3DWindow_() const;
      ///Estimates the noise by evaluating 10 random scans of MS level 1
      float estimateNoise_(const SpectrumCanvas::ExperimentType& exp);

      // Docu in base class
      virtual PreferencesDialogPage* createPreferences(QWidget* parent);

      /// Layer mangment bar
      Q3ToolBar* layer_bar_;
      /// Layer mangment widget
      Q3ListView* layer_manager_;

      /// Creates the toolbars and connects the signals and slots
      void createToolBar_();

      /** @name Toolbar members
      */
      //@{
      Q3ToolBar* tool_bar_;
      //common actions
      Q3ActionGroup* action_modes_;
      Q3Action* am_zoom_;
      Q3Action* am_translate_;
      Q3Action* am_select_;
      Q3Action* am_measure_;
      //common intensity modes
      Q3ActionGroup* intensity_modes_;
      Q3Action* im_none_;
      Q3Action* im_log_;
      Q3Action* im_percentage_;
      Q3Action* im_snap_;
      //common buttons
      QToolButton* reset_zoom_button_;
      QToolButton* grid_button_;
      QToolButton* print_button_;
      //1D specific stuff
      Q3ToolBar* tool_bar_1d_;
      Q3ActionGroup* draw_modes_;
      Q3Action* dm_peaks_1d_;
      Q3Action* dm_rawdata_1d_;
      QComboBox* link_box_;
      //2D specific stuff
      Q3ToolBar* tool_bar_2d_;
      QToolButton* dm_points_2d_;
      QToolButton* dm_surface_2d_;
      QToolButton* dm_contours_2d_;
      Q3ActionGroup* draw_modes_2d_;
      Q3Action* dm2_points_2d_;
      Q3Action* dm2_surface_2d_;
      Q3Action* dm2_contours_2d_;
      //@}

      /// Main workspace
      QWorkspace* ws_;

      ///Tab bar. The address of the corresponding window to a tab is stored as an int in tabData()
      EnhancedTabBar* tab_bar_;
      ///map between window address as int and window pointer used for toolbar and tabbar
      std::map<qlonglong,SpectrumWindow*> id_map_;

      /// Label for messages in the status bar
      QLabel* message_label_;
      /// m/z label for messages in the status bar
      QLabel* mz_label_;
      /// Intensity label for messages in the status bar
      QLabel* int_label_;
      /// RT label for messages in the status bar
      QLabel* rt_label_;

      /**
        @brief Map that stores linked pairs of 1D windows (uses int value of addresses to identify the widgets).
      	
      	Each link is stored twice (both directions).
      */
      std::map<qlonglong,qlonglong> link_map_;

      ///unlinks active spectrum (for zooming)
      void unlinkActive_();

      ///adds a Filename to the recent files
      void addRecentFile_(const String& filename);

      ///update the recent files menu
      void updateRecentMenu_();

      /// check if all avaiable preferences get set by the .ini file. If there are some missing entries fill them with default values.
      void checkPreferences_();
      /// list of the recently opened files
      std::vector<String> recent_files_;

      /// Pointer to "Tools" menu: so that derived classes can add stuff into it
      Q3PopupMenu* tools_menu_;
      /// pointer to the recent files menu
      Q3PopupMenu* recent_menu_;

  }
  ; //class

} //namespace

#endif // OPENMS_APPLICATIONS_TOPPVIEWBASE_H
