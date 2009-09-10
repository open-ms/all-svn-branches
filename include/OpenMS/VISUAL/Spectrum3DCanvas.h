// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public2
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
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DCANVAS_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>

class QPainter;
class QGLWidget;
class QResizeEvent;

namespace OpenMS
{
  class Spectrum3DOpenGLCanvas;
  
  /**
    @brief Canvas for 3D-visualization of peak map data
		
		The Spectrum3DCanvas uses the helper class Spectrum3DOpenGLCanvas for the actual 3D rendering.
		Deriving Spectrum3DCanvas directly from QGLWidget is not possible due to the "Deadly Diamond" shape
		of inheritence.
		
		@image html Spectrum3DWidget.png
		
		@htmlinclude OpenMS_Spectrum3DCanvas.parameters
		
    @ingroup SpectrumWidgets
  */	
  class OPENMS_DLLAPI Spectrum3DCanvas
  	: public SpectrumCanvas
  {
    Q_OBJECT
    
    friend class Spectrum3DOpenGLCanvas;
    
	  public:
			
	    /// Constructor
	    Spectrum3DCanvas(const Param& preferences, QWidget* parent = 0);	
	    /// Destructor
	    virtual  ~Spectrum3DCanvas();
	    
	  	///Different shade modes
	    enum ShadeModes 
	    {
				SHADE_FLAT = 0,            
				SHADE_SMOOTH = 1         
	    };
			
			///Enumerate all avaiable paint styles
			enum ViewModes 
			{
				VM_2D,
				VM_3D
			};

			///Enumerate all avaiable actions
			enum Actions 
			{
			  A_CAMERA_RESET,
				A_CAMERA_MOVEUP,
				A_CAMERA_MOVEDOWN,
				A_CAMERA_MOVELEFT,
				A_CAMERA_MOVERIGHT,
				A_CAMERA_ZOOMIN,
				A_CAMERA_ZOOMOUT,
			  A_DATA_RESET,
				A_DATA_MOVEUP,
				A_DATA_MOVEDOWN,
				A_DATA_MOVELEFT,
				A_DATA_MOVERIGHT,
				A_DATA_ZOOMIN,
				A_DATA_ZOOMOUT
			};
			
	    ///returns the Spectrum3DOpenGLcanvas     
	    Spectrum3DOpenGLCanvas* openglwidget();
	    Spectrum3DOpenGLCanvas* openglwidget() const;
	    	    
	    ///@name Remplemented Qt events
	    //@{
	    void resizeEvent(QResizeEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
	    //@}
	    /// Returns if the legend is shown
	    bool isLegendShown() const;
	    ///Shows/hides the legend
	    void showLegend(bool);
	    ///pointer to the SpectrumOpenGLCanvas implementation
	    Spectrum3DOpenGLCanvas* openglcanvas_;
			
			// docu in base class
			virtual void showCurrentLayerPreferences();
			
			// Docu in base class
			virtual void saveCurrentLayer(bool visible);
			
			/// Set view mode of the current layer
			void setViewMode(const ViewModes mode);

			/// Set draw mode of the current layer
			void setDrawMode(const LayerData::DrawModes mode);

			/// Set primitive mode of the current layer
			void setPrimitiveMode(const LayerData::PrimitiveModes mode);

			/// Set action of the current layer
			void setAction(const Actions action);

			/// Get draw mode of the current layer
			ViewModes getViewMode() const;
												
			/// Get draw mode of the current layer
			LayerData::DrawModes getDrawMode() const;

			/// Get draw mode of the current layer
			LayerData::PrimitiveModes getPrimitiveMode() const;
						
		public slots:
			
	    // Docu in base class
	    void activateLayer(Size layer_index);
	    // Docu in base class
	    void removeLayer(Size layer_index);
  	
  	protected slots:
  		
			/// Reacts on changed layer paramters
			void currentLayerParamtersChanged_();
  	
  	protected:
			
			// Docu in base class
	    bool finishAdding_();
			
  		// Reimplementation in order to update the OpenGL widget
  		virtual void update_(const char* caller_name);

			///whether the legend is shoen or not
			bool legend_shown_;

			//docu in base class
			virtual void updateLayer_(Size i);
      //docu in base class
			virtual void translateLeft_();
			//docu in base class
			virtual void translateRight_();
			//docu in base class
			virtual void translateForward_();
			//docu in base class
			virtual void translateBackward_();
  };
  
} //namespace
#endif
