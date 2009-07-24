// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H
#define OPENMS_VISUAL_SPECTRUM3DOPENGLCANVAS_H

#include <float.h>
#include <QtOpenGL/QGLWidget>

// OpenMS
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>

using std::cout;
using std::endl;
    
namespace OpenMS
{
	class Spectrum3DCanvas;
	
	/**
		@brief OpenGL Canvas for 3D-visualization of map data
		
		@note Do not use this class directly. Use Spectrum3DCanvas instead!
		
		@ingroup SpectrumWidgets
	*/
	
  class OPENMS_DLLAPI Spectrum3DOpenGLCanvas
  	: public QGLWidget
  {
    Q_OBJECT
		
			friend class Spectrum3DCanvas;

    public:
    				
			/// Container for axis ticks
			typedef std::vector<std::vector<double> > AxisTickVector;
			
	    /**
	     @brief Constructor
	     
	     @param parent The parent widget
	     @param canvas_3d The main 3d canvas
	    */
			Spectrum3DOpenGLCanvas(QWidget* parent, Spectrum3DCanvas& canvas_3d);
			/**
				@brief Destructor
			
				Destroys the OpenGLWidget and all associated data.
			*/		
	    virtual ~Spectrum3DOpenGLCanvas();
			///virtual function provided from QGLWidget	
			void initializeGL();
			/// virtual function provided from QGLWidget	
			void resizeGL(int w,int h);
	    /// virtual function provided from QGLWidget	
			void paintGL();
			/// Builds up a display list for the 3D view
			GLuint makeDataAsStick();
		  /// Builds up a display list for the axes
		  GLuint makeAxes();
			/// Builds up a display list for axis ticks
			GLuint makeAxesTicks();
			/// Builds up a display list for the birds-eye view
			GLuint makeDataAsTopView();
			/// Builds up a display list for the background
			GLuint makeGround();
			/// Builds up a display list for grid lines
			GLuint makeGridLines();
			/// Draws the axis texts (since Qt 4.3 these cannot be put into display lists anymore...)
			void drawAxesLegend();
	
      /** @name Reimplemented QT events */
      //@{
	    void mouseMoveEvent(QMouseEvent* e);
			void mouseReleaseEvent(QMouseEvent* e);
	    void mousePressEvent(QMouseEvent* e);
      void focusOutEvent(QFocusEvent* e);
      //@}
			
			/// computes the dataset supposed to be drawn when a section has been selected in zoom mode
			void computeSelection();
	
			/// updates the min and max values of the intensity
			void updateIntensityScale();
			
			/// calcualtes the zoom area , which is shown
			void dataToZoomArray(double x_1, double y_1, double x_2, double y_2);
			
			/// returns the BB-rt-coordinate :  value --> BB-coordinates
			double scaledRT(double rt);
			/// returns the rt-value : BB-coordinates  --> value
			double scaledInversRT(double mz);
		  /// returns the BB-mz-coordinate :  values --> BB-coordinates
			double scaledMZ(double mz);
			///  returns the mz-value : BB-coordinates  --> value
			double scaledInversMZ(double mz);
	    /// returns the BB-intensity -coordinate :  values --> BB-coordinates
			double scaledIntensity(Real intensity,Size layer_index);
	
			/// recalculates the dot gradient inerpolation values.
			void recalculateDotGradient_(Size layer);
			///calculate the ticks for the gridlines
			void calculateGridLines_();
		
			/// return width
			float width() const { return width_; }
			float height() const { return heigth_; }
			
	    /// return xRot_
	    int xRotation() const { return xrot_; }
	    /// return yRot_
			int yRotation() const { return yrot_; }
	    /// return zRot_
			int zRotation() const { return zrot_; }
	    /// normalize the angel
			void normalizeAngle(int* angle);
			//document me
			void setAngels(int xrot, int yrot, int zrot);
			//document me
			void resetTranslation();
			//document me
			void timeMessure();
			
			/// stores the original rotation and zoom factor (e.g. before changing into zoom mode)
			void storeRotationAndZoom();
			/// restores the original rotation and zoom factor (e.g. before changing into zoom mode)
			void restoreRotationAndZoom();
      
      /** @name Different OpenGL display lists */
      //@{
	    GLuint stickdata_;
			GLuint axes_;
			GLuint axes_ticks_;
			GLuint gridlines_;
			GLuint ground_;
			//@}
		
			/// reference to Spectrum3DCanvas
			Spectrum3DCanvas& canvas_3d_;
	  
			/// member x-variables for the rotation
	    int xrot_;
			/// member y-variables for the rotation
			int yrot_;
			/// member z-variables for the rotation
	    int zrot_;
	    
			/// member x-variable that stores the original angle during zoom mode 
	    int xrot_tmp_;
			/// member y-variable that stores the original angle during zoom mode 
			int yrot_tmp_;
			/// member z-variable that stores the original angle during zoom mode 
	    int zrot_tmp_;
	    
	    
	
			/// member variables fot the zoom-modus
	    QPoint mouse_move_end_, mouse_move_begin_;    
	
			///member variable for the x and y axis of the BB 
	    double corner_;
			/// member variable for the zoom mode
			double zoom_ ;
			/// member variable that stores original zoom factor during zoom mode
			double zoom_tmp_;
			
			/// member variable for the z- axis of the BB
			double near_;
	 		/// member variable for the z- axis of the BB
	 		double far_;
			/// the width of the viewport
			float width_;	
			/// the height of the viewport
			float heigth_;
			/// object which contains the min and max values of mz, rt and intensity
			DRange<3> overall_values_;
			///object wich contains the values of the current min and max intensity
			DRange<1> int_scale_;
			///member gridvectors which contains the data for the mz-axis-ticks
			AxisTickVector grid_mz_;
			///member gridvectors which contains the data for the rt-axis-ticks
			AxisTickVector grid_rt_;
			///member gridvectors which contains the data for the intensity-axis-ticks
			AxisTickVector grid_intensity_;
			/// x1 coordinate of the zoomselection
			double x_1_;
			/// x2 coordinate of the zoomselection
			double x_2_;
			/// y1 coordinate of the zoomselection
			double y_1_;
			/// y2 coordinate of the zoomselection
			double y_2_;
			/// x- translation
			double trans_x_;
			/// y_translation
			double trans_y_;
		  
		protected slots:
			/// Slot that reacts on action mode changes
			void actionModeChange();
	};

  struct vertexList
  {
    double x1, y1, z1, x2, y2, z2;
    QColor color1, color2;
  };
  
  class map2d
    : std::vector<double>
  {
    friend class Spectrum3DOpenGLCanvas;
    
    public:
      enum drawingMode { POINTS, BARPLOTS, PSEUDOGEL };
          
    private:        
      Size cols_;
      Size rows_;
      double mz_min_;
      double mz_max_;
      double mz_width_;
      double rt_min_;
      double rt_max_;
      double rt_width_;
      drawingMode drawing_mode_;
      Spectrum3DOpenGLCanvas* parent_;
      Size layer_index_;
      
      Size position_(const int col, const int row) const
      {
        int col_ = col < 0 ? 0 : col;
        col_ = col < (int) cols_ ? col : (int) cols_-1;
        
        int row_ = row < 0 ? 0 : row;
        row_ = row < (int) rows_ ? row : (int) rows_-1;
        
        return (Size) (col_ + cols_ * row_); 
      }

      Size position_(const double mz, const double rt) const
      {
        return position_(mzToIndex(mz), rtToIndex(rt));
      }

      Size position_(const Size col, const Size row) const
      {
        return position_((int) col, (int) row);
      }
                      
    public:
      typedef Spectrum3DCanvas::ExperimentType::ConstAreaIterator AreaIt;
            
      map2d(const Size cols, const Size rows, Spectrum3DOpenGLCanvas* parent, const Size layer_index)
      {
        parent_ = parent;
        layer_index_ = layer_index;
        drawing_mode_ = POINTS;
        cols_ = cols;
        rows_ = rows;
        
        resize(cols*rows, 0.0);
        
        mz_min_ = parent_->canvas_3d_.getVisibleArea().min_[0];
        mz_max_ = parent_->canvas_3d_.getVisibleArea().max_[0];
        mz_width_ = (mz_max_ - mz_min_) / (cols-1);
        
        rt_min_ = parent_->canvas_3d_.getVisibleArea().min_[1];
        rt_max_ = parent_->canvas_3d_.getVisibleArea().max_[1];
        rt_width_ = (rt_max_ - rt_min_) / (rows-1);
cout << "constructeur" << endl;
cout << "cols: " << cols_ << " - rows: " << rows_ << endl;
cout << "mz: " << mz_min_ << " / " << mz_max_ << " - width: " << mz_width_ << endl;
cout << "rt: " << rt_min_ << " / " << rt_max_ << " - width: " << rt_width_ << endl;
      }
      
      void setDrawingMode(const drawingMode mode)
      {
        drawing_mode_ = mode;
      }
      
      double indexToMz(const Size col) const
      {
if(col>cols_)
{
cout << "1 col: " << col << " / " << cols_ << endl;
return 100.0;
}        
        return (mz_min_ + (col<cols_ ? col : cols_-1) * mz_width_);
      }
      
      double indexToRt(const Size row) const
      {
if(row>rows_)
{
cout << "2 row: " << row << " / " << rows_ << endl;
return 100.0;
}      
        return (rt_min_ + (row<rows_ ? row : rows_-1) * rt_width_);
      }
      
      Size mzToIndex(const double mz) const
      {
        if(mz < mz_min_) return 0;
        if(mz > mz_max_) return (cols_-1);

Size col = (Size) floor(0.5 + (mz-mz_min_) / mz_width_);
if(col>cols_)
{
cout << "3 col: " << col << " / " << cols_ << endl;
return 1;
}
        return (Size) floor(0.5 + (mz-mz_min_) / mz_width_);
      }
      
      Size rtToIndex(const double rt) const
      {
        if(rt < rt_min_) return 0;
        if(rt > rt_max_) return (rows_-1);
        
Size row = (Size) floor(0.5 + (rt-rt_min_) / rt_width_);
if(row>rows_)
{
cout << "4 row: " << row << " / " << rows_ << endl;
return 1;
}        
        return (Size) floor(0.5 + (rt-rt_min_) / rt_width_);
      }
      
      void set(const double value, const Size col, const Size row)
      { 
if(col>cols_ || row>rows_)
{
cout << "5 col: " << col << " / " << cols_ << " - row: " << row << " / " << rows_ << endl;
return;
}  
      
        if(col<cols_ && row<rows_)
        {
          at(position_(col, row)) = value;
        }
      }
      
      double get(const Size col, const Size row) const
      { 
if(col>cols_ || row>rows_)
{
cout << "6 col: " << col << " / " << cols_ << " - row: " << row << " / " << rows_ << endl;
return 100.0;
}  
        return at(position_(col, row));
      }
    
      void set(const double value, const double mz, const double rt)
      { 
        if(mz>=mz_min_ && mz<=mz_max_ && rt>=rt_min_ && rt<=rt_max_)
        {
          at(position_(mz, rt)) = value;
        }
      }
      
      double get(const double mz, const double rt) const
      { 
        return at(position_(mz, rt));
      }
      
      void importData(AreaIt begin, AreaIt end)
      {
        time_t start_t,end_t;
        time (&start_t);
double mz1, mz2, rt1, rt2, int1, int2;
        
		    for(AreaIt it=begin; it!=end; ++it)
		    {
//cout << "7 mz: " << it->getMZ() << " - rt: " << it.getRT() << endl;		    
mz1 = mz1 > it->getMZ() ? it->getMZ() : mz1;
mz2 = mz2 < it->getMZ() ? it->getMZ() : mz2;
rt1 = rt1 > it.getRT() ?  it.getRT() : rt1;
rt2 = rt2 < it.getRT() ?  it.getRT() : rt2;
int1 = int1 > it->getIntensity() ? it->getIntensity() : int1;
int2 = int2 < it->getIntensity() ? it->getIntensity() : int2;

		      if(it->getMZ()>=mz_min_ && it->getMZ()<=mz_max_ && it.getRT()>=rt_min_ && it.getRT()<=rt_max_)
		      {		      
		        if(it->getIntensity() > get(it->getMZ(), it.getRT()))
		        {
		          set(it->getIntensity(), it->getMZ(), it.getRT());
		        }
		      }
		    }

cout << "ZZ mz: " << mz1 << " / " << mz2 << " rt: " << rt1 << " / " << rt2 << " int: " << int1 << " / " << int2 << endl;

		    time (&end_t);
		    double diff = difftime (end_t,start_t);
		    std::cout << "time:" << diff << std::endl;
      }

    // private:
	    double scaledRT(double rt)
	    {
		    double scaledrt = 2.0 * parent_->corner_ * (rt - rt_min_) / (rt_max_ -rt_min_);
		    return scaledrt;
	    }

	    double scaledMZ(double mz)
	    {
		    double scaledmz = 2.0 * parent_->corner_ * (mz - mz_min_) / (mz_max_ -mz_min_);
		    return scaledmz;
	    }
	
	    double scaledIntensity(Real intensity)
	    {
		    double scaledintensity = 0.0;
		    switch(parent_->canvas_3d_.getIntensityMode())
		    {
			    case SpectrumCanvas::IM_SNAP:
			    {
				    scaledintensity = 2.0 * parent_->corner_ * (intensity - parent_->int_scale_.min_[0]) / (parent_->int_scale_.max_[0] - parent_->int_scale_.min_[0]);
				    break;
				  }
				   
			    case SpectrumCanvas::IM_NONE:
			    {
			      double overall_intensity_min = parent_->canvas_3d_.getDataRange().min_[2];
			      double overall_intensity_max = parent_->canvas_3d_.getDataRange().max_[2];
			      
				    scaledintensity = 2.0 * parent_->corner_ * (intensity - overall_intensity_min) / (overall_intensity_max - overall_intensity_min);
				    break;
				  }
				    
			    case SpectrumCanvas::IM_PERCENTAGE: 
			    {
				    scaledintensity = 2.0 * parent_->corner_ * intensity / parent_->canvas_3d_.getMaxIntensity(layer_index_);	
				    break;
				  }
		    }
		    return scaledintensity;
	    }    
	      
    public:
      vertexList getVertex(const Size col, const Size row)
      {
if(col>cols_ || row>rows_)
{
cout << "8 col: " << col << " / " << cols_ << " - row: " << row << " / " << rows_ << endl;
}        
        
        vertexList vertex;
        const MultiGradient& gradient = parent_->canvas_3d_.getLayer(layer_index_).gradient;
        
        switch(drawing_mode_)
        {
          case BARPLOTS:
          {
				    switch (parent_->canvas_3d_.getIntensityMode())
			      {
					    case SpectrumCanvas::IM_PERCENTAGE:	
					    {
						    vertex.color2 = gradient.precalculatedColorAt(0);
						    break;
					    }
					    
					    case SpectrumCanvas::IM_NONE:
					    {
						    vertex.color2 = gradient.precalculatedColorAt(parent_->canvas_3d_.getDataRange().min_[2]);
						    break;
					    }
					    
					    case SpectrumCanvas::IM_SNAP:
					    {
						    vertex.color2 = gradient.precalculatedColorAt(parent_->int_scale_.min_[0]);
						    break;
						  }
			      }		            
		        vertex.x2 = -parent_->corner_ + (GLfloat) scaledMZ(indexToMz(col));
		        vertex.y2 = -parent_->corner_;
		        vertex.z2 = -parent_->near_ - 2 * parent_->corner_ - (GLfloat) scaledRT(indexToRt(row));	
		      }
		        			              
          case POINTS:
          {
				    switch (parent_->canvas_3d_.getIntensityMode())
			      {
				      case SpectrumCanvas::IM_PERCENTAGE:
				      {
				        vertex.color1 = gradient.precalculatedColorAt(get(col,row) * 100.0 / parent_->canvas_3d_.getMaxIntensity(layer_index_));
				        break;
				      }
				        
				      case SpectrumCanvas::IM_NONE:
				      case SpectrumCanvas::IM_SNAP:
				      {
				        vertex.color1 = gradient.precalculatedColorAt(get(col,row));
					      break;
					    }
			      }		       
		        vertex.x1 = -parent_->corner_ + (GLfloat) scaledMZ(indexToMz(col));
		        vertex.y1 = -parent_->corner_ + (GLfloat) scaledIntensity(get(col,row));		        
		        vertex.z1 = -parent_->near_ - 2 * parent_->corner_ - (GLfloat) scaledRT(indexToRt(row));			         
            break;
          }
            
          case PSEUDOGEL:
          {
				    switch (parent_->canvas_3d_.getIntensityMode())
			      {
				      case SpectrumCanvas::IM_PERCENTAGE:
				      {
				        vertex.color1 = gradient.precalculatedColorAt(get(col,row) * 100.0 / parent_->canvas_3d_.getMaxIntensity(layer_index_));
				        break;
				      }
				        
				      case SpectrumCanvas::IM_NONE:
				      case SpectrumCanvas::IM_SNAP:
				      {
				        vertex.color1 = gradient.precalculatedColorAt(get(col,row));
					      break;
					    }
			      }	
			      vertex.color2 = vertex.color1;
			                
            double width = scaledRT(indexToRt(row)) - scaledRT(indexToRt(row+1));
            
		        vertex.x1 = -parent_->corner_ + (GLfloat) scaledMZ(indexToMz(col));
		        vertex.y1 = -parent_->corner_ + (GLfloat) scaledIntensity(get(col,row));
		        vertex.z1 = -parent_->near_ - 2 * parent_->corner_ - (GLfloat) scaledRT(indexToRt(row)) - width / 2;	
		                  
		        vertex.x2 = vertex.x1;
		        vertex.y2 = vertex.y1;
		        vertex.z2 = -parent_->near_ - 2 * parent_->corner_ - (GLfloat) scaledRT(indexToRt(row)) + width / 2;	          
            break;
          }
        }
        return vertex;
      }
      
      vertexList getVertex(const int col, const int row)
      {
if(col>(int)cols_ || row>(int)rows_ || col<0 || row<0)
{
cout << "10 col: " << col << " / " << cols_ << " - row: " << row << " / " << rows_ << endl;
}   
        Size col_ = col<0 ? 0 : (Size) col;
        Size row_ = row<0 ? 0 : (Size) row;
        return getVertex(col_, row_);               
      }
            
    private:
      vertexList vertexToVector_(const vertexList v1, const vertexList v2)
      {
        vertexList vertex;
        vertex.x1 = v2.x1 - v1.x1;
        vertex.y1 = v2.y1 - v1.y1;
        vertex.z1 = v2.z1 - v1.z1;
        return vertex;
      }
      
      vertexList crossProduct_(const vertexList v1, const vertexList v2)
      {
        vertexList vertex;
        vertex.x1 = v1.y1 * v2.z1 + v1.z1 * v2.y1;
        vertex.y1 = v1.z1 * v2.x1 + v1.x1 * v2.z1;
        vertex.z1 = v1.x1 * v2.y1 + v1.y1 * v2.x1;
        return vertex;
      }     
      
      vertexList addVectors_(const vertexList v1, const vertexList v2)
      {
        vertexList vertex;
        vertex.x1 = v1.x1 + v2.x1;
        vertex.y1 = v1.y1 + v2.y1;
        vertex.z1 = v1.z1 + v2.z1;
        return vertex;
      }     
    
    public:
      vertexList getNormals(const Size col, const Size row)
      {
if(col>cols_ || row>rows_)
{
cout << "9 col: " << col << " / " << cols_ << " - row: " << row << " / " << rows_ << endl;
}         
        vertexList vertex;
        vertex.x1 = 0.0;
        vertex.y1 = 0.0;
        vertex.z1 = 0.0;
        
        vertexList mainVertex = getVertex(col,row);
      
        switch(drawing_mode_)
        {
          case BARPLOTS:  
          case POINTS:
          {
            int cc[] = {-1,  0, +1, +1, +1,  0, -1, -1, -1};
            int rr[] = {-1, -1, -1,  0, +1, +1, +1,  0, -1};
            
            for(Size index=0; index<8; ++index)
            {
              vertex = addVectors_(vertex, crossProduct_(
                vertexToVector_(mainVertex, getVertex(col+cc[index],row+rr[index])), 
                vertexToVector_(mainVertex, getVertex(col+cc[index+1],row+rr[index+1]))));
            }
            break;
          }
             
          case PSEUDOGEL:
          {
            vertexList v;
            v.x1 = 0.0;
            v.y1 = 0.0;
            v.z1 = 1.0;
            
            vertex = addVectors_(
              crossProduct_(
                vertexToVector_(mainVertex, getVertex(col+1,row)), 
                vertexToVector_(mainVertex, v)),
              crossProduct_(
                vertexToVector_(mainVertex, v), 
                vertexToVector_(mainVertex, getVertex(col-1,row))));
            break;
          }
        }
        
        // normalize
        double length = sqrt(vertex.x1 * vertex.x1 + vertex.y1 * vertex.y1 + vertex.z1 * vertex.z1);
        vertex.x1 /= length;
        vertex.y1 /= length;
        vertex.z1 /= length;
        
        return vertex;
      }
  };	
}
#endif
