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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/MATH/MISC/Struct3d.h>

#include <QtGui/QMouseEvent>
#include <QtGui/QKeyEvent>

namespace OpenMS
{
	  
	Spectrum3DOpenGLCanvas::Spectrum3DOpenGLCanvas(QWidget *parent,Spectrum3DCanvas& canvas_3d)
	  : QGLWidget(parent),
	    canvas_3d_(canvas_3d)
	{
		canvas_3d.rubber_band_.setParent(this);
		
		//Set focus policy and mouse tracking in order to get keyboard events
		setMouseTracking(TRUE);
		setFocusPolicy(Qt::StrongFocus);
	  
	  corner_ = 100.0;  
	  near_ = 0.0;  
	  far_ = 600.0;
	  zoom_ = 1.5;
	  xrot_ = 220;
	  yrot_ = 220;
	  zrot_ = 0;
	  trans_x_ = 0.0;
	  trans_y_ = 0.0;
		view_mode_ = Spectrum3DCanvas::VM_3D;
	}
	  
	Spectrum3DOpenGLCanvas::~Spectrum3DOpenGLCanvas()
	{
	}
	
	void Spectrum3DOpenGLCanvas::calculateGridLines_()
	{
	  switch(canvas_3d_.intensity_mode_)
	  {
		  case SpectrumCanvas::IM_SNAP:
		    updateIntensityScale();
		    AxisTickCalculator::calcGridLines(int_scale_.min_[0],int_scale_.max_[0],3,grid_intensity_,7,5); 
		    break;
		  case SpectrumCanvas::IM_NONE:
		    AxisTickCalculator::calcGridLines(canvas_3d_.overall_data_range_.min_[2],canvas_3d_.overall_data_range_.max_[2],3,grid_intensity_,7,5); 
		    break;
		  case SpectrumCanvas::IM_PERCENTAGE:
		    AxisTickCalculator::calcGridLines(0.0,100.0,3,grid_intensity_,7,5); 
		    break;
	  }
	  
	  AxisTickCalculator::calcGridLines(canvas_3d_.visible_area_.min_[1],canvas_3d_.visible_area_.max_[1],3,grid_rt_,7,5);
	  AxisTickCalculator::calcGridLines(canvas_3d_.visible_area_.min_[0],canvas_3d_.visible_area_.max_[0],3,grid_mz_,7,5);
	}
	
	void Spectrum3DOpenGLCanvas::initializeGL()
	{	
	  QColor color(canvas_3d_.param_.getValue("background_color").toQString());
	  qglClearColor(color);
	  glEnable(GL_DEPTH_TEST);
	  glEnable(GL_BLEND);
	  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	  calculateGridLines_();
		
		//abort if no layers are displayed
	  if(canvas_3d_.getLayerCount()==0) return;

	  if(canvas_3d_.action_mode_==SpectrumCanvas::AM_ZOOM)
	  {
      if(!canvas_3d_.rubber_band_.isVisible())
      {
        axes_ = makeAxes();
        if(canvas_3d_.show_grid_)
        {
          gridlines_ = makeGridLines();
        }
        xrot_ = 90*16;
        yrot_ = 0;
        zrot_ = 0;  
        zoom_ = 1.25;
        stickdata_ = makeData();
        axes_ticks_ = makeAxesTicks();   
        drawAxesLegend();
      }						
		}
		else if(canvas_3d_.action_mode_==SpectrumCanvas::AM_TRANSLATE)
	  {
			if(canvas_3d_.show_grid_)
      {
        gridlines_ = makeGridLines();
      }   
      axes_ = makeAxes();   
      ground_ = makeGround();
      x_1_ = 0.0;
      y_1_ = 0.0;
      x_2_ = 0.0;
      y_2_ = 0.0;
      stickdata_ = makeData();
      axes_ticks_ = makeAxesTicks();
      drawAxesLegend();
	  }
	}

	void Spectrum3DOpenGLCanvas::paintGL()
	{	
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glLoadIdentity();
	
		glTranslated(0.0, 0.0,-3.0*corner_);
		glRotated(xrot_ / 16.0, 1.0, 0.0, 0.0);
		glRotated(yrot_ / 16.0, 0.0, 1.0, 0.0);
		glRotated(zrot_/16.0, 0.0, 0.0, 1.0);
		glTranslated(trans_x_, trans_y_,3.0*corner_);
		
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		//abort if no layers are displayed
	  if(canvas_3d_.getLayerCount()==0) return;
	  
		glCallList(ground_);
		if(canvas_3d_.show_grid_)
		{
			glCallList(gridlines_);
		}
		glCallList(axes_);
		glCallList(axes_ticks_);
		drawAxesLegend();
		if (canvas_3d_.action_mode_==SpectrumCanvas::AM_ZOOM || canvas_3d_.action_mode_==SpectrumCanvas::AM_TRANSLATE)
		{
			glCallList(stickdata_);
		}
	}
	
	void Spectrum3DOpenGLCanvas::resizeGL(int w,int h)
	{
    width_ = (float) w;
    heigth_ = (float) h;
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);   
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-corner_*zoom_, corner_*zoom_,  -corner_*zoom_, corner_*zoom_ , near_, far_);
    glMatrixMode(GL_MODELVIEW);
	}	
	
	void Spectrum3DOpenGLCanvas::setAngels(int xrot, int yrot, int zrot)
	{
		xrot_=xrot;
		yrot_=yrot;
		zrot_=zrot;
	}
	
	void Spectrum3DOpenGLCanvas::resetTranslation()
	{
		trans_x_ = 0.0;
		trans_y_ = 0.0;
	}
	
	void Spectrum3DOpenGLCanvas::storeRotationAndZoom()
	{
		xrot_tmp_ = xrot_;
		yrot_tmp_ = yrot_;
		zrot_tmp_ = zrot_;
		zoom_tmp_ = zoom_;
	}
	
	void Spectrum3DOpenGLCanvas::restoreRotationAndZoom()
	{
		xrot_ = xrot_tmp_;
		yrot_ = yrot_tmp_;
		zrot_ = zrot_tmp_;
		zoom_ = zoom_tmp_;
	}
	
	void Spectrum3DOpenGLCanvas::drawAxesLegend()
	{
		QFont font("Typewriter");
		font.setPixelSize(10);

		QString text;
		qglColor(Qt::black);
		
		//RT axis legend
 		if(canvas_3d_.legend_shown_)
		{
			font.setPixelSize(12);
			
			static QString mz_label = (String(Peak2D::shortDimensionName(Peak2D::MZ))+" ["+String(Peak2D::shortDimensionUnit(Peak2D::MZ))+"]").toQString();
	 		renderText(0.0,  -corner_-20.0,  -near_-2*corner_+20.0, mz_label, font);
	 		
	 		static QString rt_label = (String(Peak2D::shortDimensionName(Peak2D::RT))+" ["+String(Peak2D::shortDimensionUnit(Peak2D::RT))+"]").toQString();
			renderText(-corner_-20.0, -corner_-20.0, -near_-3*corner_, rt_label, font);
			
			font.setPixelSize(10);
		}
		
		//RT numbers
		if(grid_rt_.size()>0)
		{
			for (Size i = 0;i<grid_rt_[0].size();i++)
			{
				text = QString::number(grid_rt_[0][i]);
				renderText(-corner_-15.0, -corner_-5.0, -near_-2*corner_-scaledRT(grid_rt_[0][i]), text, font);
			}
		}
		if(zoom_<3.0 && grid_rt_.size()>=2)
		{
			for (Size i = 0;i<grid_rt_[1].size();i++)
			{
				text = QString::number(grid_rt_[1][i]);
				renderText(-corner_-15.0,  -corner_-5.0, -near_-2*corner_-scaledRT(grid_rt_[1][i]), text, font);
			}
		}
		if(zoom_<2.0 && grid_rt_.size()>=3)
		{
			for (Size i = 0;i<grid_rt_[2].size();i++)
			{
				text = QString::number(grid_rt_[2][i]);
				renderText(-corner_-15.0, -corner_-5.0, -near_-2*corner_-scaledRT(grid_rt_[2][i]), text, font);
			}
		}
		
		//m/z numbers
		if(grid_mz_.size()>0)
		{
			for (Size i = 0;i<grid_mz_[0].size();i++)
			{
				text = QString::number(grid_mz_[0][i]);
				renderText(-corner_-text.length()+scaledMZ(grid_mz_[0][i]), -corner_-5.0, -near_-2*corner_+15.0, text, font);
			}
		}
		if(zoom_<3.0 && grid_mz_.size()>=2)
		{
			for (Size i = 0;i<grid_mz_[1].size();i++)
			{
				text = QString::number(grid_mz_[1][i]);
				renderText(-corner_-text.length()+scaledMZ(grid_mz_[1][i]), -corner_-5.0, -near_-2*corner_+15.0, text, font);
			}
		}
		if(zoom_<2.0 && grid_mz_.size()>=3)
		{
			for (Size i = 0;i<grid_mz_[2].size();i++)
			{
				text = QString::number(grid_mz_[2][i]);
				renderText(-corner_-text.length()+scaledMZ(grid_mz_[2][i]), -corner_-5.0, -near_-2*corner_+15.0, text, font);
			}
		}
		
		//draw intensity legend if not in zoom mode
		if(canvas_3d_.action_mode_ != SpectrumCanvas::AM_ZOOM)
		{
			switch (canvas_3d_.intensity_mode_)
				{
				case SpectrumCanvas::IM_PERCENTAGE:
			 		if(canvas_3d_.legend_shown_)
					{
						font.setPixelSize(12);
						renderText(-corner_-20.0, corner_+10.0, -near_-2*corner_+20.0, "intensity %", font);
						font.setPixelSize(10);
					}
					
					for (Size i = 0;i<grid_intensity_[0].size();i++)
					{ 
						text = QString::number(grid_intensity_[0][i]);
						renderText(-corner_-text.length()-width_/200.0-5.0, -corner_+(2.0*grid_intensity_[0][i]), -near_-2*corner_, text, font);
					}
					break;
				case SpectrumCanvas::IM_NONE:
				case SpectrumCanvas::IM_SNAP:
					int expo = 0;
					if(grid_intensity_.size()>=1)
					{
						expo = (int)ceil(log10(grid_intensity_[0][0]));
					}
					if(grid_intensity_.size()>=2)
					{
						if(expo>=ceil(log10(grid_intensity_[1][0])))
						{
							expo = (int)ceil(log10(grid_intensity_[1][0]));
						}
					}
					if(grid_intensity_.size()>=3)
					{
						if(expo>=ceil(log10(grid_intensity_[2][0])))
						{
							expo =(int) ceil(log10(grid_intensity_[2][0]));
						}	
					}
					
			 		if(canvas_3d_.legend_shown_)
					{
						font.setPixelSize(12);
						text = QString("intensity e+%1").arg((double)expo,0,'f',1);
						renderText(-corner_-20.0, corner_+10.0, -near_-2*corner_+20.0, text, font);
						font.setPixelSize(10);
					}
					
					
					if(zoom_<3.0 && grid_intensity_.size()>=2)
					{
						for (Size i = 0;i<grid_intensity_[0].size();i++)
						{ 
							double intensity = (double)grid_intensity_[0][i]/pow(10.0,expo);
							text = QString("%1").arg(intensity,0,'f',1);
							renderText(-corner_-text.length()-width_/200.0-5.0, -corner_+scaledIntensity(grid_intensity_[0][i],canvas_3d_.current_layer_), -near_-2*corner_, text, font);
						}
						for (Size i = 0;i<grid_intensity_[1].size();i++)
						{
							double intensity = (double)grid_intensity_[1][i]/pow(10.0,expo);
							text = QString("%1").arg(intensity,0,'f',1);
							renderText(-corner_-text.length()-width_/200.0-5.0, -corner_+scaledIntensity(grid_intensity_[1][i],canvas_3d_.current_layer_), -near_-2*corner_, text, font);
						}
					}
					if(width_>800 && heigth_>600&& zoom_<2.0 && grid_intensity_.size()>=3)
					{
						for (Size i = 0;i<grid_intensity_[2].size();i++)
						{
							double intensity = (double)grid_intensity_[2][i]/pow(10.0,expo);
							text = QString("%1").arg(intensity,0,'f',1);
							renderText(-corner_-text.length()-width_/200.0-5.0, -corner_+scaledIntensity(grid_intensity_[2][i],canvas_3d_.current_layer_), -near_-2*corner_, text, font);
						}
					}
				break;
			}
		}
	}
	
	GLuint Spectrum3DOpenGLCanvas::makeGround()
	{
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		glBegin(GL_QUADS);
	 	QColor color(canvas_3d_.param_.getValue("background_color").toQString());
		qglColor(color);
		glVertex3d(-corner_, -corner_-2.0, -near_-2*corner_);
		glVertex3d(-corner_, -corner_-2.0, -far_+2*corner_);
		glVertex3d(corner_, -corner_-2.0, -far_+2*corner_);
		glVertex3d(corner_, -corner_-2.0, -near_-2*corner_);
		glEnd();
		glEndList();
		return list;
	}
	
	GLuint Spectrum3DOpenGLCanvas::makeAxes()
	{
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		glLineWidth(3.0);
		glShadeModel(GL_FLAT);
		glBegin(GL_LINES);
		qglColor(Qt::black);
		//x_achse
		glVertex3d(-corner_, -corner_, -near_-2*corner_);
		glVertex3d( corner_, -corner_, -near_-2*corner_);
		//z-achse
		glVertex3d(-corner_, -corner_, -near_-2*corner_);
		glVertex3d( -corner_, -corner_, -far_+2*corner_);
		//y-achse
		glVertex3d(-corner_, -corner_, -near_-2*corner_);
		glVertex3d( -corner_, corner_, -near_-2*corner_);
		glEnd();
		glEndList();
		return list;
	}
	
/*	GLuint Spectrum3DOpenGLCanvas::makeDataAsTopView()
	{
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		glPointSize(3.0);
	
		for (Size i =0;i<canvas_3d_.getLayerCount();++i)
		{	
			const LayerData& layer = canvas_3d_.getLayer(i);
			if(layer.visible)
			{	
				if((Int)layer.param.getValue("dot:shade_mode"))
				{
					glShadeModel(GL_SMOOTH); 
				}
				else
				{
					glShadeModel(GL_FLAT); 
				}
				
				for (Spectrum3DCanvas::ExperimentType::ConstAreaIterator it = layer.peaks.areaBeginConst(canvas_3d_.visible_area_.min_[1],canvas_3d_.visible_area_.max_[1],canvas_3d_.visible_area_.min_[0],canvas_3d_.visible_area_.max_[0]); 
						 it != layer.peaks.areaEndConst(); 
						 ++it)
				{
					PeakIndex pi = it.getPeakIndex();
					if (layer.filters.passes(layer.peaks[pi.spectrum],pi.peak))
					{
						glBegin(GL_POINTS);
						double intensity = 0;
						switch (canvas_3d_.intensity_mode_)
						{
							case SpectrumCanvas::IM_NONE:
								qglColor(layer.gradient.precalculatedColorAt(it->getIntensity()));
								break;
							case SpectrumCanvas::IM_PERCENTAGE:	
								intensity = it->getIntensity() * 100.0 /canvas_3d_.getMaxIntensity(i);
								qglColor(layer.gradient.precalculatedColorAt(intensity ));
								break;
							case SpectrumCanvas::IM_SNAP:
								qglColor(layer.gradient.precalculatedColorAt(it->getIntensity()));
								break;
						}
						glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()),
											 -corner_,
											 -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						glEnd();		
					}
				}
			}
		}
		glEndList();
		return list; 
	}*/

  void Spectrum3DOpenGLCanvas::redraw()
  {
cout << "10" << endl;
	  //update the content
		canvas_3d_.update_buffer_ = true;
		canvas_3d_.update_(__PRETTY_FUNCTION__);  
  }
  
	void Spectrum3DOpenGLCanvas::setViewMode(const Spectrum3DCanvas::ViewModes mode)
	{
	  view_mode_ = mode;
	  redraw();
	}
	
	void Spectrum3DOpenGLCanvas::setMappingMode(const MappingThread::MappingModes mode)
	{
	  canvas_3d_.getCurrentLayer().setMappingMode(mode);
	  redraw();
	}
	
	void Spectrum3DOpenGLCanvas::setPrimitiveMode(const LayerData::PrimitiveModes mode)
	{
	  canvas_3d_.getCurrentLayer().setPrimitiveMode(mode);
	  redraw();
	}
	
	void Spectrum3DOpenGLCanvas::setAction(const Spectrum3DCanvas::Actions action)
	{
    switch(action)
    {                  
      case Spectrum3DCanvas::A_CAMERA_RESET :	
	      corner_ = 100.0;  
	      near_ = 0.0;  
	      far_ = 600.0;
	      zoom_ = 1.5;
	      xrot_ = 220;
	      yrot_ = 220;
	      zrot_ = 0;
	      trans_x_ = 0.0;
	      trans_y_ = 0.0;
	      break;
      case Spectrum3DCanvas::A_CAMERA_MOVEUP :
        corner_ += 20.0;
        break;
      case Spectrum3DCanvas::A_CAMERA_MOVEDOWN :
        corner_ -= 20.0;
        break;
      case Spectrum3DCanvas::A_CAMERA_MOVELEFT :
        near_ += 10.0;
        break;
      case Spectrum3DCanvas::A_CAMERA_MOVERIGHT :
        near_ -= 10.0;
        break;
      case Spectrum3DCanvas::A_CAMERA_ZOOMIN :
	      zoom_ += 0.2;
	      break;
      case Spectrum3DCanvas::A_CAMERA_ZOOMOUT :
      	zoom_ -= 0.2;
      	break;
      case Spectrum3DCanvas::A_DATA_RESET :
        far_ += 20.0;
        break;
      case Spectrum3DCanvas::A_DATA_MOVEUP :
        far_ -= 20.0;
        break;
      case Spectrum3DCanvas::A_DATA_MOVEDOWN :
        xrot_ += 20;
        break;      
      case Spectrum3DCanvas::A_DATA_MOVELEFT :
        yrot_ += 20;
        break;      
      case Spectrum3DCanvas::A_DATA_MOVERIGHT :
        zrot_ += 20;
        break;      
      case Spectrum3DCanvas::A_DATA_ZOOMIN :
        trans_x_ += 20.0;
        break;      
      case Spectrum3DCanvas::A_DATA_ZOOMOUT :
        trans_y_ += 20.0;
        break;      
    }
	  redraw();
	}
	
	Spectrum3DCanvas::ViewModes Spectrum3DOpenGLCanvas::getViewMode() const
	{
	  return view_mode_;	 	  
	}

	MappingThread::MappingModes Spectrum3DOpenGLCanvas::getMappingMode() const
	{
	  return canvas_3d_.getCurrentLayer().getMappingMode();
	}
	
	LayerData::PrimitiveModes Spectrum3DOpenGLCanvas::getPrimitiveMode() const
	{
	  return canvas_3d_.getCurrentLayer().getPrimitiveMode();
	}
			
	GLuint Spectrum3DOpenGLCanvas::makeData()
	{
cout << "Spectrum3DOpenGLCanvas::makeData()" << endl;
cout << "corner_: " << corner_ << " near_: " << near_ << " far_: " << far_ << " zoom_: " << zoom_ << endl;
cout << "xrot_: " << xrot_ << " yrot_: " << yrot_ << " zrot_: " << zrot_ << " trans_x_: " << trans_x_ << " trans_y_: " << trans_y_ << endl;

	  // init gl
		GLuint list = glGenLists(1);
	  glNewList(list, GL_COMPILE);

    if(Spectrum3DCanvas::VM_2D == getViewMode())
    {
      setAngels(220,220,0);
    }
				        
    for(Size iLayer=0; iLayer<canvas_3d_.getLayerCount(); ++iLayer)
    {
	    LayerData& layer = canvas_3d_.getLayer(iLayer);
	    if(layer.visible)
	    {
        // set drawing width
        glPointSize(layer.param.getValue("dot:line_width"));
        glLineWidth(layer.param.getValue("dot:line_width"));
        
        // set shade mode
        if((Int) layer.param.getValue("dot:shade_mode"))
        {
          glShadeModel(GL_SMOOTH);
          glEnable(GL_POINT_SMOOTH);
          glEnable(GL_LINE_SMOOTH);
        }
        else
        {
          glShadeModel(GL_FLAT);
          glDisable(GL_POINT_SMOOTH);
          glDisable(GL_LINE_SMOOTH);
        }			    

        recalculateDotGradient_(iLayer);
            
	      if(MappingThread::MM_NONE == layer.getMappingMode())
	      {
	        double currentRT = 0.0;
          glBegin(GL_POINTS);
                  	        
          for(Spectrum3DCanvas::ExperimentType::ConstAreaIterator it = layer.peaks.areaBeginConst(
		          canvas_3d_.visible_area_.min_[1],
		          canvas_3d_.visible_area_.max_[1],
		          canvas_3d_.visible_area_.min_[0],
		          canvas_3d_.visible_area_.max_[0]); 
				    it != layer.peaks.areaEndConst(); 
				    ++it)
		      {
			      PeakIndex pi = it.getPeakIndex();
			      if (layer.filters.passes(layer.peaks[pi.spectrum],pi.peak))
			      {
			        Struct3d baseVertex(- corner_ + (GLfloat) scaledMZ(it->getMZ()),
			                            - near_ - 2 * corner_ - (GLfloat) scaledRT(it.getRT()),
											            - corner_);
			        Struct3d topVertex(- corner_ + (GLfloat) scaledMZ(it->getMZ()),
			                            - near_ - 2 * corner_ - (GLfloat) scaledRT(it.getRT()),
											            - corner_ + (GLfloat) scaledIntensity(it->getIntensity(), iLayer));
									         
							QColor baseColor, topColor;
				      switch (canvas_3d_.intensity_mode_)
				      {
					      case SpectrumCanvas::IM_NONE:
						      baseColor = QColor( layer.gradient.precalculatedColorAt(canvas_3d_.overall_data_range_.min_[2]));
						      topColor = QColor( layer.gradient.precalculatedColorAt(it->getIntensity()));
						      break;				      
					      case SpectrumCanvas::IM_PERCENTAGE:	
						      baseColor = QColor( layer.gradient.precalculatedColorAt(0));
						      topColor = QColor( layer.gradient.precalculatedColorAt(it->getIntensity() * 100.0 /canvas_3d_.getMaxIntensity(iLayer)) );
						      break;
					      case SpectrumCanvas::IM_SNAP:
						      baseColor = QColor(layer.gradient.precalculatedColorAt(int_scale_.min_[0]));
						      topColor = QColor(layer.gradient.precalculatedColorAt(it->getIntensity()));
						      break;
				      }
				      
				      switch(getPrimitiveMode())
              {                  
                case LayerData::PM_LINES :
                  if(it.getRT() != currentRT)
                  {
                    glEnd();
                    glBegin(GL_LINES);
                    currentRT = it.getRT();
                  }
                  qglColor(baseColor);
                  glVertex3d(baseVertex.mz, baseVertex.intensity, baseVertex.rt);
                  qglColor(topColor);
                  glVertex3d(topVertex.mz, topVertex.intensity, topVertex.rt);
                  break;
                  
                case LayerData::PM_LINESTRIP :
                  if(it.getRT() != currentRT)
                  {
                    glEnd();
                    glBegin(GL_LINES);
                    currentRT = it.getRT();
                  }
                  qglColor(topColor);
                  glVertex3d(topVertex.mz, topVertex.intensity, topVertex.rt);
                  break;
                  
                case LayerData::PM_TRIANGLES :
                  if(it.getRT() != currentRT)
                  {
                    glEnd();
                    glBegin(GL_TRIANGLE_STRIP);
                    currentRT = it.getRT();
                  }
                  qglColor(baseColor);
                  glVertex3d(baseVertex.mz, baseVertex.intensity, baseVertex.rt);
                  qglColor(topColor);
                  glVertex3d(topVertex.mz, topVertex.intensity, topVertex.rt);
                  break;

                case LayerData::PM_POINTS :
                default :
                  if(it.getRT() != currentRT)
                  {
                    glEnd();
                    glBegin(GL_POINTS);
                    currentRT = it.getRT();
                  }
                  qglColor(topColor);
                  glVertex3d(topVertex.mz, topVertex.intensity, topVertex.rt);
                  break;             
              }
				      glEnd();
			      }
		      }
	      }
	      else
	      {
          if(layer.getMappingThread()->isValide())
          {
            switch(getPrimitiveMode())
            {
              case LayerData::PM_POINTS :
                glBegin(GL_POINTS);
                break;
              case LayerData::PM_LINES :
                glBegin(GL_LINES);
                break;
              case LayerData::PM_LINESTRIP :
                glBegin(GL_LINE_STRIP);
                break;
              case LayerData::PM_TRIANGLES :
                glBegin(GL_TRIANGLES);
                break;
              case LayerData::PM_TRIANGLESTRIP :
                glBegin(GL_TRIANGLE_STRIP);
                break;
              case LayerData::PM_QUADS :
                glBegin(GL_QUADS);
                break;
              case LayerData::PM_QUADSTRIP :
                glBegin(GL_QUAD_STRIP);
                break;
              case LayerData::PM_POLYGON :
                glBegin(GL_POLYGON);
                break;
	            default :
                glBegin(GL_POINTS);	            
	              break;                
            }
		          
            Vector3d vertexList = layer.getMappingThread()->getVertex();
            Vector3d normalsList = layer.getMappingThread()->getNormals();
            Iterator3d vertexIt = vertexList.begin();
            Iterator3d normalIt = normalsList.begin();
            for(; vertexIt!=vertexList.end(); ++vertexIt)
            {
			        switch (canvas_3d_.intensity_mode_)
			        {
				        case SpectrumCanvas::IM_PERCENTAGE :	
				          qglColor( layer.gradient.precalculatedColorAt(vertexIt->intensity * 100.0 / canvas_3d_.getMaxIntensity(iLayer)) );
				          break;
				        case SpectrumCanvas::IM_NONE :	
				          qglColor( layer.gradient.precalculatedColorAt(vertexIt->intensity) );
				          break;
				        case SpectrumCanvas::IM_SNAP :	
				          qglColor( layer.gradient.precalculatedColorAt(vertexIt->intensity) );
				          break;
				      }
				      if(Spectrum3DCanvas::VM_2D == getViewMode())
				      {
			          glVertex3d(- corner_ + (GLfloat) scaledMZ(vertexIt->mz),
								           - corner_,
								           - near_ - 2 * corner_ - (GLfloat) scaledRT(vertexIt->rt));
						  }
						  else
						  {
			          glVertex3d(- corner_ + (GLfloat) scaledMZ(vertexIt->mz),
								           - corner_ + (GLfloat) scaledIntensity(vertexIt->intensity, iLayer),
								           - near_ - 2 * corner_ - (GLfloat) scaledRT(vertexIt->rt));
								if(!normalsList.empty())
								{
								  // add normals
								  // glNormal3d(normal.x1, normal.y1, normal.z1);    
								}
						  }
            }
	          glEnd();
          }
          else
          {
  cout << "14" << endl;
					  if( !layer.getMappingThread()->isRunning() )
					  {
		          Size rows = layer.peaks.RTEnd(canvas_3d_.visible_area_.max_[1]) - layer.peaks.RTBegin(canvas_3d_.visible_area_.min_[1]);
	  cout << "15" << endl;            
		          Size cols = (Size) ceil((canvas_3d_.visible_area_.max_[0] - canvas_3d_.visible_area_.min_[0]) / 0.5);
		          if(cols > 300) 
		            cols = 300;
	  cout << "16" << endl;
		          layer.getMappingThread()->setDataSize(cols, rows);
		          layer.getMappingThread()->setRange(
		              canvas_3d_.visible_area_.min_[0],
		              canvas_3d_.visible_area_.max_[0],
		              canvas_3d_.visible_area_.min_[1],
		              canvas_3d_.visible_area_.max_[1]);
		          connect(layer.getMappingThread(), SIGNAL(finish()), this, SLOT(redraw()));
	  cout << "17" << endl;
		          layer.getMappingThread()->start();
	  cout << "18" << endl;
					  }
	  cout << "19" << endl;
          }
        }
      }		
    }

		glEndList();
		return list; 
	}   
			    
		  /*	  
	      recalculateDotGradient_(iLayer);
	       
	      if((Int)layer.param.getValue("dot:shade_mode"))
			    glShadeModel(GL_SMOOTH); 
		    else
			    glShadeModel(GL_FLAT); 
	
	      if(pSimple)
	      {			  			      
	        map_->init(cols, rows, this, iLayer);
	        map_->importData(
	          layer.peaks.areaBeginConst(
	            canvas_3d_.visible_area_.min_[1],
	            canvas_3d_.visible_area_.max_[1],
	            canvas_3d_.visible_area_.min_[0],
	            canvas_3d_.visible_area_.max_[0]),
				    layer.peaks.areaEndConst());
					
			    if(getMappingMode() == Spectrum3DCanvas::DM_POINTS)
		      {			    
		        glBegin(GL_POINTS);
		        glPointSize(10.0);
	        
		        for(Size ii=0; ii<cols; ++ii)
		        {
		          for(Size jj=0; jj<rows; ++jj)
		          {				    
		            vertexList vertex = map_->getVertex(ii,jj);  
		            qglColor(vertex.color1);
		            //glVertex3d(vertex.x1, vertex.y1, vertex.z1);
		          }
		        }
		        glEnd();				    						    
		      }
		      
			    if(getMappingMode() == Spectrum3DCanvas::DM_PEAKS)
		      {				    
		        glBegin(GL_LINES);
		        
		        for(Size ii=0; ii<cols; ++ii)
		        {
		          for(Size jj=0; jj<rows; ++jj)
		          {
		            vertexList vertex = map_->getVertex(ii,jj);  
		            qglColor(vertex.color1);
		            glVertex3d(vertex.x1, vertex.y1-1, vertex.z1);
		            qglColor(vertex.color2);
		            glVertex3d(vertex.x2, vertex.y2+1, vertex.z2);
		          }
		        }
		        glEnd();
		      }
		      
		      if(getMappingMode() == Spectrum3DCanvas::DM_LINES || getMappingMode() == Spectrum3DCanvas::DM_MAP)
		      {
		        for(Size jj=0; jj<rows-1; ++jj)
		        {
		          if(getMappingMode() == Spectrum3DCanvas::DM_MAP)
		          {
    			      glBegin(GL_TRIANGLE_STRIP);

                vertexList vertex = map_->getVertex(Size(0), jj);
                qglColor(vertex.color1);
                glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                //vertexList normal = map_->getNormals(0, jj);
                //glNormal3d(normal.x1, normal.y1, normal.z1);
                
                vertex = map_->getVertex(Size(0), jj+1);
                qglColor(vertex.color1);
                glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                //normal = map_->getNormals(0, jj+1);
                //glNormal3d(normal.x1, normal.y1, normal.z1);              	        
		              
		            for(Size ii=0; ii<cols-1; ++ii)
		            {
                  vertex = map_->getVertex(ii+1, jj);
                  qglColor(vertex.color1);
                  glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                  //vertexList normal = map_->getNormals(ii+1, jj);
                  //glNormal3d(normal.x1, normal.y1, normal.z1);                
                  
                  vertex = map_->getVertex(ii+1, jj+1);
                  qglColor(vertex.color1);
                  glVertex3d(vertex.x1, vertex.y1, vertex.z1);					
                  //normal = map_->getNormals(ii+1, jj+1);
                  //glNormal3d(normal.x1, normal.y1, normal.z1);   
                }
		            glEnd();
		          }

		          {
				        glBegin(GL_LINE_STRIP);

                vertexList vertex = map_->getVertex(Size(0), jj);
                qglColor(vertex.color1);
                glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                //vertexList normal = map_->getNormals(0, jj);
                //glNormal3d(normal.x1, normal.y1, normal.z1);
                
                vertex = map_->getVertex(Size(0), jj+1);
                qglColor(vertex.color1);
                glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                //normal = map_->getNormals(0, jj+1);
                //glNormal3d(normal.x1, normal.y1, normal.z1);              	        
		              
		            for(Size ii=0; ii<cols-1; ++ii)
		            {
                  vertex = map_->getVertex(ii+1, jj);
                  qglColor(vertex.color1);
                  glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                  //vertexList normal = map_->getNormals(ii+1, jj);
                  //glNormal3d(normal.x1, normal.y1, normal.z1);                
                  
                  vertex = map_->getVertex(ii+1, jj+1);
                  qglColor(vertex.color1);
                  glVertex3d(vertex.x1, vertex.y1, vertex.z1);					
                  //normal = map_->getNormals(ii+1, jj+1);
                  //glNormal3d(normal.x1, normal.y1, normal.z1);                   		        
		            }
		            glEnd();
		          }				        
		        }
		      }
          
          else if(pAspect=="pseudogel")
		      {
  cout << "4" << endl;				    
		        // glLineWidth(layer.param.getValue("dot:line_width"));
		        for(Size jj=0; jj<rows-1; ++jj)
		        {
		          glBegin(GL_LINE_STRIP);
		          
              vertexList vertex = map_->getVertex(Size(0), jj);
              qglColor(vertex.color1);
              glVertex3d(vertex.x1, vertex.y1, vertex.z1);
              
              vertex = map_->getVertex(Size(0), jj+1);
              qglColor(vertex.color1);
              glVertex3d(vertex.x1, vertex.y1, vertex.z1);				        
		            
		          for(Size ii=0; ii<cols-1; ++ii)
		          {
                vertex = map_->getVertex(ii+1, jj);
                qglColor(vertex.color1);
                glVertex3d(vertex.x1, vertex.y1, vertex.z1);
                
                vertex = map_->getVertex(ii+1, jj+1);
                qglColor(vertex.color1);
                glVertex3d(vertex.x1, vertex.y1, vertex.z1);					  		        
		          }
		          glEnd();
		        }
		      }
	      }
	      else
	      {
          // glLineWidth(layer.param.getValue("dot:line_width"));
		      for(
		        Spectrum3DCanvas::ExperimentType::ConstAreaIterator it = layer.peaks.areaBeginConst(
		          canvas_3d_.visible_area_.min_[1],
		          canvas_3d_.visible_area_.max_[1],
		          canvas_3d_.visible_area_.min_[0],
		          canvas_3d_.visible_area_.max_[0]); 
				    it != layer.peaks.areaEndConst(); 
				    ++it)
		      {
			      PeakIndex pi = it.getPeakIndex();
			      if (layer.filters.passes(layer.peaks[pi.spectrum],pi.peak))
			      {
				      glBegin(GL_LINES);
				      switch (canvas_3d_.intensity_mode_)
				      {
					      case SpectrumCanvas::IM_PERCENTAGE:	
						      qglColor( layer.gradient.precalculatedColorAt(0));
						      glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()), 
											       -corner_,
											       -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						      qglColor( layer.gradient.precalculatedColorAt(it->getIntensity() * 100.0 /canvas_3d_.getMaxIntensity(iLayer)) );
						      glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()),
											       -corner_+(GLfloat)scaledIntensity(it->getIntensity(),iLayer),
											       -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						      break;
					
					      case SpectrumCanvas::IM_NONE:
					
						      qglColor( layer.gradient.precalculatedColorAt(canvas_3d_.overall_data_range_.min_[2]));
						      glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()), 
											       -corner_,
											       -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						      qglColor( layer.gradient.precalculatedColorAt(it->getIntensity()));
						      glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()),
											       -corner_+(GLfloat)scaledIntensity(it->getIntensity(),iLayer),
											       -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						      break;
					
					      case SpectrumCanvas::IM_SNAP:
						
						      qglColor(layer.gradient.precalculatedColorAt(int_scale_.min_[0]));
						      glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()), 
											       -corner_,
											       -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						      qglColor(layer.gradient.precalculatedColorAt(it->getIntensity()));
						      glVertex3d(-corner_+(GLfloat)scaledMZ(it->getMZ()),
											       -corner_+(GLfloat)scaledIntensity(it->getIntensity(),iLayer),
											       -near_-2*corner_-(GLfloat)scaledRT(it.getRT()));
						
						      break;
						
				      }
				      glEnd();
			      }
		      }			  
	      }*/
	
	GLuint Spectrum3DOpenGLCanvas::makeGridLines()
	{
		GLuint list = glGenLists(1);
		glNewList(list, GL_COMPILE); 
		glEnable (GL_LINE_STIPPLE);
		glLineStipple (1, 0x0101);	
	 	glBegin(GL_LINES);
		glColor4ub(0, 0, 0, 80);	
		//mz
		if(grid_mz_.size()>=1)
		{
			for (Size i = 0;i<grid_mz_[0].size();i++)
			{
				glVertex3d(-corner_+scaledMZ(grid_mz_[0][i]), -corner_, -near_-2*corner_);
				glVertex3d(-corner_+scaledMZ(grid_mz_[0][i]), -corner_, -far_+2*corner_);
			}
		}
		if(grid_mz_.size()>=2)
		{
			for (Size i = 0;i<grid_mz_[1].size();i++)
			{	
				glVertex3d(-corner_+scaledMZ(grid_mz_[1][i]), -corner_, -near_-2*corner_);
				glVertex3d(-corner_+scaledMZ(grid_mz_[1][i]), -corner_, -far_+2*corner_);
			}
		}
		if(grid_mz_.size()>=3)
		{
			for (Size i = 0;i<grid_mz_[2].size();i++)
			{	
				glVertex3d(-corner_+scaledMZ(grid_mz_[2][i]), -corner_, -near_-2*corner_);
				glVertex3d(-corner_+scaledMZ(grid_mz_[2][i]), -corner_, -far_+2*corner_);
			}
		}
		//rt
		if(grid_rt_.size()>=1)
		{
			for (Size i = 0;i<grid_rt_[0].size();i++)
			{
				glVertex3d(-corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[0][i]));
				glVertex3d( corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[0][i]));
			}
		}
		if(grid_rt_.size()>=2)
		{
			for (Size i = 0;i<grid_rt_[1].size();i++)
			{
				glVertex3d(-corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[1][i]));
				glVertex3d( corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[1][i]));
			}
		}
		if(grid_rt_.size()>=3)
		{
			for (Size i = 0;i<grid_rt_[2].size();i++)
			{
				glVertex3d(-corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[2][i]));
				glVertex3d( corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[2][i]));
			}
		}
		glEnd();	
		glDisable (GL_LINE_STIPPLE);
		glEndList();	
		return list; 
	}
	
	GLuint Spectrum3DOpenGLCanvas::makeAxesTicks()
	{
		GLuint list = glGenLists(1);
		glNewList(list,GL_COMPILE);
		glShadeModel(GL_FLAT);
		glLineWidth(2.0);
		glBegin(GL_LINES);
		qglColor(Qt::black);
	
		//MZ
		if(grid_mz_.size()>=1)
		{
			for (Size i = 0;i<grid_mz_[0].size();i++)
			{
				glVertex3d(-corner_+scaledMZ(grid_mz_[0][i]), -corner_, -near_-2*corner_);
				glVertex3d( -corner_+scaledMZ(grid_mz_[0][i]), -corner_+4.0, -near_-2*corner_);
			}
		}
		if(grid_mz_.size()>=2)
		{
			for (Size i = 0;i<grid_mz_[1].size();i++)
			{
				glVertex3d(-corner_+scaledMZ(grid_mz_[1][i]), -corner_, -near_-2*corner_);
				glVertex3d( -corner_+scaledMZ(grid_mz_[1][i]), -corner_+3.0, -near_-2*corner_);
			}
		}
		if(grid_mz_.size()>=3)
		{
			for (Size i = 0;i<grid_mz_[2].size();i++)
			{
				glVertex3d(-corner_+scaledMZ(grid_mz_[2][i]),  -corner_,  -near_-2*corner_);
				glVertex3d( -corner_+scaledMZ(grid_mz_[2][i]), -corner_+2.0, -near_-2*corner_);
			}
		}
		
		//RT
		if(grid_rt_.size()>=1)
		{
			for (Size i = 0;i<grid_rt_[0].size();i++)
			{
				glVertex3d(-corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[0][i]));
				glVertex3d( -corner_, -corner_+4.0, -near_-2*corner_-scaledRT(grid_rt_[0][i]));
			}
		}
		if(grid_rt_.size()>=2)
		{
			for (Size i = 0;i<grid_rt_[1].size();i++)
			{
				glVertex3d(-corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[1][i]));
				glVertex3d( -corner_, -corner_+3.0, -near_-2*corner_-scaledRT(grid_rt_[1][i]));
			}
		}	
		if(grid_rt_.size()>=3)
		{
			for (Size i = 0;i<grid_rt_[2].size();i++)
			{
				glVertex3d(-corner_, -corner_, -near_-2*corner_-scaledRT(grid_rt_[2][i]));
				glVertex3d( -corner_, -corner_+2.0, -near_-2*corner_-scaledRT(grid_rt_[2][i]));
			}
		}
	
		switch(canvas_3d_.intensity_mode_)
		{
			case SpectrumCanvas::IM_PERCENTAGE:
				if(grid_intensity_.size()>=1)
				{
					for (Size i = 0;i<grid_intensity_[0].size();i++)
					{	
						glVertex3d(-corner_,  -corner_+(2.0 *grid_intensity_[0][i]),  -near_-2*corner_);
						glVertex3d( -corner_+4.0, -corner_+(2.0 * grid_intensity_[0][i]), -near_-2*corner_-4.0);
					}
				}
				break;
				
			case SpectrumCanvas::IM_NONE:
			case SpectrumCanvas::IM_SNAP:
				if(grid_intensity_.size()>=1)
				{
					for (Size i = 0;i<grid_intensity_[0].size();i++)
					{	
						glVertex3d(-corner_,  -corner_+scaledIntensity(grid_intensity_[0][i],canvas_3d_.current_layer_), -near_-2*corner_);
						glVertex3d( -corner_+4.0, -corner_+scaledIntensity(grid_intensity_[0][i],canvas_3d_.current_layer_), -near_-2*corner_-4.0);
					}
				}
				if(grid_intensity_.size()>=2)
				{
					for (Size i = 0;i<grid_intensity_[1].size();i++)
					{
						glVertex3d(-corner_, -corner_+scaledIntensity(grid_intensity_[1][i],canvas_3d_.current_layer_),  -near_-2*corner_);
						glVertex3d( -corner_+3.0, -corner_+scaledIntensity(grid_intensity_[1][i],canvas_3d_.current_layer_), -near_-2*corner_-3.0);
					}
				}
				if(grid_intensity_.size()>=3)
				{
					for (Size i = 0;i<grid_intensity_[2].size();i++)
					{ 
						glVertex3d(-corner_, -corner_+scaledIntensity(grid_intensity_[2][i],canvas_3d_.current_layer_), -near_-2*corner_);
						glVertex3d( -corner_+2.0, -corner_+scaledIntensity(grid_intensity_[2][i],canvas_3d_.current_layer_), -near_-2*corner_-2.0);
					}
				}
				break;
		}
		glEnd();
		glEndList();
		return list; 
	}
	
	double Spectrum3DOpenGLCanvas::scaledRT(double rt)
	{
		double scaledrt = rt - canvas_3d_.visible_area_.min_[1];
		scaledrt = scaledrt * 2.0 * corner_ /(canvas_3d_.visible_area_.max_[1]-canvas_3d_.visible_area_.min_[1]);
		return scaledrt;
	}
	
	double Spectrum3DOpenGLCanvas::scaledInversRT(double rt)
	{
		double i_rt =(rt* canvas_3d_.visible_area_.max_[1] -canvas_3d_.visible_area_.min_[1]*rt);
		i_rt = i_rt/200.0;
		i_rt = i_rt + canvas_3d_.visible_area_.min_[1]; 	
		//	cout<<"rt"<<rt<<"  "<<"scaledinver"<<i_rt<<endl;
		return i_rt;
	 }
	
	double Spectrum3DOpenGLCanvas::scaledMZ(double mz)
	{
		double scaledmz = mz - canvas_3d_.visible_area_.min_[0];
		scaledmz = scaledmz * 2.0 * corner_/(canvas_3d_.visible_area_.max_[0]-canvas_3d_.visible_area_.min_[0]);
		return scaledmz;
	}
	
	double Spectrum3DOpenGLCanvas::scaledInversMZ(double mz)
	{
		double i_mz =(mz *canvas_3d_.visible_area_.max_[0] - mz *canvas_3d_.visible_area_.min_[0]);
		i_mz = i_mz/200;
		i_mz = i_mz + canvas_3d_.visible_area_.min_[0]; 
		return i_mz;
	}
	
	double Spectrum3DOpenGLCanvas::scaledIntensity(Real intensity,Size layer_index)
	{
		double scaledintensity= 0;
		switch(canvas_3d_.intensity_mode_)
		{
			case  SpectrumCanvas::IM_SNAP:
				scaledintensity = intensity -int_scale_.min_[0];
				scaledintensity = ( scaledintensity * 2.0 * corner_)/(int_scale_.max_[0]-int_scale_.min_[0]);
				break;
			case  SpectrumCanvas::IM_NONE:
				scaledintensity = intensity -canvas_3d_.overall_data_range_.min_[2];
				scaledintensity = ( scaledintensity * 2.0 * corner_)/(canvas_3d_.overall_data_range_.max_[2]-canvas_3d_.overall_data_range_.min_[2]);
				break;
			case  SpectrumCanvas::IM_PERCENTAGE: 
		 		scaledintensity =  intensity * 100.0 /canvas_3d_.getMaxIntensity(layer_index);
				scaledintensity = scaledintensity * 2.0 * corner_/100.0;	
				break;
		}
		return scaledintensity;
	}
	
	void Spectrum3DOpenGLCanvas::normalizeAngle(int *angle)
	{
		while (*angle < 0)
		{
			*angle += 360 * 16;
		}
		while (*angle > 360 * 16)
		{
			*angle -= 360 * 16;
		}
	}
	
	///////////////wheel- and MouseEvents//////////////////

	void Spectrum3DOpenGLCanvas::actionModeChange()
	{
		//change from translate to zoom
		if (canvas_3d_.action_mode_== SpectrumCanvas::AM_ZOOM)
		{
			storeRotationAndZoom();
			setAngels(220,220,0);
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		}
		//change from zoom to translate
		else if (canvas_3d_.action_mode_== SpectrumCanvas::AM_TRANSLATE)
		{
			// if still in selection mode, quit selection mode first:
			if(canvas_3d_.rubber_band_.isVisible())
			{
				computeSelection();
			}
			restoreRotationAndZoom();
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		}
	}

  void Spectrum3DOpenGLCanvas::focusOutEvent(QFocusEvent* e)
  {
		canvas_3d_.focusOutEvent(e);  	
  }

	void Spectrum3DOpenGLCanvas::mousePressEvent (QMouseEvent* e)
	{
		mouse_move_begin_ = e->pos();
		mouse_move_end_ = e->pos();	

		if (canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM && e->button()==Qt::LeftButton)
		{
			canvas_3d_.rubber_band_.setGeometry(QRect(e->pos(),QSize()));
			canvas_3d_.rubber_band_.show();
			canvas_3d_.update_buffer_ = true;
			canvas_3d_.update_(__PRETTY_FUNCTION__);
		}
	}

	void Spectrum3DOpenGLCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		if (e->buttons() & Qt::LeftButton)
		{
		  if(canvas_3d_.action_mode_==SpectrumCanvas::AM_ZOOM)
		  {
				canvas_3d_.rubber_band_.setGeometry(QRect(mouse_move_begin_, e->pos()).normalized());
				canvas_3d_.update_(__PRETTY_FUNCTION__);
			}
			else if(canvas_3d_.action_mode_==SpectrumCanvas::AM_TRANSLATE)
		  {
				Int x_angle = xrot_ + 8 * ( e->y() - mouse_move_end_.y() );
				normalizeAngle(&x_angle);
				xrot_ = x_angle;
	
				Int y_angle = yrot_ + 8 * ( e->x() - mouse_move_end_.x() );
				normalizeAngle(&y_angle);
				yrot_ = y_angle;
				
				drawAxesLegend();
				
				mouse_move_end_ = e->pos();
				canvas_3d_.update_(__PRETTY_FUNCTION__);
			}
		}
	}
	
	void Spectrum3DOpenGLCanvas::mouseReleaseEvent (QMouseEvent* e)
	{
		if(canvas_3d_.action_mode_ == SpectrumCanvas::AM_ZOOM && e->button()==Qt::LeftButton)
		{				
			computeSelection();
		}
	}
	
	void Spectrum3DOpenGLCanvas::computeSelection()
	{
		QRect rect = canvas_3d_.rubber_band_.geometry();
		x_1_ = ((rect.topLeft().x()- width_/2) * corner_ *1.25* 2) / width_;
		y_1_ = -300 + (((rect.topLeft().y()-heigth_/2) * corner_*1.25* 2) / heigth_);
		x_2_ = ((rect.bottomRight().x()- width_/2) * corner_ *1.25* 2) / width_;
		y_2_ = -300 + (((rect.bottomRight().y()-heigth_/2) * corner_*1.25* 2) / heigth_);
		dataToZoomArray(x_1_, y_1_, x_2_, y_2_);
		canvas_3d_.rubber_band_.hide();
		canvas_3d_.update_buffer_ = true;
		canvas_3d_.update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DOpenGLCanvas::dataToZoomArray(double x_1, double y_1, double x_2, double y_2)
	{
		double scale_x1 = scaledInversMZ(x_1+100.0);
		double scale_x2 = scaledInversMZ(x_2+100.0);
		double scale_y1 = scaledInversRT(-200-y_1);
		double scale_y2 = scaledInversRT(-200-y_2);
		DRange<2> new_area_;
		if(scale_x1<=scale_x2)
		{
			new_area_.min_[0]= scale_x1;
			new_area_.max_[0]= scale_x2;
		} 
		else
		{
			new_area_.min_[0]= scale_x2;
			new_area_.max_[0]= scale_x1;
		}
		if(scale_y1<=scale_y2)
		{
			new_area_.min_[1]= scale_y1;
			new_area_.max_[1]= scale_y2;
		} 
		else
		{
		 new_area_.min_[1]= scale_y2;
		 new_area_.max_[1]= scale_y1;
		} 
		canvas_3d_.changeVisibleArea_(new_area_, true, true);
	}
	
	void Spectrum3DOpenGLCanvas::updateIntensityScale()
	{
		int_scale_.min_[0]= canvas_3d_.overall_data_range_.max_[2];
		int_scale_.max_[0]= canvas_3d_.overall_data_range_.min_[2];
		
		for (Size i =0;i<canvas_3d_.getLayerCount();i++)
		{
			for (SpectrumCanvas::ExperimentType::ConstIterator spec_it = canvas_3d_.getLayer(i).peaks.RTBegin(canvas_3d_.visible_area_.min_[1]); 
					 spec_it != canvas_3d_.getLayer(i).peaks.RTEnd(canvas_3d_.visible_area_.max_[1]); 
					 ++spec_it)
			{
				for (SpectrumCanvas::ExperimentType::SpectrumType::ConstIterator it = spec_it->MZBegin(canvas_3d_.visible_area_.min_[0]); it!=spec_it->MZEnd(canvas_3d_.visible_area_.max_[0]); ++it)
				{
					if(	int_scale_.min_[0]>= it->getIntensity())
					{
						int_scale_.min_[0]= it->getIntensity(); 
					}
					if(	int_scale_.max_[0]<= it->getIntensity())
					{
						int_scale_.max_[0]= it->getIntensity(); 
					} 
				}
			}
		}
	}
	
	void Spectrum3DOpenGLCanvas::recalculateDotGradient_(Size layer)
	{
		canvas_3d_.getLayer_(layer).gradient.fromString(canvas_3d_.getLayer_(layer).param.getValue("dot:gradient"));
		switch(canvas_3d_.intensity_mode_)
		{
			case SpectrumCanvas::IM_SNAP:
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(int_scale_.min_[0], int_scale_.max_[0], UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				break;
			case SpectrumCanvas::IM_NONE:
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(canvas_3d_.overall_data_range_.min_[2], canvas_3d_.overall_data_range_.max_[2], UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				break;
			case SpectrumCanvas::IM_PERCENTAGE:
				canvas_3d_.getLayer_(layer).gradient.activatePrecalculationMode(0.0, 100.0, UInt(canvas_3d_.param_.getValue("dot:interpolation_steps")));
				break;
		}
	}

}//end of namespace
