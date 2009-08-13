// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2009
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
// $Maintainer: $
// $Authors: $
// -------------------------------------------------------------------------- 

#include <OpenMS/VISUAL/Arrow3d.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
  Arrow3d::Arrow3d()
  {	
    hat = gluNewQuadric();
	  gluQuadricDrawStyle(hat, GLU_FILL);
	  gluQuadricNormals(hat, GLU_SMOOTH);
	  gluQuadricOrientation(hat, GLU_OUTSIDE);
	  
	  disk = gluNewQuadric();
	  gluQuadricDrawStyle(disk, GLU_FILL);
	  gluQuadricNormals(disk, GLU_SMOOTH);
	  gluQuadricOrientation(disk, GLU_OUTSIDE);
	  
	  base = gluNewQuadric();
	  gluQuadricDrawStyle(base, GLU_FILL);
	  gluQuadricNormals(base, GLU_SMOOTH);
	  gluQuadricOrientation(base, GLU_OUTSIDE);
	  
	  bottom = gluNewQuadric();
	  gluQuadricDrawStyle(bottom, GLU_FILL);
	  gluQuadricNormals(bottom, GLU_SMOOTH);
	  gluQuadricOrientation(bottom, GLU_OUTSIDE);
  }
  
  Arrow3d::~Arrow3d()
  {
	  gluDeleteQuadric(hat);
	  gluDeleteQuadric(disk);
	  gluDeleteQuadric(base);
	  gluDeleteQuadric(bottom);
  }
  
  void Arrow3d::draw(const Struct3d& pos, const Struct3d& dir, const double sizeFactor) const
  {
	  glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
      
    color_.glColor();

	  double coneRadius = 0.06 * sizeFactor;
	  double bodyRadius = 0.02 * sizeFactor;
	  double coneHeight = 0.4 * sizeFactor;
	  double bodyHeight = 0.6 * sizeFactor;

	  Struct3d axis;
	  double phi = rotate_(axis, dir);
	
	  pos.glTranslate();
    axis.glRotate(phi);
	  gluDisk(disk, 0.0, bodyRadius, 5, 1);
	  gluCylinder(base, bodyRadius, bodyRadius, bodyHeight, 5, 1);
	  glTranslatef(0, 0, bodyHeight);
	  gluDisk(disk, bodyRadius, coneRadius, 5, 1);	  	  	  
	  gluCylinder(hat, coneRadius, 0.0, coneHeight, 5, 1);

    glPopMatrix();
  }

  void Arrow3d::setColor(ColorRGBA color)
  {
    color_ = color;
  }

  double Arrow3d::rotate_(Struct3d& axis, const Struct3d& vector) const
  {
	  axis = normalCrossProduct(Struct3d(0.0, 0.0, 1.0), vector);
	  double cosphi = dotProduct(Struct3d(0.0, 0.0, 1.0), vector.normalize());
	
	  return 180 * acos(cosphi) / Constants::PI;
  }
}//Namespace

