#include <MarlinRecoMT/LCLine3D.h>
#include <MarlinRecoMT/LCPlane3D.h>

#include <iostream>
#include <cmath>
#include <float.h>
#include <exception>

namespace marlinreco_mt {

  LCLine3D::LCLine3D()
  {
    _reference.SetXYZ(0.,0.,0.);
    _point.SetXYZ(0.,0.,0.);
    _direction.SetXYZ(1.,0.,0.);
  }

  LCLine3D::LCLine3D(const LCVector3D & point, const LCVector3D & lineDirection)
  {
    set( point, lineDirection, LCVector3D(0.,0.,0.) );
  }

  LCLine3D::LCLine3D(const LCVector3D & point,
  		   const LCVector3D & lineDirection,
  		   const LCVector3D & reference) 
  {
    set( point, lineDirection, reference );
  }

  LCLine3D::LCLine3D(double d0, double phi0, double z0, double tanLambda) 
  {
    set( d0, phi0, z0, tanLambda, LCVector3D(0.,0.,0.) );
  }

  LCLine3D::LCLine3D(double d0, double phi0, double z0, double tanLambda,
  	 const LCVector3D & reference) 
  {
    set( d0, phi0, z0, tanLambda, reference );
  }

  LCLine3D::LCLine3D(const LCLine3D & line) 
  {
    _point     = line._point;
    _direction = line._direction;
    _reference = line._reference;
  }

  bool LCLine3D::set(const LCVector3D & point,
  		   const LCVector3D & lineDirection,
  		   const LCVector3D & reference) 
  {
    //  std::cout << "LCLine  " << _reference << " " << point << " " << direction << std::endl;

    _reference = reference;
    _direction = lineDirection.unit();
    if (_direction.mag2() == 0)
      {
        return false;
      }

  // calculate _point to be the PCA to the reference point according to the
  // definition given in LC-LC-DET-2006-004:
  // the x,y compnents have to beh teh  PCA, the z compnente is calculated 
  // after that. 

    LCVector3D p = point, d = _direction;
    p.SetZ(0.);
    d.SetZ(0.); 
    auto mag = std::sqrt(d.mag2());

    if (mag !=  0.)
      {
        double sFaktor = 1./mag;
        d = d.unit();
        double s = ( - p.Dot(d) ) / d.mag2() ;
        // x,y componentes:
        //      _point = p + s * d ;
        // z component:
        _point = ( (point + _direction*s*sFaktor) );
      }
    else
      {
        _point = point;
        _point.SetZ(0.);
      }
    //  std::cout << "LCLine: " << _reference << " " << _point << " " << _direction << std::endl;
    return true;
  }

  bool LCLine3D::set(double d0, double phi0, double z0, double tanLambda,
  	 const LCVector3D & reference) 
  {
    _reference = reference;
    _direction.SetXYZ( cos(phi0), sin(phi0), tanLambda );
    _direction = _direction.unit();
    if (d0 == 0.)
      {
        _point.SetXYZ(0.,0.,z0);
      }
    else
      {
        _point.SetXYZ( ( d0*sin(phi0) ), ( d0*cos(phi0) ), z0 );
      }
    return true;
  }

  LCLine3D & LCLine3D::operator=(const LCLine3D & rhs) 
  {
    _point     = rhs._point;
    _direction = rhs._direction;
    _reference = rhs._reference;

    return *this;
  }

  LCVector3D LCLine3D::position(const double s) const 
  {
    return (_reference+_point + s*_direction) ;
  }

  LCVector3D LCLine3D::direction() const 
  {
    return _direction;
  }

  double LCLine3D::distance(const LCVector3D & point) const 
  {
    return std::sqrt(( point - position( projectPoint( point ) ) ).mag2()) ;
  }

  double LCLine3D::projectPoint(const LCVector3D & point) const 
  {
    // the last therm : (...) / _direction.mag2() is not there becaus 
    // the _direction vector is normalised.
    //  return ( 2*point*_direction - _point*_direction ) / _direction.mag2() ;
    //  return ( point*_direction - (_reference+_point)*_direction ) / _direction.mag2() ;
    double x = ( point.Dot(_direction) - (_reference+_point).Dot(_direction) ) / _direction.mag2() ;
    //  std::cout << "x: " << x << std::endl;
    //  std::cout << "point: " << point 
    //	    << " direction: " << _direction << " _point: " << _point 
    //	    << " d.mag: " << _direction.mag2() << std::endl;


    return x;
  }

  bool LCLine3D::operator==(const LCLine3D & rhs) const 
  {
    return (_point == rhs._point &&
  	  _direction == rhs._direction && 
  	  _reference == rhs._reference);
  }

  bool LCLine3D::operator!=(const LCLine3D & rhs) const 
  {
    return (_point != rhs._point ||
  	  _direction != rhs._direction ||
  	  _reference != rhs._reference) ;
  }

  double LCLine3D::intersectionWithPlane(const LCPlane3D plane, bool& pointExists) const 
  {
    double c = direction().Dot(plane.normal()) ;

    if (c == 0)
      { // no interaction 
        pointExists = false;
        return DBL_MAX;
      }

    pointExists = true;
    return - ( position().Dot(plane.normal()) + plane.d() ) / c ;
  }

  std::ostream & operator << (std::ostream &os, const LCLine3D &l)
  {
    return os << l.position() << "+s*" << l.direction() ;
  }
  
}
