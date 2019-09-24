#include <MarlinRecoMT/Circle.h>

// -- marlin headers
#include <marlin/Exceptions.h>

// -- std headers
#include <cmath>

namespace marlinreco_mt {

  Circle::Circle( const dd4hep::rec::Vector2D &pt1, const dd4hep::rec::Vector2D &pt2, const dd4hep::rec::Vector2D &pt3 ) {  	
  	if (!this->isPerpendicular(pt1, pt2, pt3) ) {
      this->calculateCircleProperties(pt1, pt2, pt3) ;
  	}
    else if (!this->isPerpendicular(pt1, pt3, pt2) ) {
      this->calculateCircleProperties(pt1, pt3, pt2) ;
  	}
    else if (!this->isPerpendicular(pt2, pt1, pt3) ) {
      this->calculateCircleProperties(pt2, pt1, pt3) ;	
  	}
  	else if (!this->isPerpendicular(pt3, pt2, pt1) ) {
      this->calculateCircleProperties(pt3, pt2, pt1) ;	
    }
    else if (!this->isPerpendicular(pt2, pt3, pt1) ) {
      this->calculateCircleProperties(pt2, pt3, pt1) ;	
  	}
    else if (!this->isPerpendicular(pt3, pt1, pt2) ) {
      this->calculateCircleProperties(pt3, pt1, pt2) ;	
  	}
    else { 
      throw marlin::Exception( "Circle::Circle: Couldn't construct circle from input 2D vectors" ) ;
  	}
  }
  
  //--------------------------------------------------------------------------
  
  double Circle::radius() const {
    return _radius ;
  }
  
  //--------------------------------------------------------------------------
  
  const dd4hep::rec::Vector2D &Circle::center() const {
    return _center ;
  }
  
  //--------------------------------------------------------------------------

  bool Circle::isPerpendicular( const dd4hep::rec::Vector2D &pt1, const dd4hep::rec::Vector2D &pt2, const dd4hep::rec::Vector2D &pt3 ) const {
  	double yDelta_a = pt2.v() - pt1.v();
  	double xDelta_a = pt2.u() - pt1.u();
  	double yDelta_b = pt3.v() - pt2.v();
  	double xDelta_b = pt3.u() - pt2.u();
  	if ( fabs(xDelta_a) <= TOLERANCE && fabs(yDelta_b) <= TOLERANCE ) {
  	  return false;
  	}
  	if ( fabs(yDelta_a) <= TOLERANCE ) {
  		return true;
  	}
  	else if ( fabs(yDelta_b) <= TOLERANCE ) {
  		return true;
  	}
  	else if ( fabs(xDelta_a) <= TOLERANCE ) {
  		return true;
  	}
  	else if ( fabs(xDelta_b) <= TOLERANCE ) {
  		return true;
  	}
  	return false ;
  }
  
  //--------------------------------------------------------------------------

  void Circle::calculateCircleProperties( const dd4hep::rec::Vector2D &pt1, const dd4hep::rec::Vector2D &pt2, const dd4hep::rec::Vector2D &pt3 ) {
    double yDelta_a = pt2.v() - pt1.v() ;
    double xDelta_a = pt2.u() - pt1.u() ;
    double yDelta_b = pt3.v() - pt2.v() ;
    double xDelta_b = pt3.u() - pt2.u() ;
    if ( fabs(xDelta_a) <= TOLERANCE && fabs(yDelta_b) <= TOLERANCE ) {
      _center = dd4hep::rec::Vector2D( 0.5*(pt2.u() + pt3.u()), 0.5*(pt1.v() + pt2.v()) ) ;
      _radius = std::sqrt( (pt1.u()-_center.u())*(pt1.u()-_center.u()) + (pt1.v()-_center.v())*(pt1.v()-_center.v()));
      return ;
    }
    // isPerpendicular() assure that xDelta(s) are not zero
    double aSlope = yDelta_a / xDelta_a ; 
    double bSlope = yDelta_b / xDelta_b ;
    // checking whether the given points are colinear. 	
    if ( fabs(aSlope-bSlope) <= TOLERANCE ) {
      throw marlin::Exception( "Circle::calculateCircleProperties: Couldn't calculate properties, the three input 2D points are colinear" ) ;
  	}
    const double x = (aSlope*bSlope*(pt1.v() - pt3.v()) + bSlope*(pt1.u() + pt2.u()) - aSlope*(pt2.u()+pt3.u()) )/(2* (bSlope-aSlope) ) ;
    const double y = -1*(x - (pt1.u()+pt2.u())/2)/aSlope +  (pt1.v()+pt2.v())/2  ;
    _center =  dd4hep::rec::Vector2D( x, y ) ;
    _radius = std::sqrt((pt1.u()-_center.u())*(pt1.u()-_center.u()) + (pt1.v()-_center.v())*(pt1.v()-_center.v()));
  }
    
}
