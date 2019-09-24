// Circle.h: interface for the Circle class.
// Circle class.
// Purpose : Represent the circle object
// Input : 3 different points
// Process : Calcuate the radius and center
// Output : Circle
//           
// This class originally designed for representation of discretized curvature information 
// of sequential pointlist  
// KJIST CAD/CAM     Ryu, Jae Hun ( ryu@geguri.kjist.ac.kr)
// Last update : 2019.09, R.Ete, DESY

#include <DDRec/Vector2D.h>

namespace marlinreco_mt {
	
	class Circle {
	public:
		/// The tolerance used in the circle computations
		static constexpr double TOLERANCE = 0.000000001 ;
		
	public:
		/// Default constructor
		Circle() = default ;
		
		/// Default destructor
		~Circle() = default ;
		
		/// Constructor with co-planar vectors
		Circle( const dd4hep::rec::Vector2D &p1, const dd4hep::rec::Vector2D &p2, const dd4hep::rec::Vector2D &p3 ) ;
		
		/// Get the circle radius 
		double radius() const ;
		
		/// Get the circle center vector
		const dd4hep::rec::Vector2D &center() const ;

	private:
		void calculateCircleProperties( const dd4hep::rec::Vector2D &p1, const dd4hep::rec::Vector2D &p2, const dd4hep::rec::Vector2D &p3 ) ;
		bool isPerpendicular( const dd4hep::rec::Vector2D &p1, const dd4hep::rec::Vector2D &p2, const dd4hep::rec::Vector2D &p3 ) const ;
		
	private:
		/// The circle radius
		double                   _radius {0.f} ;
		/// The circle center
		dd4hep::rec::Vector2D    _center {} ;
	};
		
}
