#include <MarlinRecoMT/ErrorOfSigma.h>

// -- std headers
#include <iostream>
#include <cmath>

namespace marlinreco_mt {

  ErrorOfSigma::ErrorOfSigma( unsigned n ) :
    _n( n ) {
    if( n < 30 ) { // FIXME - implement proper errors for small n
      std::cout << "ErrorOfSigma::ErrorOfSigma: errors will be inaccurate for n = " << n
		<< " - should be > 30 "
		<< std::endl ;
    }
  }

  double ErrorOfSigma::lowerError( double sigma ) {
    return ( 1. - std::sqrt( _n - 1.0 ) / std::sqrt( getChiSquaredPlus()  ) ) * sigma ;
  }

  double ErrorOfSigma::upperError( double sigma ) {
    return ( std::sqrt( _n - 1.0 ) / std::sqrt( getChiSquaredMinus() ) - 1.  ) * sigma ;
  }


  double ErrorOfSigma::getChiSquaredPlus() {
    double t = 2./ ( 9. * _n )  ;
    double chiPlus = _n  * std::pow( ( 1. - t + std::sqrt( t) ) , 3. ) ;
    return chiPlus ;
  }

  double ErrorOfSigma::getChiSquaredMinus() {
    double t = 2./ ( 9. * _n )  ;
    double chiMinus = _n  * std::pow( ( 1. - t -  std::sqrt( t) ) , 3. ) ;
    return chiMinus ;
  }

}
