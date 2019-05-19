#include <MarlinRecoMT/SimpleClusterSmearer.h>

// -- std headers
#include <cmath>

namespace marlinreco_mt {

  SimpleClusterSmearer::SimpleClusterSmearer( const std::vector<float>& resVec ) :
    _resVec(0) {
    // copy the resolution vector parameters into a more structured vector
    const unsigned int size = resVec.size() / ( sizeof( ClusterResolution)  / sizeof(float) );  // ==3
    _resVec.reserve(size);
    int index = 0 ;
    for( unsigned int i=0 ; i <  size ; i++ ){
      float A     =  resVec[ index++ ] ;
      float B     =  resVec[ index++ ] ;
      float thMin =  resVec[ index++ ] ;
      float thMax =  resVec[ index++ ] ;
      _resVec.push_back( ClusterResolution( A, B , thMin, thMax ) );
    }
  }


  TLorentzVector SimpleClusterSmearer::smearedFourVector( const TLorentzVector& v, int /*pdgCode*/ ) {
    // find resolution for polar angle
    double theta = v.Theta() ;
    if( theta > M_PI_2 ) {
      theta = M_PI - theta ; // need to transform to [0,pi/2]
    }
    std::pair<double,double> resolution = std::make_pair( -1., -1. ) ;
    for( unsigned int i=0 ; i <  _resVec.size() ; i++ ) {
      if( theta <= _resVec[i].ThMax  &&  theta > _resVec[i].ThMin ) {
      	resolution.first =  _resVec[i].A ;
      	resolution.second =  _resVec[i].B ;
      	break ;
      }
    }
    TLorentzVector sv( 0., 0., 0., 0. ) ;

    if( resolution.first > - 1e-10  ) {
      // do the smearing ....
      double E = v.E() ;
      double Eres = std::sqrt( resolution.first * resolution.first +
				resolution.second * resolution.second / E  )  ;
      std::normal_distribution<float> gaus( 0, E*Eres ) ;
      double deltaE = gaus( _ranEngine ) ;
      // assume massless clusters ...
      auto n3v = v.Vect() ;
      n3v.SetMag( E + deltaE ) ;
      sv.SetVectM( n3v, 0. ) ;
    }
    return sv ;
  }

}
