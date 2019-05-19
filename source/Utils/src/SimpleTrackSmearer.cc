#include <MarlinRecoMT/SimpleTrackSmearer.h>

#include <cmath>
#include <cstdlib>

namespace marlinreco_mt {

  SimpleTrackSmearer::SimpleTrackSmearer( const std::vector<float>& resVec ) :
    _resVec(0) {
    const unsigned int size = resVec.size() / ( sizeof(TrackResolution)  / sizeof(float) );  // ==3
    _resVec.reserve(size);
    // copy the resolution vector parameters into a more structured vector
    int index = 0 ;
    for( unsigned int i=0 ; i < size ; i++ ){
      float dPP   =  resVec[ index++ ] ;
      float thMin =  resVec[ index++ ] ;
      float thMax =  resVec[ index++ ] ;
      _resVec.push_back( TrackResolution( dPP, thMin, thMax ) );
    }
  }


  TLorentzVector SimpleTrackSmearer::smearedFourVector( const TLorentzVector& v, int pdgCode ){
    // find resolution for polar angle
    double theta = v.Theta() ;
    if( theta > M_PI_2 ) {
      theta = M_PI - theta ; // need to transform to [0,pi/2]
    }
    double resolution = -1. ;
    for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){
      if( theta <= _resVec[i].ThMax  &&  theta > _resVec[i].ThMin ) {
      	resolution =  _resVec[i].DPP ;
      	break ;
      }
    }
    TLorentzVector sv( 0., 0. , 0., 0. ) ;
    if( resolution > - 1e-10  ) {
      // do the smearing ....
      double P = v.Vect().Mag() ;
      std::normal_distribution<float> gaus ( 0. , P*P*resolution ) ;
      double deltaP = gaus( _ranEngine ) ;
      auto n3v = v.Vect() ;
      n3v.SetMag( P + deltaP ) ;
      // assume perfect electron and muon ID and
      // assign pion mass to everything else
      double mass = PION_MASS ;
      if( std::abs( pdgCode ) == 12 ) { // electron
	       mass = ELECTRON_MASS ;
      }
      else if( std::abs( pdgCode ) == 13 ) { // muon
	       mass = MUON_MASS ;
      }
      sv.SetVectM(  n3v  , mass  ) ;
    }
    return sv ;
  }

}
