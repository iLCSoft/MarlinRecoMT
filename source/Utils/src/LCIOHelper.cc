#include <MarlinRecoMT/LCIOHelper.h>

// -- lcio headers
#include <LCIOTypes.h>

namespace marlinreco_mt {

  void LCIOHelper::mergeLCParameters( const EVENT::LCParameters &src, EVENT::LCParameters &dst ) {
    // merge string parameters
  	EVENT::StringVec strKeys ;
  	src.getStringKeys( strKeys ) ;
  	for( auto key : strKeys ) {
  	  EVENT::StringVec vals ;
  	  src.getStringVals( key , vals ) ;
  	  dst.setValues( key , vals ) ;
  	}
    // merge int parameters
  	EVENT::StringVec intKeys ;
  	src.getIntKeys( intKeys ) ;
  	for( auto key : intKeys ) {
  	  EVENT::IntVec vals ;
  	  src.getIntVals( key , vals ) ;
  	  dst.setValues( key , vals ) ;
  	}
    // merge float parameters
  	EVENT::StringVec floatKeys ;
  	src.getFloatKeys( floatKeys ) ;
  	for( auto key : floatKeys ){
  	  EVENT::FloatVec vals ;
  	  src.getFloatVals( key , vals ) ;
  	  dst.setValues( key , vals ) ;
  	}
  }

}
