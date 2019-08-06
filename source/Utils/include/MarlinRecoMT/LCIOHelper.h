#ifndef MARLINRECOMT_LCIOHELPER_h
#define MARLINRECOMT_LCIOHELPER_h 1

// -- lcio headers
#include <EVENT/LCParameters.h>
#include <EVENT/LCObject.h>

// -- std headers
#include <utility>
#include <stdexcept>

namespace marlinreco_mt {

  class LCIOHelper {
    // static API only
    LCIOHelper() = delete ;

  public:
    static void mergeLCParameters( const EVENT::LCParameters &src, EVENT::LCParameters &dst ) ;
    
    /// Combine cellID0 and cellID1 into a single long variable
    template <typename T>
    static long long cellIDToLong( const EVENT::LCObject *obj ) ;

    /// Combine cellID0 and cellID1 into a single long variable
    static long long cellIDToLong( int cellID0, int cellID1 ) ;
      
    /// Split cellID0 and cellID1 from a single long variable
    static std::pair<int, int> longToCellID( long long l ) ;
  };
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  template <typename T>
  inline long long LCIOHelper::cellIDToLong( const EVENT::LCObject *obj ) {
    auto objCast = dynamic_cast<const T*>( obj ) ;
    if( nullptr == objCast ) {
      throw std::runtime_error("LCIOHelper::cellIDToLong<T>: invalid object cast") ;
    }
    return LCIOHelper::cellIDToLong( objCast->getCellID0(), objCast->getCellID1() ) ;
  }

}

#endif
