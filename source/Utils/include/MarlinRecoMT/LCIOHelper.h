#ifndef MARLINRECOMT_LCIOHELPER_h
#define MARLINRECOMT_LCIOHELPER_h 1

#include <EVENT/LCParameters.h>

namespace marlinreco_mt {

  class LCIOHelper {
    // static API only
    LCIOHelper() = delete ;

  public:
    static void mergeLCParameters( const EVENT::LCParameters &src, EVENT::LCParameters &dst ) ;

  };

}

#endif
