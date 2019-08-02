#ifndef MARLINRECOMT_OVERLAYMERGING_H
#define MARLINRECOMT_OVERLAYMERGING_H 1

// -- std headers
#include <string>
#include <map>

// -- lcio headers
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>

namespace marlinreco_mt {

  class OverlayMerging {
    // static API only
    OverlayMerging() = delete ;
  public:
    using CollectionMap = std::map<std::string, std::string> ;

  public:

    static void merge( const EVENT::LCEvent *src, EVENT::LCEvent *dst, const CollectionMap &mergeMap ) ;

    static void mergeCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;

    static void mergeMCParticleCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;

    static void mergeSimCalorimeterHitCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;
    
    static void mergeCalorimeterHitCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;

    static void mergeAnyCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;

  };

}

#endif
