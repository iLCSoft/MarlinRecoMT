#ifndef MARLINRECOMT_OVERLAYMERGING_H
#define MARLINRECOMT_OVERLAYMERGING_H 1

// -- std headers
#include <string>
#include <map>

// -- lcio headers
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>

namespace marlinreco_mt {

  /**
   *  @brief  OverlayMerging class
   *          Helper class to merge collections for Overlay processors
   */
  class OverlayMerging {
    // static API only
    OverlayMerging() = delete ;
  public:
    using CollectionMap = std::map<std::string, std::string> ;

  public:
    /**
     *  @brief  Merge two events. Only collections appearing in both events are merged
     *  
     *  @param  src the source event
     *  @param  dst the destination event
     */
    static void mergeEvents( const EVENT::LCEvent *src, EVENT::LCEvent *dst ) ;

    /**
     *  @brief  Merge two events. Only the collections appearing in the provided map
     *          are merged, if of course they appear at least in the source event.
     *          If the collection doesn't exists in the destination event, a new collection 
     *          is added 
     *  
     *  @param  src the source event
     *  @param  dst the destination event
     *  @param  mergeMap the map of collection to merge
     */
    static void mergeEvents( const EVENT::LCEvent *src, EVENT::LCEvent *dst, const CollectionMap &mergeMap ) ;

    /**
     *  @brief  Merge two collections. The merging strategy differs depending 
     *          on the collection type:
     *          - CalorimeterHit: hit energies are added
     *          - SimCalorimeterHit: MCContributions are added
     *          - MCParticle: mc particles are added to the destination colelction and flagged as 'overlay' 
     *          For any other collection types, the element are moved from one 
     *          collection to the other.
     * 
     *  @param  src the source collection
     *  @param  dst the destination collection
     */
    static void mergeCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;

  private:
    // Collection merging for different collection types
    static void mergeMCParticleCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;
    static void mergeSimCalorimeterHitCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;  
    static void mergeCalorimeterHitCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;
    static void mergeAnyCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) ;
  };

}

#endif
