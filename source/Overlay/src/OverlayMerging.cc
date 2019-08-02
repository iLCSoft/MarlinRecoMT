#include <MarlinRecoMT/OverlayMerging.h>
#include <MarlinRecoMT/LCIOHelper.h>

// -- marlin headers
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- std headers
#include <string>
#include <map>

// -- lcio headers
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <Exceptions.h>

namespace marlinreco_mt {

  inline long long cellIDToLong( int cellID0, int cellID1 ) {
    return ((long long) cellID0 << 32) | cellID1 ;
  }

  //--------------------------------------------------------------------------

  void OverlayMerging::merge( const EVENT::LCEvent *src, EVENT::LCEvent *dst, const CollectionMap &mergeMap ) {
    for( auto iter : mergeMap ) {
      EVENT::LCCollection *srcCollection = nullptr ;
      EVENT::LCCollection *dstCollection = nullptr ;
      // get the source collection
      try {
        srcCollection = src->getCollection( iter.first ) ;
      } 
      catch ( EVENT::DataNotAvailableException&) {
        continue ;
      }
      // get or create the destination collection
      try {
        dstCollection = dst->getCollection( iter.second ) ;
      } 
      catch ( EVENT::DataNotAvailableException&) {
        streamlog_out( DEBUG ) << "destination collection " << iter.second  << " was created." << std::endl ; 
        dstCollection = new IMPL::LCCollectionVec( srcCollection->getTypeName() ) ;
        // merge collection parameters
        LCIOHelper::mergeLCParameters( srcCollection->getParameters(), dstCollection->parameters() ) ;
        dst->addCollection( dstCollection, iter.second ) ;
      }
      dstCollection->setFlag( srcCollection->getFlag() ) ;
      OverlayMerging::mergeCollections( srcCollection, dstCollection ) ;
    }
  }

  //--------------------------------------------------------------------------

  void OverlayMerging::mergeCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) {
    auto dstType = dst->getTypeName() ;
    // check if collections have the same type
    if ( dstType != src->getTypeName() ) {
      throw EVENT::Exception( "OverlayMerging::mergeCollections: collection types are different" ) ;
    }
    if( dstType == EVENT::LCIO::MCPARTICLE ) {
      OverlayMerging::mergeMCParticleCollections( src, dst ) ;
    }
    else if( dstType == EVENT::LCIO::SIMCALORIMETERHIT ) {
      OverlayMerging::mergeSimCalorimeterHitCollections( src, dst ) ;
    }
    else if( dstType == EVENT::LCIO::CALORIMETERHIT ) {
      OverlayMerging::mergeCalorimeterHitCollections( src, dst ) ;
    }
    else {
      OverlayMerging::mergeAnyCollections( src, dst ) ;
    }
  }

  //--------------------------------------------------------------------------

  void OverlayMerging::mergeMCParticleCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) {
    if( ( src->getTypeName() != EVENT::LCIO::MCPARTICLE ) or ( dst->getTypeName() != EVENT::LCIO::MCPARTICLE ) ) {
      throw EVENT::Exception( "OverlayMerging::mergeMCParticleCollections: not MCParticle collections" ) ;
    }
    int nelts = src->getNumberOfElements();
    for( int i=nelts-1 ; i>=0 ; i-- ) {
      IMPL::MCParticleImpl* p =  dynamic_cast<IMPL::MCParticleImpl*>( src->getElementAt(i) ) ;
      p->setOverlay( true ) ;
      dst->addElement( p ) ;
      src->removeElementAt( i ) ;
    }
  }

  //--------------------------------------------------------------------------

  void OverlayMerging::mergeSimCalorimeterHitCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) {
    int neltsSrc = src->getNumberOfElements();
    int neltsDst = dst->getNumberOfElements();
    // create a map of dest Collection
    std::map<long long, IMPL::SimCalorimeterHitImpl*> dstMap {} ;
    for ( int i=0 ; i<neltsDst ; i++ ) {
      auto dstHit = dynamic_cast<IMPL::SimCalorimeterHitImpl*> ( dst->getElementAt(i) );
      dstMap.insert( std::pair<long long, IMPL::SimCalorimeterHitImpl*>(
        cellIDToLong( dstHit->getCellID0(), dstHit->getCellID1()),
        dstHit )
      ) ;
    }
    // process the src collection and merge with dest
    for ( int i=neltsSrc-1 ; i>=0 ; i-- ) {
      auto srcHit = dynamic_cast<IMPL::SimCalorimeterHitImpl*> ( src->getElementAt(i) ) ;
      auto findIter = dstMap.find( cellIDToLong( srcHit->getCellID0(), srcHit->getCellID1() ) ) ;
      if ( findIter == dstMap.end() ) {
        dst->addElement( srcHit ) ;
      }
      else {
        int numMC = srcHit->getNMCContributions() ;
        for( int j=0 ; j<numMC ; j++ ) {
          findIter->second->addMCParticleContribution( srcHit->getParticleCont(j), srcHit->getEnergyCont(j), srcHit->getTimeCont(j), srcHit->getPDGCont(j));
        }
        delete srcHit;
      }
      src->removeElementAt( i ) ;
    }
  }
  
  //--------------------------------------------------------------------------
  
  void OverlayMerging::mergeCalorimeterHitCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) {
    int neltsSrc = src->getNumberOfElements();
    int neltsDst = dst->getNumberOfElements();
    // create a map of dest Collection
    std::map<long long, IMPL::CalorimeterHitImpl*> dstMap {} ;
    for ( int i=0 ; i<neltsDst ; i++ ) {
      auto dstHit = dynamic_cast<IMPL::CalorimeterHitImpl*> ( dst->getElementAt(i) );
      dstMap.insert( std::pair<long long, IMPL::CalorimeterHitImpl*>(
        cellIDToLong( dstHit->getCellID0(), dstHit->getCellID1()),
        dstHit )
      ) ;
    }
    // process the src collection and merge with dest
    for ( int i=neltsSrc-1 ; i>=0 ; i-- ) {
      auto srcHit = dynamic_cast<IMPL::CalorimeterHitImpl*> ( src->getElementAt(i) ) ;
      auto findIter = dstMap.find( cellIDToLong( srcHit->getCellID0(), srcHit->getCellID1() ) ) ;
      if ( findIter == dstMap.end() ) {
        dst->addElement( srcHit ) ;
      }
      else {
        findIter->second->setEnergy( findIter->second->getEnergy() + srcHit->getEnergy() );
      }
      src->removeElementAt( i ) ;
    }
  }

  //--------------------------------------------------------------------------

  void OverlayMerging::mergeAnyCollections( EVENT::LCCollection* src, EVENT::LCCollection* dst ) {
    int nelts = src->getNumberOfElements() ;
    for( int i=nelts-1 ; i>=0 ; i-- ) {
      auto elt = src->getElementAt(i) ;
      dst->addElement( elt ) ;
      src->removeElementAt( i ) ;
    }
  }

}
