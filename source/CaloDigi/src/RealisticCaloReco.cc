// Calorimeter digitiser
#include <MarlinRecoMT/RealisticCaloReco.h>

// -- marlin headers
#include <marlin/ProcessorApi.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- lcio headers
#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCRelationNavigator.h>

// -- marlinrecomt headers
#include <MarlinRecoMT/CalorimeterHitType.h>

// -- std headers
#include <iostream>
#include <string>
#include <assert.h>
#include <cmath>

namespace marlinreco_mt {

  RealisticCaloReco::RealisticCaloReco( const std::string &pname ) : 
    marlin::Processor( pname ) {
    _description = "Performs simple reconstruction of calo hits..." ;
  }
  
  //--------------------------------------------------------------------------

  void RealisticCaloReco::init() {
    printParameters() ;
    if( _inputRelationCollections.get().size() != _inputCollections.get().size()
     || _outputCollections.get().size() != _inputCollections.get().size() 
     || _outputRelationCollections.get().size() != _inputCollections.get().size() 
     || _calibrationCoefficients.get().empty() 
     || _calibrationCoefficients.get().size() != _calibrationLayers.get().size() ) {
      marlin::ProcessorApi::abort( this, "Invalid parameters from steering file. Please check your inputs!" ) ;
    }
  }
  
  //--------------------------------------------------------------------------

  void RealisticCaloReco::processEvent( EVENT::LCEvent *evt ) {
    // common collection flags for all output collections
    IMPL::LCFlagImpl collectionFlag {} ;
    collectionFlag.setBit( EVENT::LCIO::CHBIT_LONG);
    collectionFlag.setBit( EVENT::LCIO::RCHBIT_TIME ) ; // store timing on output hits.
    IMPL::LCFlagImpl relationFlag {} ;
    relationFlag.setBit( EVENT::LCIO::LCREL_WEIGHTED ) ; // for the hit relations
    // * Reading Collections of digitised calorimeter Hits *
    for (unsigned int i(0); i < _inputCollections.get().size(); ++i) {
      std::string colName =  _inputCollections.get()[i] ;
      std::string relName =  _inputRelationCollections.get()[i] ;
      log<DEBUG>() << "looking for hit, relation collection: " << colName << " " << relName << std::endl ;
      try {
        auto collection = evt->getCollection( colName ) ;
        auto relationCollection = evt->getCollection( relName ) ;
        auto cellIDString = collection->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ) ;
        UTIL::LCRelationNavigator navigator( relationCollection ) ;
        UTIL::CellIDDecoder<EVENT::CalorimeterHit> cellIDDecoder ( collection ) ;
        // create new collection
        auto outputCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::CALORIMETERHIT ) ;
        outputCollection->setFlag( collectionFlag.getFlag() ) ;
        // relation from reconstructed to sim hits
        auto relationOutputCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::LCRELATION ) ;
        relationOutputCollection->setFlag( relationFlag.getFlag() ) ;
        relationOutputCollection->parameters().setValue( RELATIONFROMTYPESTR , EVENT::LCIO::CALORIMETERHIT ) ;
        relationOutputCollection->parameters().setValue( RELATIONTOTYPESTR   , EVENT::LCIO::SIMCALORIMETERHIT ) ;

        int numElements = collection->getNumberOfElements();
        log<DEBUG>() << colName << " number of elements = " << numElements << std::endl ;

        for ( int j=0 ; j<numElements ; ++j ) {
          auto hit = static_cast<EVENT::CalorimeterHit*>( collection->getElementAt( j ) ) ;
          // make new hit
        	auto newhit = new IMPL::CalorimeterHitImpl(); 
        	newhit->setCellID0( hit->getCellID0() ) ;
        	newhit->setCellID1( hit->getCellID1() ) ;
        	newhit->setEnergy( this->reconstructEnergy( cellIDDecoder, hit ) ) ; // overloaded method, technology dependent
        	newhit->setRawHit( hit->getRawHit() ) ;
        	newhit->setTime( hit->getTime() ) ;
        	newhit->setPosition( hit->getPosition() ) ;
        	newhit->setType( hit->getType() ) ;
        	outputCollection->addElement( newhit ) ;
        	// get the simcalohit corresponding to this digitised hit
          auto relatedObjects = navigator.getRelatedToObjects( hit ) ;
        	if ( not relatedObjects.empty() ) {
        	  auto simhit = static_cast<EVENT::SimCalorimeterHit*>( relatedObjects[0] ); // assume the first one (should be only one)
        	  // make a relation, add to collection - keep relations from reco to sim hits
        	  relationOutputCollection->addElement( new IMPL::LCRelationImpl( newhit , simhit , 1.0 ) ) ;
        	} 
          else {
        	  log<WARNING>() << "could not find relation to sim calo hit!" << std::endl ;
        	}
        }
        // add collection to event
        outputCollection->parameters().setValue( EVENT::LCIO::CellIDEncoding, cellIDString ) ;
        evt->addCollection( outputCollection.release(),         _outputCollections.get()[i] ) ;
        evt->addCollection( relationOutputCollection.release(), _outputRelationCollections.get()[i] ) ;
      }
      catch( EVENT::DataNotAvailableException &e ) {
        log<DEBUG>() << "could not find input ECAL collection " << colName << std::endl ;
      }
    }
  }
  
  //--------------------------------------------------------------------------

  float RealisticCaloReco::getLayerCalib( int ilayer ) const {
    float calibrationCoefficient = 0 ;
    // retrieve calibration constants
    // Fixed the following logic (DJeans, June 2016)
    int min(0),max(0) ;
    for (unsigned int k(0); k < _calibrationLayers.get().size(); ++k ) {
      if ( k > 0 ) {
        min += _calibrationLayers.get()[k-1] ;
      }
      max += _calibrationLayers.get()[k] ;
      if (ilayer >= min && ilayer < max) {
        calibrationCoefficient = _calibrationCoefficients.get()[k] ;
        break ;
      }
    }
    return calibrationCoefficient ;
  }
  
}

