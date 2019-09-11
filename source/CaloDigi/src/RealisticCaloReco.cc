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

  RealisticCaloReco::RealisticCaloReco( const std::string &pname ) : marlin::Processor( pname ) {

    _description = "Performs simple reconstruction of calo hits..." ;

    registerInputCollections( EVENT::LCIO::CALORIMETERHIT,
  			    "inputHitCollections",
  			    "input hit collection names",
  			    _inputCollections,
  			    _inputCollections ) ;

    registerInputCollections( EVENT::LCIO::LCRELATION,
  			    "inputRelationCollections",
  			    "input relation collection names (digi<->sim), one per inputHitCollection",
  			    _inputRelationCollections,
  			    _inputRelationCollections ) ;

    // output collection names
    registerProcessorParameter( "outputHitCollections",
  			      "output hit collection names",
  			      _outputCollections,
  			      _outputCollections ) ;

    registerProcessorParameter( "outputRelationCollections",
  			      "output hit collection names",
  			      _outputRelationCollections,
  			      _outputRelationCollections ) ;

    registerProcessorParameter("calibration_layergroups" ,
                               "grouping of calo layers" ,
                               _calibrationLayers,
                               _calibrationLayers ) ;

    registerProcessorParameter("calibration_factorsMipGev" ,
                               "Calibration coefficients (MIP->shower GeV) of layers groups" ,
                               _calibrationCoefficients,
                               _calibrationCoefficients ) ;

    registerProcessorParameter("CellIDLayerString" ,
                               "name of the part of the cellID that holds the layer" , 
                               _cellIDLayerString , 
                               std::string("K-1") ) ;
  }
  
  //--------------------------------------------------------------------------

  void RealisticCaloReco::init() {
    printParameters() ;

    // if no output collection names specified, set some default based on the input collection names
    if ( _outputCollections.empty() ) {
      for ( size_t i=0 ; i<_inputCollections.size() ; i++ ) {
        _outputCollections.push_back( _inputCollections[i] + "Reco" ) ;
      }
    }
    if ( _outputRelationCollections.size()==0 ) {
      for (size_t i=0; i<_inputCollections.size(); i++) {
        _outputRelationCollections.push_back( _inputCollections[i] + "DigiRelation" ) ;
      }
    }
    if( _inputRelationCollections.size() != _inputCollections.size()
     || _outputCollections.size() != _inputCollections.size() 
     || _outputRelationCollections.size() != _inputCollections.size() 
     || _calibrationCoefficients.empty() 
     || _calibrationCoefficients.size() != _calibrationLayers.size() ) {
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
    for (unsigned int i(0); i < _inputCollections.size(); ++i) {
      std::string colName =  _inputCollections[i] ;
      std::string relName =  _inputRelationCollections[i] ;
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
        evt->addCollection( outputCollection.release(),         _outputCollections[i] ) ;
        evt->addCollection( relationOutputCollection.release(), _outputRelationCollections[i] ) ;
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
    for (unsigned int k(0); k < _calibrationLayers.size(); ++k ) {
      if ( k > 0 ) {
        min += _calibrationLayers[k-1] ;
      }
      max += _calibrationLayers[k] ;
      if (ilayer >= min && ilayer < max) {
        calibrationCoefficient = _calibrationCoefficients[k] ;
        break ;
      }
    }
    return calibrationCoefficient ;
  }
  
}

