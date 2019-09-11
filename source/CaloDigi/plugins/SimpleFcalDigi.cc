// -- marlin reco mtheaders
#include <MarlinRecoMT/CalorimeterHitType.h>

// -- marlin headers
#include <marlin/Logging.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorApi.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- lcio headers
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// -- std headers
#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>

namespace marlinreco_mt {
  
  /** === SimpleFCalDigi Processor === <br>
   *  Simple calorimeter digitizer for the Fcal Processor. <br>
   *  Converts SimCalorimeterHit collection to one 
   *  CalorimeterHit collection applying a threshold and an calibration constant...
   */
  class SimpleFCalDigi : public marlin::Processor {
  public:
    static constexpr const char *RELATIONFROMTYPESTR = "FromType" ;
    static constexpr const char *RELATIONTOTYPESTR = "ToType" ;
    marlin::Processor *newProcessor() { return new SimpleFCalDigi() ; }
    SimpleFCalDigi() ;
    
    // marlin::Processor
    void processEvent( EVENT::LCEvent * evt ) ;
    
  private:
    std::string              _inputCollection {} ;
    std::string              _outputCollection {} ;
    std::string              _outputRelCollection {} ;
    std::string              _cellIDLayerString {"K-1"} ;
    std::string              _caloLayout {"endcap"} ;
    std::string              _caloID {"fcal"} ;
    std::string              _caloType {"had"} ;
    float                    _thresholdFcal {0.f} ;
    float                    _calibrCoeffFcal {31.f} ;
  };

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  SimpleFCalDigi::SimpleFCalDigi() : 
    marlin::Processor("SimpleFCalDigi") {
      
    _description = "Performs simple digitization of SimCalorimeter hits in forward calorimeters ..." ;
    
    registerInputCollection( EVENT::LCIO::SIMCALORIMETERHIT, 
                            "FCALCollection" , 
                            "Fcal Collection Name" ,
                            _inputCollection ,
                            _inputCollection ) ;
    
    registerOutputCollection( EVENT::LCIO::CALORIMETERHIT, 
                              "FCALOutputCollection" , 
                              "Fcal Collection of real Hits" , 
                              _outputCollection , 
                              _outputCollection ) ; 
    
    registerOutputCollection( EVENT::LCIO::LCRELATION, 
                              "RelationOutputCollection" , 
                              "CaloHit Relation Collection" , 
                              _outputRelCollection , 
                              _outputRelCollection ) ; 
    
    registerProcessorParameter("FcalThreshold" , 
                               "Threshold for Fcal Hits in GeV" ,
                               _thresholdFcal,
                               _thresholdFcal ) ;
    
    registerProcessorParameter("CalibrFCAL" , 
                             "Calibration coefficients for FCAL" ,
                             _calibrCoeffFcal,
                             _calibrCoeffFcal ) ;

    registerProcessorParameter("CellIDLayerString" ,
                               "name of the part of the cellID that holds the layer" , 
                               _cellIDLayerString , 
                               _cellIDLayerString ) ;
                                                            
    registerProcessorParameter("CaloType" ,
                               "type of calorimeter: em, had, muon" , 
                               _caloType , 
                               _caloType ) ;

    registerProcessorParameter("CaloID" ,
                               "ID of calorimeter: lcal, fcal, bcal", 
                               _caloID , 
                               _caloID ) ;

    registerProcessorParameter("CaloLayout" ,
                               "subdetector layout: barrel, endcap, plug, ring", 
                               _caloLayout , 
                               _caloLayout ) ;
  }
  
  //--------------------------------------------------------------------------
  
  void SimpleFCalDigi::processEvent( EVENT::LCEvent *evt ) {
    try {
      // get input collection and parameters
      auto inputCollection = evt->getCollection( _inputCollection ) ;
      auto cellIDString = inputCollection->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ) ;
      int numElements = inputCollection->getNumberOfElements() ;
      UTIL::CellIDDecoder<EVENT::SimCalorimeterHit> cellIDDecoder( inputCollection ) ;
      // create output collections
      auto outputCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::CALORIMETERHIT ) ;
      outputCollection->parameters().setValue( EVENT::LCIO::CellIDEncoding, cellIDString ) ;
      auto outputRelationCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::LCRELATION ) ;
      outputRelationCollection->parameters().setValue( RELATIONFROMTYPESTR , EVENT::LCIO::CALORIMETERHIT ) ;
      outputRelationCollection->parameters().setValue( RELATIONTOTYPESTR   , EVENT::LCIO::SIMCALORIMETERHIT ) ;
      IMPL::LCFlagImpl flag {} ;
      flag.setBit( EVENT::LCIO::CHBIT_LONG ) ;
      outputCollection->setFlag( flag.getFlag() ) ;
      
      for ( int j=0 ; j<numElements; ++j ) {
        auto hit = static_cast<EVENT::SimCalorimeterHit*>( inputCollection->getElementAt( j ) ) ;
        if ( hit->getEnergy() > _thresholdFcal ) {
          // create a digitized hit
          auto calhit = new IMPL::CalorimeterHitImpl() ;
          calhit->setCellID0( hit->getCellID0() ) ;
          calhit->setCellID1( hit->getCellID1() ) ;
          calhit->setEnergy( _calibrCoeffFcal * hit->getEnergy() ) ;
          calhit->setPosition( hit->getPosition() ) ;
          calhit->setType( CHT ( caloTypeFromString(_caloType), caloIDFromString(_caloID), layoutFromString(_caloLayout),  cellIDDecoder( hit )[ _cellIDLayerString ] ) );
          calhit->setRawHit( hit ) ;
          outputCollection->addElement( calhit ) ;
          /// create a relation object
          auto rel = new IMPL::LCRelationImpl( calhit, hit, 1.) ;
          outputRelationCollection->addElement( rel ) ;
        }
      }
      evt->addCollection( outputCollection.release(), _outputCollection ) ;
      evt->addCollection( outputRelationCollection.release(), _outputRelCollection );
    }
    catch(EVENT::DataNotAvailableException &e){
      log<DEBUG3>() << "FCal input collection " << _inputCollection << " not available: " << e.what() << std::endl ;
    }
  }

  // processor declaration
  SimpleFCalDigi aSimpleFCalDigi ;
}

