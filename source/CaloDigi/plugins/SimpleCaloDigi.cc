
// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/ProcessorApi.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- marlinreco mt headers
#include <MarlinRecoMT/CalorimeterHitType.h>

// -- lcio headers
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// #include <algorithm>
// #include <string>
#include <cctype> 
#include <cstdlib>  // abs

// -- dd4hep headers
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

namespace marlinreco_mt {
  
  /** === SimpleCaloDigi Processor === <br>
   *  Simple calorimeter digitizer for calorimeter detectors.
   *  Converts SimCalorimeterHit collection to a 
   *  CalorimeterHit collection applying a threshold and an calibration constant... 
   *  Works for muon chambers, standard calorimeters and FCal calorimeters as well
   *  @version $Id$
   */
  class SimpleCaloDigi : public marlin::Processor {
  public:
    static constexpr const char *RELATIONFROMTYPESTR = "FromType" ;
    static constexpr const char *RELATIONTOTYPESTR = "ToType" ;
    
  public:
    SimpleCaloDigi() ;
    marlin::Processor *newProcessor() { return new SimpleCaloDigi() ; }
    void init() ;
    void processEvent( EVENT::LCEvent * evt ) ; 
    
  private:
    bool useLayer( unsigned int layer ) const ;
    
  protected:
    marlin::InputCollectionsProperty _inputCollections {this, EVENT::LCIO::SIMCALORIMETERHIT, "InputCollections" , 
            "Sim calo hit collection names" } ;
     
    marlin::OutputCollectionProperty _outputCollection {this, EVENT::LCIO::CALORIMETERHIT, "OutputCollection" , 
            "Calo hit output collection of real Hits" } ;
    
    marlin::OutputCollectionProperty _outputRelCollection {this, EVENT::LCIO::LCRELATION, "RelationOutputCollection" , 
            "CaloHit Relation Collection" } ;
    
    marlin::Property<float> _energyThreshold {this, "EnergyThreshold" , 
             "Threshold for sim calo hit hits in GeV (raw deposited energy, not calibrated)" , 0.f } ;

    marlin::Property<float> _calibrationCoefficient {this, "CalibrCoeff" , 
             "Calibration coefficient for calo hits" , 1.f } ;

    marlin::Property<float> _maxHitEnergy {this, "MaxHitEnergy", 
             "maximum hit energy for a calo hit" , std::numeric_limits<float>::max() } ;

    marlin::Property<std::vector<unsigned int>> _layersToKeep {this, "KeepLayers" , 
             "Vector of layers to be kept. Layers start at 1!" } ;

    marlin::Property<std::string> _cellIDLayerString {this, "CellIDLayerString" ,
             "Name of the part of the cellID that holds the layer", "layer" } ;
    
    marlin::Property<std::string> _detectorName {this, "DetectorName" ,
             "Name of the subdetector" } ;
             
    marlin::Property<std::string> _caloType {this, "CaloType" ,
            "type of calorimeter: em, had, muon" } ;

    marlin::Property<std::string> _caloID {this, "CaloID" ,
            "ID of calorimeter: lcal, fcal, bcal" } ;

    marlin::Property<std::string> _caloLayout {this, "CaloLayout" ,
            "subdetector layout: barrel, endcap, plug, ring" } ;
            
    std::vector<bool>            _useLayers {} ;    
  };
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  SimpleCaloDigi::SimpleCaloDigi() : 
    marlin::Processor("SimpleCaloDigi") {
    _description = "Performs simple digitization of sim hits..." ;
  }
  
  //--------------------------------------------------------------------------

  void SimpleCaloDigi::init() {
    printParameters() ;
    // Get the number of Layers in detector
    unsigned int nLayers = 0 ;
    try {
      dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance() ;
      dd4hep::DetElement theDetector = mainDetector.detector( _detectorName ) ;
      auto calorimeterParameters  = theDetector.extension<dd4hep::rec::LayeredCalorimeterData>() ;
      nLayers =  calorimeterParameters->layers.size() ;
    }
    catch( std::exception& e ) {
      marlin::ProcessorApi::abort( this, "No detector available: " + std::string(e.what()) ) ;
    }
    // If the vectors are empty, we are keeping everything 
    if(_layersToKeep.get().size() > 0) {
      // Layers start at 0
      for(unsigned int i = 0; i < nLayers; ++i) {
        _useLayers.push_back(false) ;
        for(auto iter = _layersToKeep.get().begin(); iter < _layersToKeep.get().end(); ++iter) {
        	if (i == *iter-1) {
        	  _useLayers[i] = true ; 
            break;
        	}
        }
      }
    }
  }
  
  //--------------------------------------------------------------------------

  void SimpleCaloDigi::processEvent( EVENT::LCEvent * evt ) { 
    auto outputCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::CALORIMETERHIT ) ;
    auto relationCollection  = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::LCRELATION ) ;
    relationCollection->parameters().setValue( RELATIONFROMTYPESTR , EVENT::LCIO::CALORIMETERHIT ) ;
    relationCollection->parameters().setValue( RELATIONTOTYPESTR   , EVENT::LCIO::SIMCALORIMETERHIT ) ;
    IMPL::LCFlagImpl flag ;
    flag.setBit( EVENT::LCIO::CHBIT_LONG ) ;
    flag.setBit( EVENT::LCIO::CHBIT_ID1 ) ;
    outputCollection->setFlag( flag.getFlag() ) ;
    std::string initString ;
    // layout information
    auto caloLayout = layoutFromString( _caloLayout ) ; 
    auto caloID = caloIDFromString( _caloID ) ; 
    auto caloType = caloTypeFromString( _caloType ) ;
    // loop over input collections
    for (unsigned int i(0); i < _inputCollections.get().size(); ++i) {
      std::string colName =  _inputCollections.get()[i] ;
      try {
        auto collection = evt->getCollection( colName ) ;
        initString = collection->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ) ;
        int numElements = collection->getNumberOfElements() ;
        UTIL::CellIDDecoder<EVENT::SimCalorimeterHit> idDecoder( collection ) ;
        log<DEBUG3>() << "Number of hits: " << numElements << std::endl ;
        // Loop over hits in the current collection
        for (int j(0); j < numElements; ++j) {
        	auto hit = static_cast<EVENT::SimCalorimeterHit*>( collection->getElementAt( j ) ) ;
          if( nullptr == hit ) {
            continue ;
          }
        	float energy = hit->getEnergy() ;
        	int cellid = hit->getCellID0() ;
        	int cellid1 = hit->getCellID1() ;
        	unsigned int layer = std::abs( idDecoder(hit)[ _cellIDLayerString ] ) ;
        	//Check if we want to use this layer, else go to the next hit
        	if( not useLayer( layer ) ) {
            log<DEBUG3>() << "  Skipping hit '" << hit->id() << "' in layer " << layer << std::endl ;
            continue ;
          }
        	float calibratedEnergy = _calibrationCoefficient * energy ;
        	if( calibratedEnergy > _maxHitEnergy ) {
            calibratedEnergy = _maxHitEnergy ;
          }
        	if ( energy > _energyThreshold ) {
            log<DEBUG3>() << "  Accepting hit " << hit->id() << std::endl ;
        	  auto calhit = new IMPL::CalorimeterHitImpl();
        	  calhit->setCellID0( cellid ) ;
        	  calhit->setCellID1( cellid1 ) ;
        	  calhit->setEnergy( calibratedEnergy ) ;
        	  calhit->setPosition( hit->getPosition() ) ;
        	  calhit->setType( CHT( caloType, caloID, caloLayout, layer ) );
        	  calhit->setRawHit( hit ) ;
        	  outputCollection->addElement( calhit ) ;
            // create a calo hit <-> sim calo hit relation
        	  auto rel = new IMPL::LCRelationImpl( calhit, hit, 1. ) ;
        	  relationCollection->addElement( rel ) ;
        	}
        }
      }
      catch(EVENT::DataNotAvailableException &e) {
        log<WARNING>() << "Collection " << colName << " not available: " << e.what() << std::endl ;
      }
    }
    outputCollection->parameters().setValue( EVENT::LCIO::CellIDEncoding, initString ) ;
    evt->addCollection( outputCollection.release(), _outputCollection ) ;
    evt->addCollection( relationCollection.release(), _outputRelCollection ) ;
  }
  
  //--------------------------------------------------------------------------

  bool SimpleCaloDigi::useLayer( unsigned int layer ) const {
    if( layer > _useLayers.size() || _useLayers.size() == 0 ) {
      return true ;
    }
    return _useLayers[layer] ;
  }
  
  // processor declaration
  SimpleCaloDigi aSimpleCaloDigi ;
}
