#ifndef MARLINRECOMT_REALISTICCALORECO_H
#define MARLINRECOMT_REALISTICCALORECO_H 1

// -- marlin headers
#include <marlin/Processor.h>

// -- lcio headers
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCEvent.h>
#include <UTIL/CellIDDecoder.h>

#include <string>
#include <vector>

namespace marlinreco_mt {

  /** === RealisticCaloReco Processor === <br>
      realistic reconstruction of calorimeter hits
      e.g. apply sampling fraction correction
      virtual class, technology indenpendent
      D.Jeans 02/2016.

      24 March 2016: removed gap corrections - to be put into separate processor
      changed relations: now keep relation between reconstructed and simulated hits.
  */
  class RealisticCaloReco : public marlin::Processor {
  public:
    static constexpr const char *RELATIONFROMTYPESTR = "FromType" ;
    static constexpr const char *RELATIONTOTYPESTR = "ToType" ;
    
  public:
    RealisticCaloReco( const std::string &pname ) ;
    RealisticCaloReco ( const RealisticCaloReco& ) = delete;
    RealisticCaloReco& operator=(const RealisticCaloReco&) = delete;
    virtual ~RealisticCaloReco() = default ;

    // from marlin::Processor
    virtual void init() ;
    virtual void processEvent( EVENT::LCEvent * evt ) ;

   protected:
    float getLayerCalib( int ilayer ) const ;
    // to be overloaded, technology-specific
    virtual float reconstructEnergy( UTIL::CellIDDecoder<EVENT::CalorimeterHit> &decoder, const EVENT::CalorimeterHit *hit ) const = 0 ;


  protected:
    marlin::InputCollectionsProperty _inputCollections {this, EVENT::LCIO::CALORIMETERHIT, "inputHitCollections",
                              "input hit collection names" } ;

    marlin::InputCollectionsProperty _inputRelationCollections {this, EVENT::LCIO::LCRELATION, "inputRelationCollections",
                              "input relation collection names (digi<->sim), one per inputHitCollection" } ;

    marlin::Property<std::vector<std::string>> _outputCollections{this, "outputHitCollections",
                              "output hit collection names" } ;

    marlin::Property<std::vector<std::string>> _outputRelationCollections {this, "outputRelationCollections",
                              "output hit collection names" } ;

    marlin::Property<std::vector<float>> _calibrationLayers {this, "calibration_layergroups" ,
                               "grouping of calo layers" } ;

    marlin::Property<std::vector<float>> _calibrationCoefficients {this, "calibration_factorsMipGev" ,
                               "Calibration coefficients (MIP->shower GeV) of layers groups" } ;

    marlin::Property<std::string> _cellIDLayerString {this, "CellIDLayerString" ,
                               "name of the part of the cellID that holds the layer", "K-1" } ;
  };
  
}

#endif



