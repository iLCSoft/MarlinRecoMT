#include <MarlinRecoMT/RealisticCaloReco.h>

#include <marlin/PluginManager.h>

namespace marlinreco_mt {
  
  /** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of silicon calorimeter hits
    D.Jeans 02/2016.
    24 March 2016: removed gap corrections - to be put into separate processor
  */
  class RealisticCaloRecoSilicon : public RealisticCaloReco {
  public:
    RealisticCaloRecoSilicon() ;

  private:
    float reconstructEnergy( UTIL::CellIDDecoder<EVENT::CalorimeterHit> &decoder, const EVENT::CalorimeterHit *hit ) const override ;
  };

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  RealisticCaloRecoSilicon::RealisticCaloRecoSilicon() : 
    RealisticCaloReco("RealisticCaloRecoSilicon") {
    _description = "Performs fist reconstruction of silicon ECAL hits" ;
  }
  
  //--------------------------------------------------------------------------
  
  float RealisticCaloRecoSilicon::reconstructEnergy( UTIL::CellIDDecoder<EVENT::CalorimeterHit> &decoder, const EVENT::CalorimeterHit *hit ) const {
    // here the input energy should be in MIPs
    float energy = hit->getEnergy() ;
    // what layer is this hit in?
    int layer   = decoder( hit ) [ _cellIDLayerString ] ;
    // now correct for sampling fraction
    energy *= this->getLayerCalib( layer ) ;
    return energy ;
  }

  // processor declaration
  MARLIN_DECLARE_PROCESSOR( RealisticCaloRecoSilicon )
}

