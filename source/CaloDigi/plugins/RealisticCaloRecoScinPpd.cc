#include <MarlinRecoMT/RealisticCaloReco.h>

#include <marlin/PluginManager.h>

namespace marlinreco_mt {
  
  /** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of scint+PPD calorimeter hits
    D.Jeans 02/2016.
  */
  class RealisticCaloRecoScinPpd : public RealisticCaloReco {
  public:
    RealisticCaloRecoScinPpd() ;

  private:
    float reconstructEnergy( UTIL::CellIDDecoder<EVENT::CalorimeterHit> &decoder, const EVENT::CalorimeterHit *hit ) const override ;

  private:
    marlin::Property<float> _photoelectronsPerMIP {this, "ppd_mipPe",
                            "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0", 10.f } ;
                            
    marlin::Property<int> _nPixels {this, "ppd_npix",
                            "total number of MPPC/SiPM pixels for implementation of saturation effect", 10000 } ;
  };

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  RealisticCaloRecoScinPpd::RealisticCaloRecoScinPpd() : 
    RealisticCaloReco("RealisticCaloRecoScinPpd") {
    _description = "Performs fist reconstruction of scintillator calo hits" ;
  }
  
  //--------------------------------------------------------------------------
  
  float RealisticCaloRecoScinPpd::reconstructEnergy( UTIL::CellIDDecoder<EVENT::CalorimeterHit> &decoder, const EVENT::CalorimeterHit *hit ) const {
    // here the input energy should be in NPE
    float energy = hit->getEnergy() ;
    // first de-saturate PPD response
    // this is the fraction of SiPM pixels fired above which a linear continuation of the saturation-reconstruction function is used. 
    // 0.95 of nPixel corresponds to a energy correction of factor ~3.
    const float r = 0.95 ;
    if (energy < r*_nPixels) { //current hit below linearisation threshold, reconstruct energy normally:
      energy = -_nPixels * std::log ( 1. - ( energy / _nPixels ) ) ;
    } 
    else { //current hit is aove linearisation threshold, reconstruct using linear continuation function:
      energy = 1 / ( 1 - r ) * ( energy - r*_nPixels ) - _nPixels * std::log( 1 - r ) ;
    }
    // then go back to MIP scale
    energy /= _photoelectronsPerMIP ;
    // what layer is this hit in?
    int layer   = decoder(hit)[_cellIDLayerString];
    // now correct for sampling fraction (calibration from MIP -> shower GeV)
    energy *= this->getLayerCalib( layer ) ;

    return energy ;
  }

  // processor declaration
  MARLIN_DECLARE_PROCESSOR( RealisticCaloRecoScinPpd )
}

