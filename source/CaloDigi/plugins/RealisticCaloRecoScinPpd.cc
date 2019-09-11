#include <MarlinRecoMT/RealisticCaloReco.h>

namespace marlinreco_mt {
  
  /** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of scint+PPD calorimeter hits
    D.Jeans 02/2016.
  */
  class RealisticCaloRecoScinPpd : public RealisticCaloReco {
  public:
    marlin::Processor *newProcessor() { return new RealisticCaloRecoScinPpd() ; }
    RealisticCaloRecoScinPpd() ;

  private:
    float reconstructEnergy( UTIL::CellIDDecoder<EVENT::CalorimeterHit> &decoder, const EVENT::CalorimeterHit *hit ) const override ;

  private:
    // # photoelectrons/MIP for MPPC
    float         _photoelectronsPerMIP {10.f} ;
    // # pixels in MPPC
    int           _nPixels {10000} ;
  };

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  RealisticCaloRecoScinPpd::RealisticCaloRecoScinPpd() : 
    RealisticCaloReco("RealisticCaloRecoScinPpd") {

    _description = "Performs fist reconstruction of scintillator calo hits";

    registerProcessorParameter("ppd_mipPe" ,
                               "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0" ,
                               _photoelectronsPerMIP,
                               _photoelectronsPerMIP) ;

    registerProcessorParameter("ppd_npix" ,
                               "total number of MPPC/SiPM pixels for implementation of saturation effect" ,
                               _nPixels,
                               _nPixels) ;
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
  RealisticCaloRecoScinPpd aRealisticCaloRecoScinPpd ;
}

