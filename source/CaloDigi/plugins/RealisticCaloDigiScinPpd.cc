// Calorimeter digitiser for the IDC ECAL and HCAL
// For other detectors/models SimpleCaloDigi should be used

// base processor
#include <MarlinRecoMT/RealisticCaloDigi.h>

// -- marlin headers
#include <marlin/Logging.h>
#include <marlin/ProcessorApi.h>

// -- std headers
#include <string>
#include <algorithm>

namespace marlinreco_mt {

  class RealisticCaloDigiScinPpd : public RealisticCaloDigi {
    
   public:
    marlin::Processor *newProcessor() { return new RealisticCaloDigiScinPpd ; }
    RealisticCaloDigiScinPpd() ;

   protected:
     // from RealisticCaloDigi
    RealisticCaloDigi::EnergyScale getMyUnit() const ;
    float digitiseDetectorEnergy( RandomGenerator &gen, float energy ) const ;
    float convertEnergy( float energy, RealisticCaloDigi::EnergyScale inputUnit ) const ;

  private:
    marlin::Property<float> _PPD_pe_per_mip {this, "ppd_mipPe" ,
                               "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0" , 10. } ;

    marlin::Property<int> _PPD_n_pixels {this, "ppd_npix" ,
                               "total number of MPPC/SiPM pixels for implementation of saturation effect" , 10000 } ;

    marlin::Property<float> _misCalibNpix {this, "ppd_npix_uncert" ,
                               "fractional uncertainty of effective total number of MPPC/SiPM pixels" , 0.05 } ;

    marlin::Property<float> _pixSpread {this, "ppd_pix_spread",
                               "variation of PPD pixel signal (as a fraction: 0.01=1%)", 0.05 } ;
  };

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  RealisticCaloDigiScinPpd::RealisticCaloDigiScinPpd() : 
    RealisticCaloDigi( "RealisticCaloDigiScinPpd" ) {
    _description = "Performs digitization of sim calo hits..." ;
  }
  
  //--------------------------------------------------------------------------
  
  RealisticCaloDigi::EnergyScale RealisticCaloDigiScinPpd::getMyUnit() const {
    return RealisticCaloDigi::EnergyScale::NPE;
  }

  //--------------------------------------------------------------------------

  float RealisticCaloDigiScinPpd::convertEnergy( float energy, RealisticCaloDigi::EnergyScale inUnit ) const {
    if ( inUnit == RealisticCaloDigi::EnergyScale::NPE ) {
      return energy;
    }
    else if ( inUnit == RealisticCaloDigi::EnergyScale::MIP ) {
      return _PPD_pe_per_mip * energy ;
    }
    else if ( inUnit == RealisticCaloDigi::EnergyScale::GEVDEP ) {
      return _PPD_pe_per_mip * energy / _calib_mip ;
    }
    else {
      marlin::ProcessorApi::abort( this, "convertEnergy: Unknown conversion unit!" ) ;
    }
    // just to make the compiler happy ...
    return 0.f ;
  }
  
  //--------------------------------------------------------------------------

  float RealisticCaloDigiScinPpd::digitiseDetectorEnergy( RandomGenerator &gen, float energy ) const {
    // input energy in deposited GeV
    // output in npe
    float npe = energy*_PPD_pe_per_mip / _calib_mip; // convert to pe scale

    if ( _PPD_n_pixels > 0 ) {
      // apply average sipm saturation behaviour
      npe = _PPD_n_pixels*(1.0 - exp( -npe/_PPD_n_pixels ) ) ;
      //apply binomial smearing
      float p = npe / _PPD_n_pixels ; // fraction of hit pixels on SiPM
      std::binomial_distribution<int> binom( _PPD_n_pixels, p ) ;
      npe = binom( gen ) ; //npe now quantised to integer pixels

      if ( _pixSpread > 0) {
        // variations in pixel capacitance
        std::normal_distribution<float> norm( 1., _pixSpread / std::sqrt(npe) ) ;
        npe *= norm( gen ) ;
      }
    }
    return npe;
  }

  // processor declaration
  RealisticCaloDigiScinPpd aRealisticCaloDigiScinPpd ;
}
