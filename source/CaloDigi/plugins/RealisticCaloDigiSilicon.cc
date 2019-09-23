// Calorimeter digitiser for the IDC ECAL and HCAL
// For other detectors/models SimpleCaloDigi should be used
#include <MarlinRecoMT/RealisticCaloDigi.h>

// -- marlin headers
#include <marlin/Logging.h>
#include <marlin/ProcessorApi.h>
#include <marlin/PluginManager.h>

// -- std headers
#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>

namespace marlinreco_mt {

  class RealisticCaloDigiSilicon : public RealisticCaloDigi {
  public:
    RealisticCaloDigiSilicon() ;

  protected:
     // from RealisticCaloDigi
    RealisticCaloDigi::EnergyScale getMyUnit() const ;
    float digitiseDetectorEnergy( RandomGenerator &gen, float energy ) const ;
    float convertEnergy( float energy, RealisticCaloDigi::EnergyScale inputUnit ) const ;
    
  private:
    marlin::Property<float> _ehEnergy {this, "silicon_pairEnergy",
                               "energy required to create e-h pair in silicon (in eV)", 3.6 } ;
  };

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  RealisticCaloDigiSilicon::RealisticCaloDigiSilicon() : 
    RealisticCaloDigi("RealisticCaloDigiSilicon") {
    _description = "Digitization of silicon simcalo hits" ;
  }
  
  //--------------------------------------------------------------------------
  
  RealisticCaloDigi::EnergyScale RealisticCaloDigiSilicon::getMyUnit() const {
    return RealisticCaloDigi::EnergyScale::MIP;
  }
  
  //--------------------------------------------------------------------------

  float RealisticCaloDigiSilicon::convertEnergy( float energy, RealisticCaloDigi::EnergyScale inUnit ) const {
    // converts input energy to MIP scale
    if( inUnit == RealisticCaloDigi::EnergyScale::MIP ) {
      return energy ;
    }
    else if ( inUnit == RealisticCaloDigi::EnergyScale::GEVDEP ) {
      return energy / _calib_mip;
    }
    else {
      marlin::ProcessorApi::abort( this, "convertEnergy: Unknown conversion unit!" ) ;
    }
    // just to make the compiler happy ...
    return 0.f ;
  }

  //--------------------------------------------------------------------------

  float RealisticCaloDigiSilicon::digitiseDetectorEnergy( RandomGenerator &gen, float energy ) const {
    // applies extra digitisation to silicon hits
    //  input energy in deposited GeV
    //  output is MIP scale
    float smeared_energy(energy) ;
    if ( _ehEnergy > 0 ) {
      // calculate #e-h pairs
      float nehpairs = 1e9*energy / _ehEnergy; // check units of energy! _ehEnergy is in eV, energy in GeV
      // fluctuate it by Poisson (actually an overestimate: Fano factor actually makes it smaller, however even this overstimated effect is tiny for our purposes)
      std::poisson_distribution<int> poiss( nehpairs ) ;
      smeared_energy *= poiss( gen ) ;
    }
     // convert to MIP units
    return smeared_energy / _calib_mip;
  }

  // processor declaration
  MARLIN_DECLARE_PROCESSOR( RealisticCaloDigiSilicon )
}

