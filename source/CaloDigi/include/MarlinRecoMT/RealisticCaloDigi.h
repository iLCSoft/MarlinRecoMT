#ifndef MARLINRECOMT_REALISTICCALODIGI_H
#define MARLINRECOMT_REALISTICCALODIGI_H 1

// -- marlin headers
#include <marlin/Processor.h>

// -- lcio headers
#include <IMPL/LCFlagImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>

// -- std headers
#include <string>
#include <vector>
#include <random>
#include <map> // for pair

namespace marlinreco_mt {

  /** === RealisticCaloDigi Processor === <br>
      Digitisation of calorimeter hits
      e.g. timing, dead cells, miscalibrations
      this is virtual class, technology-blind
      technology-specific classes can inherit from this one.
      D. Jeans 02/2016, rewrite of parts of ILDCaloDigi, DDCaloDigi
      R. Ete 05/2019, rewrite for MT use
   */
  class RealisticCaloDigi : public marlin::Processor {
  public:
    using RandomGenerator = std::mt19937 ;
    static constexpr const char *RELATIONFROMTYPESTR = "FromType" ;
    static constexpr const char *RELATIONTOTYPESTR = "ToType" ;
    
  public:
    virtual ~RealisticCaloDigi() = default ;
    RealisticCaloDigi ( const RealisticCaloDigi& ) = delete;
    RealisticCaloDigi& operator=(const RealisticCaloDigi&) = delete;

    /**
     *  @brief  Constructor
     * 
     *  @param  pname the processor name (implementation)
     */
    RealisticCaloDigi( const std::string &pname ) ;

    /**
     *  @brief  Initialize base parameters
     */
    virtual void init() ;
    
    /**
     *  @brief  Process an event
     * 
     *  @param  evt the event to process
     */
    void processEvent( EVENT::LCEvent *evt ) ; 

   protected:
    /**
     *  @brief  EnergyScale enumerator
     */
    enum class EnergyScale : unsigned int {
      MIP,        /// Energy deposit in MIP
      GEVDEP,     /// Energy deposit in Gev
      NPE         /// Number of photo-electrons
    };
    
    /**
     * 
     */
    struct EventData {
      RandomGenerator          _generator {} ;
      float                    _eventCorrelMiscalib {} ;
    };

    /**
     *  @brief  From inout energy, returns the digitized energy with correction factors applied
     *
     *  @param  evtdata the additional event data 
     *  @param  energy the input sim hit energy
     */
    float energyDigi( EventData &evtdata, float energy ) const ;
    
    /**
     *  @brief  Apply timing cuts on the sim hit
     * 
     *  @param  hit the input sim hit
     */
    std::vector < std::pair < float , float > > applyTimingCuts( const EVENT::SimCalorimeterHit * hit ) const ;

    /**
     *  @brief  Get the energy unit
     */
    virtual EnergyScale getMyUnit() const = 0 ;
    
    /**
     *  @brief  Digitize the detector energy
     *
     *  @param  gen the random number generator to use
     *  @param  energy the input energy
     */
    virtual float digitiseDetectorEnergy( RandomGenerator &gen, float energy ) const = 0 ;
    
    /**
     *  @brief  Convert the input energy with the specified scale
     *  
     *  @param  energy the input energy
     *  @param  inScale the energy scale
     */
    virtual float convertEnergy( float energy, EnergyScale inScale ) const = 0 ;

  protected:
        
    marlin::InputCollectionsProperty _inputCollections {this, EVENT::LCIO::SIMCALORIMETERHIT, "inputHitCollections" ,
                            "Input simcalhit Collection Names" , {"SimCalorimeterHits"} } ;

    marlin::Property<EVENT::StringVec> _outputCollections {this, "outputHitCollections",
  			                    "Output calorimeterhit Collection Names" } ;

    marlin::Property<EVENT::StringVec> _outputRelCollections {this, "outputRelationCollections",
  			                    "Output hit relation Collection Names" } ;

    marlin::Property<float> _threshold_value {this, "threshold",
  			                    "Threshold for Hit", 0.5 } ;
                                                        
    marlin::Property<std::string> _threshold_unit {this, "thresholdUnit",
                            "Unit for threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants", "MIP" } ;
                            
    marlin::Property<bool> _time_apply {this, "timingCut",
                            "Use hit times", false } ;
                                                        
    marlin::Property<bool> _time_correctForPropagation {this, "timingCorrectForPropagation",
                            "Correct hit times for propagation: radial distance/c", false } ;
                            
    marlin::Property<float> _time_windowMin {this, "timingWindowMin",
                            "Time Window minimum time in ns", -10. } ;
    
    marlin::Property<float> _time_windowMax {this, "timingWindowMax",
                            "Time Window maximum time in ns", 100. } ;

    marlin::Property<float> _calib_mip {this, "calibration_mip",
                            "average G4 deposited energy by MIP for calibration", 1.e-4 } ;
                            
    marlin::Property<float> _misCalib_uncorrel {this, "miscalibration_uncorrel",
                            "uncorrelated random gaussian miscalibration (as a fraction: 1.0 = 100%)", 0. } ;

    marlin::Property<float> _misCalib_correl {this, "miscalibration_correl",
                            "correlated random gaussian miscalibration (as a fraction: 1.0 = 100%)", 0. } ;

    marlin::Property<float> _deadCell_fraction {this, "deadCell_fraction",
                            "random dead cell fraction (as a fraction: 0->1)", 0. } ;
                            
    marlin::Property<float> _elec_noiseMip {this, "elec_noise_mip",
                            "typical electronics noise (in MIP units)", 0. } ;
                            
    marlin::Property<float> _elec_rangeMip {this, "elec_range_mip",
                            "maximum of dynamic range of electronics (in MIPs)", 2500 } ;
                            
    marlin::Property<std::string> _cellIDLayerString {this, "CellIDLayerString",
                            "name of the part of the cellID that holds the layer", "K-1" } ;
    
    // internal variables (const by usage in processEvent)
    EnergyScale _threshold_iunit {} ;
    IMPL::LCFlagImpl _flag {} ;
    IMPL::LCFlagImpl _flag_rel {} ;
  };
  
}

#endif



