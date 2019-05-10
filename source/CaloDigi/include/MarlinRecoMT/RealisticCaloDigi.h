#ifndef MARLINRECOMT_REALISTICCALODIGI_H
#define MARLINRECOMT_REALISTICCALODIGI_H 1

// -- marlin headers
#include <marlin/Processor.h>

// -- lcio headers
#include <IMPL/LCFlagImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCEvent.h>

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
     *  @param  energy the input energy
     */
    virtual float digitiseDetectorEnergy( float energy ) const = 0 ;
    
    /**
     *  @brief  Convert the input energy with the specified scale
     *  
     *  @param  energy the input energy
     *  @param  inScale the energy scale
     */
    virtual float convertEnergy( float energy, EnergyScale inScale ) const = 0 ;

  protected:
    // processor parameters
    ///< The input collection list
    std::vector<std::string> _inputCollections {} ;
    ///< The output collection list
    std::vector<std::string> _outputCollections {} ;
    ///< The output relation collection list
    std::vector<std::string> _outputRelCollections {} ;
    ///< The hit energy threshold
    float _threshold_value {} ;
    ///< The hit energy threshold unit
    std::string _threshold_unit {} ;
    ///< apply timing cuts?
    int   _time_apply {} ; 
    ///< correct times for propagation?             
    int   _time_correctForPropagation {} ; 
    ///< timing window minimum
    float _time_windowMin {} ;  
    ///< timing window minimum        
    float _time_windowMax{ } ;
    ///< MIP calibration factor (most probable energy deposit by MIP in active material of one layer)
    float _calib_mip {} ;              
    ///< general miscalibration (uncorrelated between channels)
    float _misCalib_uncorrel {} ;       
    ///< if true, use the same cell miscalibs over events (requires more memory)
    bool  _misCalib_uncorrel_keep {} ;  
    ///< general miscalibration (100% uncorrelated between channels)
    float _misCalib_correl {} ;           
    ///< fraction of random dead channels
    float _deadCell_fraction {} ;       
    ///< keep same cells dead between events? (requires more memory)
    bool  _deadCell_keep {} ;           
    ///< electronics noise (as fraction of MIP)
    float _elec_noiseMip {} ;
    ///< electronics dynamic range (in terms of MIPs)
    float _elec_rangeMip {} ;           
    ///< The cell id string for layers
    std::string _cellIDLayerString {} ;
    
    // internal variables (const by usage in processEvent)
    EnergyScale _threshold_iunit {} ;
    IMPL::LCFlagImpl _flag {} ;
    IMPL::LCFlagImpl _flag_rel {} ;
    
    
    float _event_correl_miscalib {} ;
    std::map < std::pair <int, int> , float > _cell_miscalibs {} ;
    std::map < std::pair <int, int> , bool  > _cell_dead {} ;
  };
  
}

#endif



