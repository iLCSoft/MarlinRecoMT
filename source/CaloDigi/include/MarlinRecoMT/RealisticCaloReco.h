#ifndef MARLINRECOMT_REALISTICCALORECO_H
#define MARLINRECOMT_REALISTICCALORECO_H 1

// -- marlin headers
#include <marlin/Processor.h>

// -- lcio headers
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCEvent.h>

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
    
    virtual float reconstructEnergy( const EVENT::CalorimeterHit *hit ) const = 0 ;  // to be overloaded, technology-specific

    // processor parameters
    std::vector <std::string>        _inputCollections {} ;
    std::vector <std::string>        _inputRelationCollections {} ;
    std::vector <std::string>        _outputCollections {} ;
    std::vector <std::string>        _outputRelationCollections {} ;
    std::vector <float>              _calibrationCoefficients {} ;
    std::vector <int>                _calibrationLayers {} ;
    std::string                      _cellIDLayerString {} ;
  } ;
  
}

#endif



