// Calorimeter digitiser
#include <MarlinRecoMT/RealisticCaloDigi.h>

// -- marlin headers
#include <marlin/ProcessorApi.h>

// -- lcio headers
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// -- marlinrecomt headers
#include <MarlinRecoMT/CalorimeterHitType.h>

// -- std headers
#include <iostream>
#include <string>
#include <assert.h>
#include <cmath>

namespace marlinreco_mt {

  RealisticCaloDigi::RealisticCaloDigi( const std::string &pname ) : Processor( pname ) {

    _description = "Performs digitization of sim calo hits. Virtual class." ;
  }
  
  //--------------------------------------------------------------------------

  void RealisticCaloDigi::init() {
    // usually a good idea to
    printParameters() ;
    // check that number of input and output collections names are the same
    if( _outputCollections.get().size() != _inputCollections.get().size() ) {
      marlin::ProcessorApi::abort( this, "Input/output collection list sizes are different" ) ;
    }
    if( _outputRelCollections.get().size() != _inputCollections.get().size() ) {
      marlin::ProcessorApi::abort( this, "Input/output collection list sizes are different" ) ;
    }
    // unit in which threshold is specified
    if (_threshold_unit.get().compare("MIP") == 0) {
      _threshold_iunit = EnergyScale::MIP ;
    } 
    else if (_threshold_unit.get().compare("GeV") == 0) {
      _threshold_iunit = EnergyScale::GEVDEP ;
    } 
    else if (_threshold_unit.get().compare("px") == 0) {
      _threshold_iunit = EnergyScale::NPE ;
    } 
    else {
      marlin::ProcessorApi::abort( this, "Could not identify threshold unit. Please use \"GeV\", \"MIP\" or \"px\"!" ) ;
    }
    // convert the threshold to the approriate units (i.e. MIP for silicon, NPE for scint)
    _threshold_value = convertEnergy( _threshold_value, _threshold_iunit ) ;
    // setup output collection flags
    _flag.setBit( EVENT::LCIO::CHBIT_LONG ) ;
    _flag.setBit( EVENT::LCIO::RCHBIT_TIME ) ; //store timing on output hits.
    _flag_rel.setBit( EVENT::LCIO::LCREL_WEIGHTED ) ; // for the hit relations
    // register for random seed usage
    marlin::ProcessorApi::registerForRandomSeeds( this ) ;
  }
  
  //--------------------------------------------------------------------------

  void RealisticCaloDigi::processEvent( EVENT::LCEvent * evt ) {
    // deal with random numbers there
    auto randomSeed = marlin::ProcessorApi::getRandomSeed( this, evt ) ;
    EventData eventData ;
    eventData._generator.seed( randomSeed ) ;
    std::normal_distribution<float> miscalDistCorel( 1.0, _misCalib_correl );
    // decide on this event's correlated miscalibration
    if ( _misCalib_correl > 0 ) {
      eventData._eventCorrelMiscalib = miscalDistCorel( eventData._generator ) ;
    }
    // loop over simulated hit collections
    for ( unsigned int i=0 ; i<_inputCollections.get().size() ; ++i ) {
      auto colName = _inputCollections.get().at( i ) ;
      log<marlin::DEBUG1>() << "Looking for collection: " << colName << std::endl ;
      try {
        EVENT::LCCollection * col = evt->getCollection( colName.c_str() ) ;
        std::string initString = col->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ) ;
        CHT::CaloType cht_type = caloTypeFromString( colName ) ;
        CHT::CaloID   cht_id   = caloIDFromString( colName ) ;
        CHT::Layout   cht_lay  = layoutFromString( colName ) ;
        UTIL::CellIDDecoder<EVENT::SimCalorimeterHit> idDecoder( col );
        const auto numElements = col->getNumberOfElements();
        log<marlin::DEBUG1>() << colName << " number of elements = " << numElements << std::endl ;
        // don't go further if no hits
        if ( numElements==0 ) {
          continue ;
        }
        // create new collection: hits
        IMPL::LCCollectionVec *newcol = new IMPL::LCCollectionVec( EVENT::LCIO::CALORIMETERHIT );
        newcol->setFlag(_flag.getFlag()) ;
        // hit relations to simhits [calo -> sim]
        IMPL::LCCollectionVec *relcol  = new IMPL::LCCollectionVec( EVENT::LCIO::LCRELATION );
        relcol->setFlag(_flag_rel.getFlag());
        relcol->parameters().setValue( RELATIONFROMTYPESTR , EVENT::LCIO::CALORIMETERHIT ) ;
        relcol->parameters().setValue( RELATIONTOTYPESTR   , EVENT::LCIO::SIMCALORIMETERHIT ) ;
        // loop over input hits
        for ( int j=0 ; j<numElements ; ++j ) {
          EVENT::SimCalorimeterHit * simhit = dynamic_cast<EVENT::SimCalorimeterHit*>( col->getElementAt( j ) ) ;
          // deal with timing aspects
          std::vector<std::pair<float,float>> timeClusteredHits ; // vector of (time, energy)
          if( _time_apply ) {
            timeClusteredHits = applyTimingCuts( simhit ) ;
          } 
          else { // just take full energy, assign to time 0
            timeClusteredHits.push_back( std::pair<float,float>( 0, simhit->getEnergy() ) );
          }
          // loop over all hits
          for ( std::size_t jj=0 ; jj<timeClusteredHits.size() ; jj++ ) {
            float hittime   = timeClusteredHits[jj].first ;
            float energyDep = timeClusteredHits[jj].second ;
            // apply extra energy digitisation onto the energy
            float energyDig = energyDigi( eventData, energyDep ) ;

  	        log<marlin::DEBUG0>() << " hit " << jj << " time: " << hittime << " eDep: " << energyDep << " eDigi: " << energyDig << " " << _threshold_value << std::endl ;

            if ( energyDig > _threshold_value ) { // write out this hit
              IMPL::CalorimeterHitImpl* newhit = new IMPL::CalorimeterHitImpl() ;
              newhit->setCellID0( simhit->getCellID0() ) ;
              newhit->setCellID1( simhit->getCellID1() ) ;
              newhit->setTime( hittime ) ;
              newhit->setPosition( simhit->getPosition() ) ;
  	          newhit->setEnergy( energyDig ) ;
  	          int layer = idDecoder( simhit )[ _cellIDLayerString ] ;
  	          newhit->setType( CHT( cht_type, cht_id, cht_lay, layer ) ) ;
              newhit->setRawHit( simhit ) ;
              newcol->addElement( newhit ) ; // add hit to output collection
  	          log<marlin::DEBUG1>() << "orig/new hit energy: " << simhit->getEnergy() << " " << newhit->getEnergy() << std::endl ;
              // add a relation reco <-> sim
              IMPL::LCRelationImpl *rel = new IMPL::LCRelationImpl( newhit, simhit, 1.0 ) ;
              relcol->addElement( rel ) ;
            } // threshold
          } // time sliced hits
        } // input hits
        // add collection to event
        newcol->parameters().setValue( EVENT::LCIO::CellIDEncoding, initString );
        evt->addCollection( newcol, _outputCollections.get()[i] );
        // add relation collection to event
        evt->addCollection( relcol, _outputRelCollections.get()[i] );
      } 
      catch(EVENT::DataNotAvailableException &e) {
        log<marlin::DEBUG1>() << "Could not find input collection " << colName << std::endl;
      }
    }
    log<marlin::MESSAGE>() << "End of event " << evt->getEventNumber() << std::endl ;
  }

  //--------------------------------------------------------------------------

  std::vector<std::pair<float,float>> RealisticCaloDigi::applyTimingCuts( const EVENT::SimCalorimeterHit * hit ) const {
    // apply timing cuts on simhit contributions
    //  outputs a vector of (time,energy) pairs
    //  for now, only get one output hit per input hit, however we keep the possibility to have more
    std::vector<std::pair<float,float>> timedhits ;
    float timeCorrection(0);
    if ( _time_correctForPropagation ) { // time of flight from IP to this point
      float r(0);
      for (int i=0; i<3; i++) {
        r += pow( hit->getPosition()[i], 2 ) ;
      }
      timeCorrection = std::sqrt(r) / 299.79; // [speed of light in mm/ns]
    }
    // this is Oskar's simple (and probably the most correct) method for treatment of timing
    //  - collect energy in some predefined time window around collision time (possibly corrected for TOF)
    //  - assign time of earliest contribution to hit
    float energySum = 0.f ;
    float earliestTime = std::numeric_limits<float>::max() ;
    for(int i = 0; i<hit->getNMCContributions(); i++){ // loop over all contributions
      float timei   = hit->getTimeCont(i) ; //absolute hit timing of current subhit
      float energyi = hit->getEnergyCont(i) ; //energy of current subhit
      float relativetime = timei - timeCorrection; // wrt time of flight
      if (relativetime > _time_windowMin && relativetime < _time_windowMax ) {
        energySum += energyi ;
        if (relativetime < earliestTime){
  	      earliestTime = relativetime; //use earliest hit time for simpletimingcut
        }
      }
    }      
    //accept this hit ?
    if( earliestTime > _time_windowMin && earliestTime < _time_windowMax ) { 
      timedhits.push_back( std::pair<float,float> ( earliestTime, energySum ) ) ;
    }
    return timedhits ;
  }

  //--------------------------------------------------------------------------

  float RealisticCaloDigi::energyDigi( EventData &evtdata, float energy ) const {
    // some extra digi effects
    // controlled by _applyDigi = 0 (none), 1 (apply)
    // input parameters: hit energy ( in any unit: effects are all relative )
    //                   id0,1 - cell IDs (used to memorise miscalibrations/dead cells between events if requested)
    // returns energy ( in units determined by the overloaded digitiseDetectorEnergy )

    // this is an overloaded method, provides energy in technology-dependent units
    float e_out = digitiseDetectorEnergy( evtdata._generator, energy ) ;
    // the following make only relative changes to the energy
    std::normal_distribution<float> miscalDistUncorrel( 1.0, _misCalib_uncorrel );
    // random miscalib, uncorrelated in cells
    if ( _misCalib_uncorrel > 0 ) {
      e_out *= miscalDistUncorrel( evtdata._generator ) ;
    }
    // random miscalib, correlated across cells in one event
    if ( _misCalib_correl > 0 ) {
      e_out *= evtdata._eventCorrelMiscalib ;
    }
    // get MIP scale in local unit
    float oneMipInMyUnits = convertEnergy( 1.0, EnergyScale::MIP ) ;
    // limited electronics dynamic range
    if ( _elec_rangeMip > 0 ) {
      e_out = std::min ( e_out, _elec_rangeMip * oneMipInMyUnits ) ;
    }
    // add electronics noise
    if ( _elec_noiseMip > 0 ) {
      std::normal_distribution<float> gauss( 0., _elec_noiseMip * oneMipInMyUnits ) ;
      e_out += gauss( evtdata._generator ) ;
    }
    // random cell kill
    if ( _deadCell_fraction > 0 ) {
      std::uniform_real_distribution<float> flat(0., 1.) ;
      if ( flat( evtdata._generator ) < _deadCell_fraction ) { 
        e_out = 0 ;
      }
    }
    return e_out;
  }

}



