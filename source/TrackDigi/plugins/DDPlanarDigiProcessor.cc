// -- lcio headers
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCIO.h>
#include <UTIL/CellIDEncoder.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/ProcessorApi.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- dd4hep headers
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

// -- std headers
#include <random>

namespace marlinreco_mt {

  /** ======= DDPlanarDigiProcessor ========== <br>
   * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
   * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
   * of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions. 
   * The geometry of the surface is retreived from DDRec::Surface associated to the hit via cellID.
   * 
   * 
   * <h4>Input collections and prerequisites</h4> 
   * Processor requires a collection of SimTrackerHits <br>
   * <h4>Output</h4>
   * Processor produces collection of smeared TrackerHits<br>
   * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
   * (default name VXDCollection) <br>
   * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
   * (default name VTXTrackerHits) <br>
   * @param ResolutionU resolution in direction of u (in mm) <br>
   * (default value 0.004) <br>
   * @param ResolutionV Resolution in direction of v (in mm) <br>
   * (default value 0.004) <br>
   * @param IsStrip whether the hits are 1 dimensional strip measurements <br>
   * (default value false)
   * @param Ladder_Number_encoded_in_cellID ladder number has been encoded in the cellID <br>
   * (default value false) <br>
   * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio <br>
   * (default value lcio::ILDDetID::VXD) <br>
   * <br>
   * 
   * @author F.Gaede CERN/DESY, S. Aplin DESY
   * @date Dec 2014
   */
  class DDPlanarDigiProcessor : public marlin::Processor {
    using RandomGenerator = std::mt19937 ;
    
    static constexpr unsigned int SmearingNMaxTries = 10 ; 
    
  public:
    ~DDPlanarDigiProcessor() = default ;
    DDPlanarDigiProcessor(const DDPlanarDigiProcessor&) = delete ;
    DDPlanarDigiProcessor& operator=(const DDPlanarDigiProcessor&) = delete ;
    
    marlin::Processor*  newProcessor() { return new DDPlanarDigiProcessor ; }

    /**
     *  @brief  Constructor
     */
    DDPlanarDigiProcessor() ;

    /** Called at the begin of the job before anything is read.
     * Use to initialize the processor, e.g. book histograms.
     */
    void init() ;

    /** Called for every event - the working horse.
     */
    void processEvent( EVENT::LCEvent * evt ) ;

  protected:
    // processor parameters
    std::string                _inputCollectionName {} ;
    std::string                _outputCollectionName {} ;
    std::string                _outputRelCollectionName {} ;
    std::string                _subDetectorName {} ;
    EVENT::FloatVec            _resolutionU {} ;
    EVENT::FloatVec            _resolutionV {} ;
    bool                       _isStrip {false} ;
    bool                       _forceHitsOntoSurface {false} ;
    double                     _minEnergy {0.} ;
    // to be replaced by std random stuff
    // gsl_rng* _rng ;
    const dd4hep::rec::SurfaceMap* _surfaceMap {nullptr} ;
  };

  //--------------------------------------------------------------------------

  DDPlanarDigiProcessor::DDPlanarDigiProcessor() :
    Processor("DDPlanarDigiProcessor") {
    // modify processor description
    _description = "DDPlanarDigiProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters."
      "The geoemtry of the surface is taken from the DDRec::Surface asscociated to the hit via the cellID" ;

    registerProcessorParameter( "ResolutionU" ,
                                "resolution in direction of u - either one per layer or one for all layers "  ,
                                _resolutionU ,
                                EVENT::FloatVec {0.004} ) ;

    registerProcessorParameter( "ResolutionV" , 
                                "resolution in direction of v - either one per layer or one for all layers " ,
                               _resolutionV ,
                                EVENT::FloatVec {0.004} );

    registerProcessorParameter( "IsStrip",
                                "whether hits are 1D strip hits",
                                _isStrip,
                                _isStrip );
    
    
    registerProcessorParameter( "SubDetectorName" , 
                               "Name of dub detector" ,
                               _subDetectorName ,
                                std::string("VXD") );
      
    // Input collections
    registerInputCollection( EVENT::LCIO::SIMTRACKERHIT,
                            "SimTrackHitCollectionName" , 
                            "Name of the Input SimTrackerHit collection"  ,
                            _inputCollectionName ,
                            std::string("VXDCollection") ) ;
    
    // Output collections
    registerOutputCollection( EVENT::LCIO::TRACKERHITPLANE,
                             "TrackerHitCollectionName" , 
                             "Name of the TrackerHit output collection"  ,
                             _outputCollectionName ,
                             std::string("VTXTrackerHits") ) ;
    
    registerOutputCollection( EVENT::LCIO::LCRELATION,
                             "SimTrkHitRelCollection",
                             "Name of TrackerHit SimTrackHit relation collection",
                             _outputRelCollectionName,
                             std::string("VTXTrackerHitRelations"));
    
    registerProcessorParameter( "ForceHitsOntoSurface" , 
                                "Project hits onto the surface in case they are not yet on the surface (default: false)" ,
                                _forceHitsOntoSurface ,
                                _forceHitsOntoSurface );

    registerProcessorParameter( "MinimumEnergyPerHit" ,
                                "Minimum Energy (in GeV!) to accept hits, other hits are ignored",
                                _minEnergy,
                                _minEnergy );
  }

  //--------------------------------------------------------------------------

  void DDPlanarDigiProcessor::init() {
    // usually a good idea to
    printParameters() ;
    
    // initalisation of random number generator
    marlin::ProcessorApi::registerForRandomSeeds( this ) ;
    
    if( _resolutionU.size() != _resolutionV.size() ) {
      std::stringstream ss ;
      ss << name() << "::init() - Inconsistent number of resolutions given for U and V coordinate: " 
         << "ResolutionU  :" <<   _resolutionU.size() << " != ResolutionV : " <<  _resolutionV.size() ;
      marlin::ProcessorApi::abort( this, ss.str() ) ;
    }
    
    //===========  get the surface map from the SurfaceManager ================
    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
    dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;
    dd4hep::DetElement det = theDetector.detector( _subDetectorName ) ;
    _surfaceMap = surfMan.map( det.name() ) ;
    if( nullptr == _surfaceMap ) {   
      std::stringstream err  ; 
      err << " Could not find surface map for detector: " << _subDetectorName << " in SurfaceManager " ;
      marlin::ProcessorApi::abort( this, err.str() ) ;
    }
    
    log<DEBUG3>() << " DDPlanarDigiProcessor::init(): found " << _surfaceMap->size() 
                            << " surfaces for detector:" <<  _subDetectorName << std::endl ;
  }

  //--------------------------------------------------------------------------

  void DDPlanarDigiProcessor::processEvent( EVENT::LCEvent * evt ) {
    // initalisation of random number generator
    auto eventSeed = marlin::ProcessorApi::getRandomSeed( this, evt ) ;
    log<DEBUG4>() << "seed set to " << eventSeed << std::endl ;
    RandomGenerator generator {} ;
    generator.seed( eventSeed ) ;
    std::normal_distribution<double> gaussian {} ;
    // get the input collection
    EVENT::LCCollection *inputCollection = nullptr ;
    try {
      inputCollection = evt->getCollection( _inputCollectionName ) ;
    }
    catch( EVENT::DataNotAvailableException &) {
      log<DEBUG4>() << "Collection " << _inputCollectionName.c_str() << " is unavailable in event " << evt->getEventNumber() << std::endl ;
      return ;
    }
    // output collections
    auto outputCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::TRACKERHITPLANE ) ;
    auto outputRelCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::LCRELATION ) ;
    // to store the weights
    IMPL::LCFlagImpl lcFlag( 0 ) ;
    lcFlag.setBit( EVENT::LCIO::LCREL_WEIGHTED ) ;
    outputRelCollection->setFlag( lcFlag.getFlag() ) ;
    // cellID utils
    UTIL::CellIDEncoder<IMPL::TrackerHitPlaneImpl> cellid_encoder( UTIL::LCTrackerCellID::encoding_string() , outputCollection.get() ) ;
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> cellid_decoder( inputCollection ) ;
    
    int nSimHits = inputCollection->getNumberOfElements() ;
    log<DEBUG4>() << " processing collection " << _inputCollectionName  << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    unsigned nCreatedHits = 0 ;
    unsigned nDismissedHits = 0 ;
    
    for( int i=0 ; i<nSimHits ; ++i ) {
      auto simTHit = dynamic_cast<EVENT::SimTrackerHit*>( inputCollection->getElementAt( i ) ) ;

      if( simTHit->getEDep() < _minEnergy ) {
        log<DEBUG>() << "Hit with insufficient energy " << simTHit->getEDep()*1e6 << " keV" << std::endl ;
        continue;
      }
      const int cellID0 = simTHit->getCellID0() ;

      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************

      auto findIter = _surfaceMap->find( cellID0 ) ;

      if( findIter == _surfaceMap->end() ) {
        std::stringstream err ; 
        err << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " << cellid_decoder( simTHit ).valueString() ;
        marlin::ProcessorApi::abort( this, err.str() ) ;
      }
      
      const auto surf = findIter->second ;
      int layer = cellid_decoder( simTHit )["layer"] ;
      dd4hep::rec::Vector3D oldPos( simTHit->getPosition()[0], simTHit->getPosition()[1], simTHit->getPosition()[2] ) ;
      dd4hep::rec::Vector3D newPos ;
      
      //************************************************************
      // Check if Hit is inside senstive 
      //************************************************************ 
                                   
      if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {      
        if( _forceHitsOntoSurface ) {
          dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
          dd4hep::rec::Vector3D oldPosOnSurf = (1./dd4hep::mm) * surf->localToGlobal( lv ) ; 
          log<DEBUG3>() << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf-oldPos).r() << std::endl ;       
          oldPos = oldPosOnSurf ;
        } 
        else {
          ++nDismissedHits;       
          continue; 
        }
      }
      //**************************************************************************
      // Try to smear the hit but ensure the hit is inside the sensitive region
      //**************************************************************************
      dd4hep::rec::Vector3D u = surf->u() ;
      dd4hep::rec::Vector3D v = surf->v() ;
      // get local coordinates on surface
      dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
      double uL = lv[0] / dd4hep::mm ;
      double vL = lv[1] / dd4hep::mm ;
      bool accept_hit = false ;
      unsigned  tries   =  0 ;              
      float resU = ( _resolutionU.size() > 1 ?   _resolutionU.at(  layer )     : _resolutionU.at(0)   )  ;
      float resV = ( _resolutionV.size() > 1 ?   _resolutionV.at(  layer )     : _resolutionV.at(0)   )  ; 
      
      while( tries <  DDPlanarDigiProcessor::SmearingNMaxTries ) {
      
        if( tries > 0 ) {
          log<DEBUG0>() << "retry smearing for " <<  cellid_decoder( simTHit ).valueString() << " : retries " << tries << std::endl ;
        } 
        double uSmear = gaussian( generator, std::normal_distribution<double>::param_type( 0., resU ) ) ;
        double vSmear = gaussian( generator, std::normal_distribution<double>::param_type( 0., resV ) ) ;
        dd4hep::rec::Vector3D newPosTmp = 1./dd4hep::mm  * 
        ( ! _isStrip  ? surf->localToGlobal( dd4hep::rec::Vector2D (  ( uL + uSmear ) * dd4hep::mm, ( vL + vSmear )  *dd4hep::mm ) )  :
                        surf->localToGlobal( dd4hep::rec::Vector2D (  ( uL + uSmear ) * dd4hep::mm,          0.                  ) ) ) ;
        log<DEBUG1>() << " hit at    : " << oldPos 
                                << " smeared to: " << newPosTmp
                                << " uL: " << uL 
                                << " vL: " << vL 
                                << " uSmear: " << uSmear
                                << " vSmear: " << vSmear
                                << std::endl ;
        if ( surf->insideBounds( dd4hep::mm * newPosTmp ) ) { 
          accept_hit = true ;
          newPos     = newPosTmp ;
          break;  
        } 
        else {   
          log<DEBUG1>() << "  hit at " << newPosTmp 
                                  << " " << cellid_decoder( simTHit).valueString() 
                                  << " is not on surface " 
                                  << " distance: " << surf->distance( dd4hep::mm * newPosTmp ) 
                                  << std::endl;        
        }
        ++tries;
      }
      if( not accept_hit ) {
        log<DEBUG4>() << "hit could not be smeared within ladder after " <<  DDPlanarDigiProcessor::SmearingNMaxTries << "  tries: hit dropped"  << std::endl ;
        ++nDismissedHits ;
        continue ; 
      }
      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************
      const int cellID1 = simTHit->getCellID1() ;
      float u_direction[2] ;
      u_direction[0] = u.theta();
      u_direction[1] = u.phi();    
      float v_direction[2] ;
      v_direction[0] = v.theta();
      v_direction[1] = v.phi();
      auto trkHit = std::make_unique<IMPL::TrackerHitPlaneImpl>() ;
      trkHit->setCellID0( cellID0 ) ;
      trkHit->setCellID1( cellID1 ) ;
      trkHit->setPosition( newPos.const_array()  ) ;
      trkHit->setTime( simTHit->getTime() ) ;
      trkHit->setEDep( simTHit->getEDep() ) ;
      trkHit->setU( u_direction ) ;
      trkHit->setV( v_direction ) ;
      trkHit->setdU( resU ) ;    
      log<DEBUG0>() << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1] 
                    << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
                    << std::endl ;
      if( _isStrip ) {
        // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
        double stripRes = (surf->length_along_v() / dd4hep::mm ) / std::sqrt( 12. ) ;
        trkHit->setdV( stripRes ); 
      } 
      else {
        trkHit->setdV( resV ) ;
      }
      if( _isStrip ) {
        trkHit->setType( UTIL::set_bit( trkHit->getType(), UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) ) ;
      }
      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************           
      auto rel = new IMPL::LCRelationImpl() ;
      rel->setFrom ( trkHit.get() ) ;
      rel->setTo ( simTHit );
      rel->setWeight( 1.0 ) ;
      outputRelCollection->addElement( rel ) ;
      //**************************************************************************
      // Add hit to collection
      //**************************************************************************    
      outputCollection->addElement( trkHit.release() ) ; 
      ++nCreatedHits ;
      log<DEBUG3>() << "-------------------------------------------------------" << std::endl ;
    }
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
    evt->addCollection( outputCollection.release()    , _outputCollectionName    ) ;
    evt->addCollection( outputRelCollection.release() , _outputRelCollectionName ) ;
    log<DEBUG4>() << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed as not on sensitive element" << std::endl ;
  }

  // processor declaration
  DDPlanarDigiProcessor aDDPlanarDigiProcessor ;
}
