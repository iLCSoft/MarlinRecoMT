// -- marlin reco headers
#include <MarlinRecoMT/OverlayFileHandler.h>
#include <MarlinRecoMT/OverlayMerging.h>

// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/ProcessorApi.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- std headers
#include <random>

namespace marlinreco_mt {

  /** Overlay processor allows to overlay an event with background events from 
   *  additional LCIO files based on different criteria.
   *
   *  A typical use case would be the overlay of gamma gamma -> hadrons background events
   *  with a number drawn from a poissonian distribution with a given mean 'expBG' (NumberOverlayEvents=0).
   *
   *  See Merger.cc for the collection types that can be merged.
   * 
   * @author N. Chiapolini, DESY
   * @author F. Gaede, DESY
   * @author R. Ete, DESY
   * @version $Id$
   * 
   * @param InputFileNames (StringVec) The names (with absolute or relative pathes) of the files from which the background should be read.
   *                                   Multiple files can be given by a white spaces separated list or by setting this parameter multiple times. 
   *                              
   * @param CollectionMap (StringVec)  Pairs of collection names. The input collection (given first) will be merged into the output collection.
   *                                   If the output collection does not exist, it will be created. It is recommended to set this parameter once
   *                                   for each collection pair. If this parameter is not set, all collections with the same name and type will be merged. 
   *
   * @param NumberOverlayEvents (int)  Fixed number of background events that should be added to each physics event. (default 0)
   *
   * @param expBG (double)             If this value is set, a random number of background will be added to each physics event. 
   *                                   The Random numbers will be thrown according to a Poisson distribution with this expectation value.
   *                                   If set, NumberOverlayEvents will be added to the random number.  
   * @param ExcludeCollectionMap (StringVec) List of collection to exclude for merging. This is particularly useful when you just want to exclude a few collections.
   *                                   One doesn't have to specify all collections to overlay in the CollectionMap parameter minus the collection to avoid, 
   *                                   but just the ones to exclude. Priority is given to this list over the CollectionMap.                             
   */
  class OverlayProcessor : public marlin::Processor {
    using RandomGenerator = std::mt19937 ;
    
   public:
    marlin::Processor*  newProcessor() { return new OverlayProcessor ; }

    /** Constructor
     */
    OverlayProcessor() ;

    /** Called at the begin of the job before anything is read.
     *  Use to initialize the processor, e.g. book histograms.
     */
    void init() override ;

    /** Called for every event - the working horse.
     */
    void processEvent( EVENT::LCEvent * evt ) override ;

  private:
    /// Randomly read the next event from the available file
    std::shared_ptr<EVENT::LCEvent> readNextEvent( RandomGenerator &generator ) ;
    /// Get the number of available events in all files
    unsigned int getNAvailableEvents() ;

  protected:
    
    marlin::Property<std::vector<std::string>> _fileNames {this, "InputFileNames" , 
        "Name of the lcio input file(s)", {"undefined.slcio"} } ;
  
    marlin::Property<int> _numOverlay {this, "NumberOverlayEvents" , 
        "Overlay each event with this number of background events. (default 0)" , 0 } ; 
 
    marlin::Property<double> _expBG {this, "expBG" , 
        "Add additional background events according to a poisson distribution with this expectation value. (non, if parameter not set)" , 1. } ;
    
    marlin::Property<std::vector<std::string>> _overlayCollections {this, "CollectionMap" , 
        "Pairs of collection to be merged", {"MCParticle", "MCParticle"} } ;
        
    marlin::Property<std::vector<std::string>> _excludeCollections {this, "ExcludeCollections" , 
        "List of collections to exclude for merging" } ;
    
    // internal members
    /// The total number of available overlay events from input files
    unsigned int                          _nAvailableEvents {0} ;     
    /// The map of collections to overlay, built from _overlayCollections
    std::map<std::string, std::string>    _overlayCollectionMap {} ;  
    /// The total number of processed runs
    int                                   _nRun {0} ;                 
    /// The total number of processed events
    int                                   _nEvt {0} ;                 
    /// The total number of overlaid events when processor ends
    int                                   _nTotalOverlayEvents {0} ;  
    /// The list of file handler to manage overlay input files (see LCFileHandler class)
    OverlayFileHandlerList                _fileHandlerList {} ;     
  };

  //--------------------------------------------------------------------------

  OverlayProcessor::OverlayProcessor() :
    Processor("Overlay") {
    // modify processor description
    _description = "Opens a second (chain of) lcio file(s) and overlays events..." ;
    
    // clone processors as they have to open files on first call
    forceRuntimeOption( Processor::RuntimeOption::Critical, false ) ;
    forceRuntimeOption( Processor::RuntimeOption::Clone, true ) ;
  }

  //--------------------------------------------------------------------------

  void OverlayProcessor::init() {
    // usually a good idea to
    printParameters() ;
    
    // prepare the lcio file handlers
    _fileHandlerList.resize( _fileNames.get().size() ) ;
    
    for ( unsigned int i=0 ; i<_fileNames.get().size() ; i++ ) {
      _fileHandlerList.at( i ).setFileName( _fileNames.get().at( i ) ) ;
    }
  
    // initalisation of random number generator
    marlin::ProcessorApi::registerForRandomSeeds( this ) ;
  
    if ( _overlayCollections.get().size() % 2 ) {
      marlin::ProcessorApi::abort( this, "Odd number of collection names, can't make a correct mapping" ) ;
    }

    // preparing collection map for merge.  
    // treating pairs of collection names
    for ( auto iter = _overlayCollections.get().begin() ; iter != _overlayCollections.get().end() ; ++iter ) {  
      std::string key = *iter ;
      ++iter ;
      _overlayCollectionMap[key] = *iter ;
    }
    
    _nAvailableEvents = getNAvailableEvents() ;
    log<MESSAGE>() << "Overlay::modifyEvent: total number of available events to overlay: " << _nAvailableEvents << std::endl ;
  }

  //--------------------------------------------------------------------------

  void OverlayProcessor::processEvent( EVENT::LCEvent * evt ) {
    // initalisation of random number generator
    auto eventSeed = marlin::ProcessorApi::getRandomSeed( this, evt ) ;
    // local random number generator
    RandomGenerator generator {} ;
    generator.seed( eventSeed ) ;
    std::poisson_distribution<int> poissonDistribution { _expBG } ;
    // number of bkg events to overlay
    unsigned int nEventsToOverlay = _numOverlay ;
    if ( parameterSet("expBG") ) {
      nEventsToOverlay += poissonDistribution( generator ) ;
    }
    log<DEBUG6>() << "** Processing event nr " << evt->getEventNumber() << " run " <<  evt->getRunNumber() 
     			    << "\n**  overlaying " << nEventsToOverlay << " background events. \n " 
			    << " ( seeded CLHEP::HepRandom with seed = " << eventSeed  << ") " 
			    << std::endl ;
  
    int nOverlaidEvents(0) ;
    EVENT::FloatVec overlaidEventIDs, overlaidRunIDs ;
    
    for(unsigned int i=0 ; i < nEventsToOverlay ; i++ ) {

      auto overlayEvent = readNextEvent( generator ) ;

      if( nullptr == overlayEvent ) {
	       log<ERROR>() << "loop: " << i << " ++++++++++ Nothing to overlay +++++++++++ \n " ;
	       continue ;
      } 
      
      overlaidEventIDs.push_back( overlayEvent->getEventNumber() );
      overlaidRunIDs.push_back( overlayEvent->getRunNumber() );
      
      ++nOverlaidEvents ;

      log<DEBUG6>() << "loop: " << i << " will overlay event " << overlayEvent->getEventNumber() << " - run " << overlayEvent->getRunNumber() << std::endl ;

      std::map<std::string, std::string> collectionMap;
      
      if ( ( _overlayCollectionMap.empty() || ! parameterSet("CollectionMap") ) ) {
        auto collectionNames = overlayEvent->getCollectionNames();
        for ( auto collection : *collectionNames ) {
          collectionMap[collection] = collection;
          log<DEBUG6>() << "Collection map -> " << collection << std::endl;
        }
      }
      else {
        collectionMap = _overlayCollectionMap ;
      }
      
      // Remove collections to exclude from the collection map
      if(not _excludeCollections.get().empty()) {
        for ( auto excludeCol : _excludeCollections.get() ) {
          auto findIter = collectionMap.find( excludeCol );
          if( collectionMap.end() != findIter ) {
            collectionMap.erase( findIter );
          }
        }
      }
	     OverlayMerging::mergeEvents( overlayEvent.get(), evt, collectionMap );
    }
    
    _nTotalOverlayEvents += nOverlaidEvents ;
    
    // Write info to event parameters
    std::string paramName = "Overlay." + this->name() + ".nEvents" ;
    evt->parameters().setValue(paramName, nOverlaidEvents) ;
    paramName = "Overlay." + this->name() + ".eventIDs" ;
    evt->parameters().setValues(paramName, overlaidEventIDs) ;
    paramName = "Overlay." + this->name() + ".runIDs" ;
    evt->parameters().setValues(paramName, overlaidRunIDs) ;
    int totalOverlay = 0 ;
    try {
      // for now: returns 0 if the key doesn't exists.
      // in future: might throw an exception if not found...
      totalOverlay = evt->parameters().getIntVal("Overlay.nTotalEvents") ; 
    }
    catch(...) {}
    totalOverlay += nOverlaidEvents ;
    evt->parameters().setValue("Overlay.nTotalEvents", totalOverlay) ;
  }

  //--------------------------------------------------------------------------

  std::shared_ptr<EVENT::LCEvent> OverlayProcessor::readNextEvent( RandomGenerator &generator ) {
    // get the event index to random pick an event among the possible files
    std::uniform_int_distribution<int> flatDistribution( 0, _nAvailableEvents ) ;
    const unsigned int eventIndex = flatDistribution( generator ) ;
    unsigned int currentEventIndex(0);    
    log<DEBUG>() << "Overlay::readNextEvent: index = " << eventIndex  << " over " << _nAvailableEvents << std::endl ;
    for ( auto &handler : _fileHandlerList ) {
      if ( currentEventIndex <= eventIndex && eventIndex < currentEventIndex + handler.getNumberOfEvents() ) {        
        const auto eventNumber = handler.getEventNumber( eventIndex-currentEventIndex ) ;
        const auto runNumber = handler.getRunNumber( eventIndex-currentEventIndex ) ;
        return handler.readEvent( runNumber, eventNumber ) ;
      }
      currentEventIndex += handler.getNumberOfEvents() ;
    }    
    return nullptr ;
  }
  
  //--------------------------------------------------------------------------
  
  unsigned int OverlayProcessor::getNAvailableEvents() {
    unsigned int totalNEvents {0} ;
    for ( auto &handler : _fileHandlerList ) {
      totalNEvents += handler.getNumberOfEvents() ;
    }
    return totalNEvents;
  }

  // processor declaration
  OverlayProcessor anOverlayProcessor ;
}
