// -- lcio headers
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/BitField64.h>

// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/ProcessorApi.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- std headers
#include <functional>
#include <string>
#include <memory>

// -- marlin reco mt headers
#include <MarlinRecoMT/LCIOHelper.h>

namespace marlinreco_mt {

  /** Utility processor that allows to split a collection of Hits into 
   *  several collections based on the layer information in the cellID word.
   *  Works for all four lcio hit classes.
   *
   *  @parameter InputCollection name of the hit collection with (Sim)TrackerHits/(Sim)CalorimeterHits
   *  @parameter OutputCollections ( ColName  StartLayer EndLayer )    
   * 
   * @author F. Gaede, CERN/DESY
   * @author R.Ete , DESY
   * @date  30 Oct 2014
   * @date  06 Aug 2019
   * @version $Id: $
   */
  class SplitCollectionByLayerProcessor : public marlin::Processor {
    /// helper struct
    struct OutputCollectionInfo {
      std::string                                 _name {} ;
      unsigned                                    _layerStart {0} ;
      unsigned                                    _layerEnd {0} ;
      std::unique_ptr<IMPL::LCCollectionVec>      _collection {nullptr} ;
    };

   public:
    marlin::Processor*  newProcessor() { return new SplitCollectionByLayerProcessor ; }

    /** Constructor
     */
    SplitCollectionByLayerProcessor() ;

    /** Called at the begin of the job before anything is read.
     *  Use to initialize the processor, e.g. book histograms.
     */
    void init() ;

    /** Called for every event - the working horse.
     */
    void processEvent( EVENT::LCEvent * evt ) ;

  protected:
    marlin::Property<std::string> _inputCollectionName {this, "InputCollection" , 
             "Name of the input collection with hits", "FTDCollection" } ;

    marlin::Property<std::vector<std::string>> _collectionsAndLayers {this, "OutputCollections" , 
             "Name of the output collection with start and end layer number" , { "FTD_PIXELCollection", "0", "1", "FTD_STRIPCollection", "2", "6" } } ;
    
    std::vector<OutputCollectionInfo>         _outputCollections {} ;
  };

  //--------------------------------------------------------------------------

  SplitCollectionByLayerProcessor::SplitCollectionByLayerProcessor() :
    Processor("SplitCollectionByLayer") {
    // modify processor description
    _description = "split a hit collection based on the layer number of the hits " ;
  }

  //--------------------------------------------------------------------------

  void SplitCollectionByLayerProcessor::init() {
    // usually a good idea to
    printParameters() ;
    
    // chekck consistency
    if( 0 != _collectionsAndLayers.get().size() % 3 ) {
      marlin::ProcessorApi::abort( this, "The OutputCollections parameter length should be a multiple of 3 (CollectionName layer0 layer1)." ) ;
    }
    std::size_t len = _collectionsAndLayers.get().size() % 3 ;
    _outputCollections.resize( len ) ;
    
    std::size_t index = 0 ;
    for( std::size_t i=0 ; i<len ; i++ ) {
      _outputCollections[i]._name        = _collectionsAndLayers.get()[ index ] ; index ++ ;
      _outputCollections[i]._layerStart  = std::atoi( _collectionsAndLayers.get()[ index ].c_str() ) ; index ++ ;
      _outputCollections[i]._layerEnd    = std::atoi( _collectionsAndLayers.get()[ index ].c_str() ) ; index ++ ;
    }
  }

  //--------------------------------------------------------------------------

  void SplitCollectionByLayerProcessor::processEvent( EVENT::LCEvent * evt ) {
    EVENT::LCCollection *collection = nullptr ;
    try{   
      collection =  evt->getCollection( _inputCollectionName )  ; 
    } 
    catch( EVENT::DataNotAvailableException& ) {
      log<DEBUG5>() <<  " input collection not in event : " << _inputCollectionName.get() << "   - nothing to do    !!! " << std::endl ;  
      return ;
    }
    // cellID converter function
    std::function<long long( const EVENT::LCObject* )> cellIDFunction = nullptr ;
    
    // remember the type of the hit collection
    if( collection->getTypeName() == lcio::LCIO::SIMTRACKERHIT ) {
      cellIDFunction = std::bind(&LCIOHelper::cellIDToLong<EVENT::SimTrackerHit>, std::placeholders::_1 ) ;
    }
    else if( collection->getTypeName() == lcio::LCIO::TRACKERHIT ) {
      cellIDFunction = std::bind(&LCIOHelper::cellIDToLong<EVENT::TrackerHit>, std::placeholders::_1 ) ;
    }
    else if( collection->getTypeName() == lcio::LCIO::SIMCALORIMETERHIT ) {
      cellIDFunction = std::bind(&LCIOHelper::cellIDToLong<EVENT::SimCalorimeterHit>, std::placeholders::_1 ) ;
    }
    else if( collection->getTypeName() == lcio::LCIO::CALORIMETERHIT ) {
      cellIDFunction = std::bind(&LCIOHelper::cellIDToLong<EVENT::CalorimeterHit>, std::placeholders::_1 ) ;
    }
    else {
      log<WARNING>() <<  " input collection unexpected type : " << collection->getTypeName() << "   - skipping    !!! " << std::endl ;  
      return ;
    }
    std::string encoderString = collection->getParameters().getStringVal( "CellIDEncoding" ) ;
    UTIL::BitField64 encoder( encoderString ) ;
    unsigned layerIndex = encoder.index("layer") ;
    
    // create output collections
    for( auto &outcol : _outputCollections ) {
      outcol._collection = std::make_unique<IMPL::LCCollectionVec>( collection->getTypeName() ) ;
      outcol._collection->setSubset( true ) ;
      outcol._collection->parameters().setValue( "CellIDEncoding", encoderString ) ;
      log<DEBUG5>() << " create new output collection " << outcol._name << " of type " <<  collection->getTypeName() << std::endl ;
    }
    // loop over hits
    int nHit = collection->getNumberOfElements()  ;
    for( int iHit=0; iHit< nHit ; iHit++ ) {
      auto h = collection->getElementAt( iHit ) ;
      auto id = cellIDFunction( h ) ;
      encoder.setValue( id ) ;
      unsigned int layerID = encoder[ layerIndex ] ;
      // check if we have an output collection for this layer
      for( auto &outcol : _outputCollections ) {
        if( ( outcol._layerStart <= layerID )  && ( layerID <= outcol._layerEnd ) ) {
          outcol._collection->addElement( h ) ;
          log<DEBUG0>() << " adding hit for layer " << layerID << " to collection : " << outcol._name << std::endl ;
        }    
      }
    }
    // add non empty collections to the event
    for( auto &outcol : _outputCollections ) {
      if( outcol._collection->getNumberOfElements() > 0 ) {
        evt->addCollection( outcol._collection.release(), outcol._name ) ;
        log<DEBUG5>() << " output collection " << outcol._name << " of type " <<  collection->getTypeName() << " added to the event  " << std::endl ;
      }
    }
  }

  // processor declaration
  SplitCollectionByLayerProcessor aSplitCollectionByLayerProcessor ;
}
