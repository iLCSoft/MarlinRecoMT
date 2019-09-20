// -- lcio headers
#include <IMPL/LCCollectionVec.h>
// #include <LCRTRelations.h>

// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/Logging.h>
#include <marlin/ProcessorApi.h>

namespace marlinreco_mt {

  /** Helper processor that merges several input collections into a transient subset collections.
   *  The names and optionally the IDs of the merged collections
   *  are stored in collection parameters MergedCollectionNames and MergedCollectionIDs.
   *
   * @param InputCollections    Name of the input collections
   * @param InputCollectionIDs  Optional IDs for input collections - if given, IDs will be added to all objects in merged collections as ext<CollID>()"
   *                            - it is the users responsibility to ensure uniqueness of the IDs across the event ( and that ID != 0 )
   * @param OutputCollection    Name of the output collection
   *
   * @author R. Ete, DESY
   * @author F. Gaede, DESY
   * @author B. Vormwald, DESY
   */
  class MergeCollections : public marlin::Processor {
   public:
    marlin::Processor*  newProcessor() { return new MergeCollections ; }

    /**
     *  @brief  Constructor
     */
    MergeCollections() ;

    /** Called at the begin of the job before anything is read.
     * Use to initialize the processor, e.g. book histograms.
     */
    void init() ;

    /** Called for every event - the working horse.
     */
    void processEvent( EVENT::LCEvent * evt ) ;

  private:
    ///< Helper function to get collection safely
    EVENT::LCCollection* getCollection( EVENT::LCEvent* evt, const std::string name ) const ;

  protected:
    marlin::Property<std::vector<std::string>> _inColNames {this, "InputCollections" ,
              "Names of all input collections" } ;

    marlin::Property<std::vector<int>> _inColIDs {this, "InputCollectionIDs" ,
              "IDs for input collections - if given id will be added to all objects in merged collections as ext<CollID)" } ;

    marlin::OutputCollectionProperty _outColName  {this, "OutputCollection" ,
              "Name of output collection", "MergedCollection" } ;

    marlin::Property<int> _collectionParameterIndex {this, "CollectionParameterIndex",
              "Index of input collection  that is used to copy the  collection parameters from " , 0 } ;
  };

  //--------------------------------------------------------------------------

  MergeCollections::MergeCollections() :
    Processor("MergeCollections") {
    // modify processor description
    _description = "MergeCollections creates a transient subset collection that merges all input collections " ;
  }

  //--------------------------------------------------------------------------

  void MergeCollections::init() {
    // usually a good idea to
    printParameters() ;
  }

  //--------------------------------------------------------------------------

  void MergeCollections::processEvent( EVENT::LCEvent * evt ) {
    std::vector<std::string> colNamesPresent;
    std::vector<int> colIDsPresent;
    std::vector<int> colNElements;
    std::vector<int> colNIntParam;
    std::vector<int> colNFloatParam;
    std::vector<int> colNStringParam;

    std::size_t nCol = _inColNames.get().size() ;
    std::size_t nColID = _inColIDs.get().size() ;

    if( marlin::ProcessorApi::isFirstEvent( evt ) && nColID != nCol ) {
      log<marlin::WARNING>() << " MergeCollections::processEvent : incompatible parameter vector sizes : InputCollections: " << nCol
                             << " <->  InputCollectionIDs " << nColID << std::endl;
      log<marlin::WARNING>() << " MergeCollections::processEvent : standard numbering (0,1,2,...) used." << std::endl;
    }

    //--- copy existing collections to a vector first
    std::vector<EVENT::LCCollection*> colVec( nCol ) ;

    for( std::size_t i=0 ; i<nCol ; ++i ) {
      EVENT::LCCollection *col = getCollection( evt , _inColNames.get()[i] ) ;
      if( col != 0 ) {
        colVec[i] = col  ;
        colNamesPresent.push_back(_inColNames.get()[i]);
        if( nColID == nCol ) {
          colIDsPresent.push_back(_inColIDs.get()[i]);
        }
        else {
          colIDsPresent.push_back(i) ;
        }
      }
      else {
        log<marlin::DEBUG2>() << " input collection missing : " << _inColNames.get()[i] << std::endl ;
      }
    }

    //--- now loop over collections
    IMPL::LCCollectionVec* outCol = 0 ;
    bool first = true ;

    for( std::size_t k=0 ; k<nCol ; ++k ) {
      EVENT::LCCollection* col = colVec[k] ;
      if( ! col ) {
        continue ;
      }
      if( first ){
        // copy collection flags from first collections
        outCol = new IMPL::LCCollectionVec( col->getTypeName() )  ;
        outCol->setFlag( col->getFlag() ) ;
        first = false ;
      }
      int nEle = col->getNumberOfElements() ;
      for( int j=0 ; j<nEle ; ++j ) {
        EVENT::LCObject* elem = col->getElementAt(j) ;
        outCol->addElement(  elem ) ;
      }

      int intParams = 0;
      int floatParams = 0;
      int stringParams = 0;

      std::vector<std::string> intKeys ;
      int nIntParameters = col->getParameters().getIntKeys( intKeys ).size() ;
      for(int i=0; i< nIntParameters ; i++ ){
        std::vector<int> intVec ;
        col->getParameters().getIntVals(  intKeys[i], intVec ) ;
        const std::string newIntKey = _inColNames.get()[k]+"_"+intKeys[i];
        outCol->parameters().setValues(newIntKey,intVec);
        intParams++;
        if( unsigned(_collectionParameterIndex) == k )
  	      outCol->parameters().setValues(intKeys[i],intVec);
      }

      std::vector<std::string> floatKeys ;
      int nFloatParameters = col->getParameters().getFloatKeys( floatKeys ).size() ;
      for(int i=0; i< nFloatParameters ; i++ ){
        std::vector<float> floatVec ;
        col->getParameters().getFloatVals(  floatKeys[i], floatVec ) ;
        const std::string newFloatKey = _inColNames.get()[k]+"_"+floatKeys[i];
        outCol->parameters().setValues(newFloatKey,floatVec);
        floatParams++;
        if( unsigned(_collectionParameterIndex) == k )
  	      outCol->parameters().setValues(floatKeys[i],floatVec);
      }

      std::vector<std::string> stringKeys ;
      int nStringParameters = col->getParameters().getStringKeys( stringKeys ).size() ;
      for(int i=0; i< nStringParameters ; i++ ){
        std::vector<std::string> stringVec ;
        col->getParameters().getStringVals(  stringKeys[i], stringVec ) ;
        const std::string newStringKey = _inColNames.get()[k]+"_"+stringKeys[i];
        outCol->parameters().setValues(newStringKey,stringVec);
        stringParams++;
        if( unsigned(_collectionParameterIndex) == k )
  	      outCol->parameters().setValues(stringKeys[i],stringVec);
      }
      colNElements.push_back(nEle);
      colNIntParam.push_back(intParams);
      colNFloatParam.push_back(floatParams);
      colNStringParam.push_back(stringParams);
    }
    if( outCol ) {
      outCol->parameters().setValues("MergedCollection_Names",_inColNames);
      outCol->parameters().setValues("MergedCollection_IDs",_inColIDs);
      outCol->parameters().setValues("MergedCollection_NamesPresent",colNamesPresent);
      outCol->parameters().setValues("MergedCollection_IDsPresent",colIDsPresent);
      outCol->parameters().setValues("MergedCollection_NElements",colNElements);
      outCol->parameters().setValues("MergedCollection_NIntParameters",colNIntParam);
      outCol->parameters().setValues("MergedCollection_NFloatParameters",colNFloatParam);
      outCol->parameters().setValues("MergedCollection_NStringParameters",colNStringParam);
      outCol->setTransient( false ) ;
      outCol->setSubset( true ) ;
      evt->addCollection( outCol, _outColName   ) ;
    }
  }

  //--------------------------------------------------------------------------

  EVENT::LCCollection* MergeCollections::getCollection( EVENT::LCEvent* evt, const std::string name ) const {
    if( name.size() == 0 )
      return 0 ;
    try {
      return evt->getCollection( name ) ;
    }
    catch( EVENT::DataNotAvailableException& e ) {
      log<marlin::DEBUG2>() << "getCollection :  DataNotAvailableException : " << name <<  std::endl ;
      return 0 ;
    }
  }

  // processor declaration
  MergeCollections aMergeCollections ;
}
