// -- lcio headers
#include <EVENT/LCIO.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/ILDConf.h>

// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/Logging.h>
#include <marlin/ProcessorApi.h>
using namespace marlin::loglevel ;

// -- dd4hep headers
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Vector2D.h"
#include "DDRec/Vector3D.h"
#include "DDRec/ISurface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/SurfaceHelper.h"
#include "DDRec/DetectorData.h"

// -- std headers
#include <memory>

// -- root headers
#include <Math/Cartesian3D.h>
#include <Math/RotationZYX.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TRotation.h>

namespace marlinreco_mt {

  /** ================= FTD Space Point Builder =================
   * 
   * Builds space points for pairs of silicon strip detectors. 
   * The digitisers create TrackerHitPlanars for the front and the back strips. 
   * In order to get a spacepoint (as is needed by track reconstruction) those strip measurements need to be combined into one space point.
   * 
   * This is done by this processor. 
   * 
   *  <h4>Input - Prerequisites</h4>
   *  
   * The TrackerHitPlanars as created by the Digitisers of FTD, SET or SIT. 
   * This could of course be used for other detectors as well, but as information about the detector
   * is acquired from DDRec, and as different detectors (at the moment) are stored differently
   * in DDRec, used detectors have to be taken care of in the code.
   *
   *  <h4>Output</h4> 
   *  
   * A collection of TrackerHits containing the constructed space points.
   * The created hits will store the original strip hits in their rawHits. 
   * 
   * @param TrackerHitCollection The name of the input collection of TrackerHits coming from strip detectors on the FTD, SIT or SET <br>
   * (default name FTDTrackerHits) <br>
   * 
   * @param SpacePointsCollection The name of the output collection of the created spacepoints <br>
   * (default name FTDSpacePoints) <br>
   * 
   * @param TrackerHitSimHitRelCollection The name of the input collection of the relations of the TrackerHits to SimHits<br>
   * (default name FTDTrackerHitRelations)<br>
   * 
   * @param SimHitSpacePointRelCollection The name of the SpacePoint SimTrackerHit relation output collection <br>
   * (default name VTXTrackerHitRelations) <br>
   * 
   * @author Robin Glattauer HEPHY, Vienna
   *
   */
  class DDSpacePointBuilderProcessor : public marlin::Processor {
  public:
    static constexpr float crossingPointEpsilon = 0.00001 ;
    using RotationXYZ = ROOT::Math::RotationZYX ;
    using PositionXYZ = ROOT::Math::XYZPoint ;
    using VectorXYZ   = ROOT::Math::XYZVectorF ;
    
  public:
    struct EventStatistics {
      unsigned int _createdSpacePoints {0} ;
      unsigned int _rawStripHits {0} ;
      unsigned int _possibleSpacePoints {0} ;
      unsigned int _nOutOfBoundary {0} ;
      unsigned int _nStripsTooParallel {0} ;
      unsigned int _nPlanesNotParallel {0} ;
    };
  public:
    DDSpacePointBuilderProcessor( const DDSpacePointBuilderProcessor & ) = delete ;
    DDSpacePointBuilderProcessor& operator=( const DDSpacePointBuilderProcessor & ) = delete ;
     
    marlin::Processor*  newProcessor() { return new DDSpacePointBuilderProcessor ; }

    /**
     *  @brief  Constructor
     */
    DDSpacePointBuilderProcessor() ;

    /** Called at the begin of the job before anything is read.
     * Use to initialize the processor, e.g. book histograms.
     */
    void init() ;

    /** Called for every event - the working horse.
     */
    void processEvent( EVENT::LCEvent * evt ) ;

  private:
    EVENT::LCCollection *getCollection( EVENT::LCEvent* evt, const std::string name ) const ;
    std::shared_ptr<UTIL::LCRelationNavigator> createNavigator( EVENT::LCEvent* evt, const std::string name ) const ;
    std::vector<int> getCellID0sAtBack( int cellID0 ) const ;
    std::string getCellID0Info( int cellID0 ) const ;
    int calculateCrossingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y ) const ;
    int calculatePointBetweenTwoLines( const PositionXYZ& P1, const PositionXYZ& V1, const PositionXYZ& P2, const PositionXYZ& V2, PositionXYZ& point ) const ;
    IMPL::TrackerHitImpl* createSpacePoint( EVENT::TrackerHitPlane* a , EVENT::TrackerHitPlane* b, double stripLength, EventStatistics &statistics ) const ;
    int calculatePointBetweenTwoLinesUsingVertex( const TVector3& pa, const TVector3& pb, const TVector3& pc, const TVector3& pd, const TVector3& vertex, TVector3& point) const ;
    // Shipped from CLHEP ...
    double cos2Theta(const TVector3 &p, const TVector3 & q) const ;
    
  protected:
    // processor parameters
    std::string                              _inputCollectionName {} ;
    std::string                              _inputRelCollectionName {} ;
    std::string                              _outputCollectionName {} ;
    std::string                              _outputRelCollectionName {} ;
    float                                    _nominalVertexX {0.f} ;
    float                                    _nominalVertexY {0.f} ;
    float                                    _nominalVertexZ {0.f} ;
    dd4hep::rec::Vector3D                    _nominalVertex {} ;
    float                                    _stripLengthTolerance {0.1f} ;
    double                                   _stripLength {0.f} ;
    std::string                              _subDetectorName {} ;
    const dd4hep::rec::SurfaceMap           *_surfaceMap {nullptr} ;
  };

  //--------------------------------------------------------------------------

  DDSpacePointBuilderProcessor::DDSpacePointBuilderProcessor() :
    Processor("DDSpacePointBuilder") {
    // modify processor description
    _description = "DDSpacePointBuilder combine si-strip measurements into 3D spacepoints (1TrackerHitPlanar+1TrackHitPlanar = 1 TrackerHit), that can be used by reconstruction" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection(EVENT::LCIO::TRACKERHIT,
                           "TrackerHitCollection",
                           "TrackerHitCollection",
                           _inputCollectionName,
                           std::string("FTDTrackerHits")); 

    registerInputCollection(EVENT::LCIO::LCRELATION,
                           "TrackerHitSimHitRelCollection",
                           "The name of the input collection of the relations of the TrackerHits to SimHits",
                           _inputRelCollectionName,
                           std::string("FTDTrackerHitRelations")); 

    registerOutputCollection(EVENT::LCIO::TRACKERHIT,
                            "SpacePointsCollection",
                            "SpacePointsCollection",
                            _outputCollectionName,
                            std::string("FTDSpacePoints"));

    registerOutputCollection(EVENT::LCIO::LCRELATION,
                            "SimHitSpacePointRelCollection",
                            "Name of the SpacePoint SimTrackerHit relation collection",
                            _outputRelCollectionName,
                            std::string("FTDSimHitSpacepointRelations"));
      
     
    registerProcessorParameter("NominalVertexX",
                               "The global x coordinate of the nominal vertex used for calculation of strip hit intersections",
                               _nominalVertexX,
                               float(0.0));


    registerProcessorParameter("NominalVertexY",
                               "The global x coordinate of the nominal vertex used for calculation of strip hit intersections",
                               _nominalVertexY,
                               float(0.0));


    registerProcessorParameter("NominalVertexZ",
                               "The global x coordinate of the nominal vertex used for calculation of strip hit intersections",
                               _nominalVertexZ,
                               float(0.0));
     
    // YV added
    registerProcessorParameter("StripLength",
                               "The length of the strips of the subdetector in mm",
                               _stripLength,
                               double(0.0));


    registerProcessorParameter("StriplengthTolerance",
                               "Tolerance added to the strip length when calculating strip hit intersections",
                               _stripLengthTolerance,
                               float(0.1));


    registerProcessorParameter( "SubDetectorName" , 
                               "Name of dub detector" ,
                               _subDetectorName ,
                                std::string("SIT") );
  }

  //--------------------------------------------------------------------------

  void DDSpacePointBuilderProcessor::init() {
    // usually a good idea to
    printParameters() ;
    _nominalVertex.fill( _nominalVertexX, _nominalVertexY, _nominalVertexZ ) ;
    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance() ;
    dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;
    dd4hep::DetElement det = theDetector.detector( _subDetectorName ) ;
    _surfaceMap = surfMan.map( det.name() ) ;
  }

  //--------------------------------------------------------------------------

  void DDSpacePointBuilderProcessor::processEvent( EVENT::LCEvent * evt ) {
    auto inputCollection = this->getCollection(   evt, _inputCollectionName    ) ;
    auto navigator = this->createNavigator( evt, _inputRelCollectionName ) ;
    if( (nullptr == inputCollection) or (nullptr == navigator) ) {
      return ;
    }
    auto spacePointCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::TRACKERHIT ) ;
    auto outputRelationCollection = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::LCRELATION ) ;
    // to store the weights
    IMPL::LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( EVENT::LCIO::LCREL_WEIGHTED ) ;
    outputRelationCollection->setFlag( lcFlag.getFlag() ) ;
    // local hit counters
    const unsigned int nHits = inputCollection->getNumberOfElements() ;
    EventStatistics statistics {} ;
    // store hits in map according to their CellID0
    std::map< int , std::vector< EVENT::TrackerHitPlane* > > cellID0HitMap {} ;
    for( unsigned int i=0 ; i<nHits ; i++ ) {
      auto trkHit = dynamic_cast<EVENT::TrackerHitPlane*>( inputCollection->getElementAt( i ) ) ;
      if( nullptr != trkHit ) {
        log<DEBUG3>() << "Add hit with CellID0 = " << trkHit->getCellID0() << " " << getCellID0Info( trkHit->getCellID0() ) << std::endl ;
        cellID0HitMap[ trkHit->getCellID0() ].push_back( trkHit ) ;
      }
    }
    UTIL::CellIDEncoder<IMPL::TrackerHitImpl> cellIDEncoder( UTIL::LCTrackerCellID::encoding_string() , spacePointCollection.get() ) ;
    const double stripLength = _stripLength * ( 1.0 + _stripLengthTolerance ) ;
    for( auto iter : cellID0HitMap ) {
      statistics._rawStripHits += iter.second.size() ;
      auto cellID0 = iter.first ;
      auto cellID0sBack = this->getCellID0sAtBack( cellID0 ) ;
      for( auto cellID0Back : cellID0sBack ) { 
        auto findIter = cellID0HitMap.find( cellID0Back ) ;
        if( cellID0HitMap.end() == findIter ) {
          continue ;
        }
        const auto nCombinations = iter.second.size() * findIter->second.size() ;
        log<DEBUG3>() << "strips: CellID0 " << cellID0  << " " << getCellID0Info( cellID0 )  << "(" << iter.second.size()
		      << " hits) <---> CellID0 " << cellID0Back << getCellID0Info( cellID0Back )
		      << "(" << findIter->second.size() << " hits)\n"
		      << "--> " << nCombinations << " possible combinations\n";
        statistics._possibleSpacePoints += nCombinations ;
        
        for( auto hitBack : iter.second ) {
          for( auto hitFront : findIter->second ) {
            auto& simHitsFront = navigator->getRelatedToObjects( hitFront ) ;
            auto& simHitsBack  = navigator->getRelatedToObjects( hitBack ) ;
            log<DEBUG3>() << "attempt to create space point from:" << std::endl ;
            log<DEBUG3>() << " front hit: " << hitFront << " no. of simhit = " << simHitsFront.size() ;
            if( not simHitsFront.empty() ) { 
              auto simhit = static_cast<const EVENT::SimTrackerHit*>( simHitsFront.at(0) ) ;
              log<DEBUG3>() << " first simhit = " << simhit << " mcp = "<< simhit->getMCParticle() << " ( " << simhit->getPosition()[0] << " " << simhit->getPosition()[1] << " " << simhit->getPosition()[2] << " ) " ; 
            }
            log<DEBUG3>() << std::endl;            
            log<DEBUG3>() << "  rear hit: " << hitBack << " no. of simhit = " << simHitsBack.size() ;
            if( not simHitsBack.empty() ) { 
              auto simhit = static_cast<const EVENT::SimTrackerHit*>( simHitsBack.at(0) ) ;
              log<DEBUG3>() << " first simhit = " << simhit << " mcp = "<< simhit->getMCParticle() << " ( " << simhit->getPosition()[0] << " " << simhit->getPosition()[1] << " " << simhit->getPosition()[2] << " ) " ; 
            }            
            log<DEBUG3>() << std::endl ;
            bool ghostHit = true ;
            if ( (simHitsFront.size() == 1) && (simHitsBack.size() == 1) ) {
              log<DEBUG3>() << "SpacePoint creation from two good hits:" << std::endl ;
              ghostHit = static_cast<EVENT::SimTrackerHit*>(simHitsFront[0])->getMCParticle() != static_cast<EVENT::SimTrackerHit*>(simHitsBack[0])->getMCParticle() ;
            }
            if ( ghostHit ) {
              log<DEBUG3>() << "SpacePoint Ghosthit!" << std::endl ;
            }
            auto spacePoint = this->createSpacePoint( hitFront, hitBack, stripLength, statistics ) ;
            if( nullptr == spacePoint ) {
              if ( ghostHit ) {
                log<DEBUG3>() << "Ghosthit correctly rejected" << std::endl ;
              }
              else {
                log<DEBUG3>() << "True hit rejected!" << std::endl ;
              }
              continue ;
            }
            cellIDEncoder.setValue( cellID0 ); //give the new hit, the CellID0 of the front hit
            cellIDEncoder.setCellID( spacePoint ) ;
            // store the hits it's composed of:
            spacePoint->rawHits().push_back( hitFront ) ;
            spacePoint->rawHits().push_back( hitBack ) ;
            spacePoint->setType( UTIL::set_bit( spacePoint->getType(), UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ) ) ;            
            spacePointCollection->addElement( spacePoint ) ; 
            statistics._createdSpacePoints++;
            ///////////////////////////////
            // make the relations
            if( simHitsFront.size() == 1 ) {              
              auto simHit = static_cast< EVENT::SimTrackerHit* >( simHitsFront[0] ) ;
              if( nullptr != simHit ) {
                auto rel = new IMPL::LCRelationImpl() ;
                rel->setFrom ( spacePoint ) ;
                rel->setTo  ( simHit ) ;
                rel->setWeight( 0.5 ) ;
                outputRelationCollection->addElement( rel ) ;
              }
            }            
            if( simHitsBack.size() == 1 ) {
              auto simHit = static_cast< EVENT::SimTrackerHit* >( simHitsBack[0] );
              
              if( nullptr != simHit ) {
                auto rel = new IMPL::LCRelationImpl() ;
                rel->setFrom ( spacePoint ) ;
                rel->setTo  ( simHit ) ;
                rel->setWeight( 0.5 ) ;
                outputRelationCollection->addElement( rel ) ;
              }
            }
          }  
        }
      }
    }
    // add collections
    evt->addCollection( spacePointCollection.release() ,     _outputCollectionName ) ;
    evt->addCollection( outputRelationCollection.release() , _outputRelCollectionName ) ;
    // log stats
    log<DEBUG3>() << "\n";
    log<DEBUG3>() << "Created " << statistics._createdSpacePoints << " space points ( raw strip hits: " << statistics._rawStripHits << ")\n";
    log<DEBUG3>() << "  There were " << statistics._rawStripHits << " strip hits available, giving " << statistics._possibleSpacePoints << " possible space points\n";
    log<DEBUG3>() << "  " << statistics._nStripsTooParallel << " space points couldn't be created, because the strips were too parallel\n";
    log<DEBUG3>() << "  " << statistics._nPlanesNotParallel << " space points couldn't be created, because the planes of the measurement surfaces where not parallel enough\n";
    log<DEBUG3>() << "  " << statistics._nOutOfBoundary     << " space points couldn't be created, because the result was outside the sensor boundary\n"; 
    log<DEBUG3>() << "\n";
  }

  //--------------------------------------------------------------------------

  EVENT::LCCollection* DDSpacePointBuilderProcessor::getCollection( EVENT::LCEvent* evt, const std::string name ) const {
    if( name.size() == 0 ) {
      return nullptr ;
    }
    try {
      return evt->getCollection( name ) ;
    }
    catch( EVENT::DataNotAvailableException& e ) {
      log<marlin::DEBUG2>() << "getCollection :  DataNotAvailableException : " << name <<  std::endl ;
      return nullptr ;
    }
  }
  
  //--------------------------------------------------------------------------

  std::shared_ptr<UTIL::LCRelationNavigator> DDSpacePointBuilderProcessor::createNavigator( EVENT::LCEvent* evt, const std::string name ) const {
    if( name.size() == 0 ) {
      return nullptr ;
    }
    try {
      auto collection = evt->getCollection( name ) ;
      return std::make_shared<UTIL::LCRelationNavigator>( collection ) ;
    }
    catch( EVENT::DataNotAvailableException& e ) {
      log<marlin::DEBUG2>() << "createNavigator :  DataNotAvailableException : " << name <<  std::endl ;
      return nullptr ;
    }
  }
  
  //--------------------------------------------------------------------------
  
  IMPL::TrackerHitImpl* DDSpacePointBuilderProcessor::createSpacePoint( EVENT::TrackerHitPlane* a , EVENT::TrackerHitPlane* b, double stripLength, EventStatistics &statistics ) const {
    const auto mmInverse = ( 1. / dd4hep::mm ) ;
    // point A
    dd4hep::rec::Vector3D positionA( a->getPosition()[0] * dd4hep::mm, a->getPosition()[1]* dd4hep::mm, a->getPosition()[2] * dd4hep::mm ) ;
    auto surfaceA = _surfaceMap->find( a->getCellID0() )->second ;
    auto normalA = mmInverse * surfaceA->normal().to<TVector3>() ;
    auto uA = mmInverse * surfaceA->u().to<TVector3>() ;
    auto vA = mmInverse * surfaceA->v().to<TVector3>() ;
    // point B
    dd4hep::rec::Vector3D positionB( b->getPosition()[0] * dd4hep::mm, b->getPosition()[1]* dd4hep::mm, b->getPosition()[2] * dd4hep::mm ) ;
    auto surfaceB = _surfaceMap->find( b->getCellID0() )->second ;
    auto normalB = mmInverse * surfaceB->normal().to<TVector3>() ;
    auto uB = mmInverse * surfaceB->u().to<TVector3>() ;
    auto vB = mmInverse * surfaceB->v().to<TVector3>() ;
    // First: check if the two measurement surfaces are parallel (i.e. the w are parallel or antiparallel)
    double angle = std::fabs( normalB.Angle( normalA ) ) ;
    static const double angleLimit = 1.*M_PI/180.;
    if( ( angle > angleLimit ) && ( angle < M_PI-angleLimit ) ) {
      statistics._nPlanesNotParallel++;
      log<DEBUG3>() << "\tThe planes of the measurement surfaces are not parallel enough, the angle between the W vectors is " << angle
      << " where the angle has to be smaller than " << angleLimit << " or bigger than " << M_PI-angleLimit << "\n\n";
      return nullptr ;
    }
    // Next: check if the angle between the strips is not 0
    angle = std::fabs( vB.Angle( vA ) ) ;
    if(( angle < angleLimit )||( angle > M_PI-angleLimit )) {      
      statistics._nStripsTooParallel++;
      log<DEBUG3>() << "\tThe strips (V vectors) of the measurement surfaces are too parallel, the angle between the V vectors is " << angle
      << " where the angle has to be between " << angleLimit << " or bigger than " << M_PI-angleLimit << "\n\n";
      return nullptr ;
    }
    // Next we want to calculate the crossing point.
    auto ddLocalDirA = surfaceA->globalToLocal( positionA ) ;
    auto ddLocalDirB = surfaceB->globalToLocal( positionB ) ;
    dd4hep::rec::Vector2D ddStartVecA( ddLocalDirA.u(), (-stripLength * dd4hep::mm)/2.0 ) ;
    dd4hep::rec::Vector2D ddEndVecA(   ddLocalDirA.u(),  (stripLength * dd4hep::mm)/2.0 ) ;
    dd4hep::rec::Vector2D ddStartVecB( ddLocalDirB.u(), (-stripLength * dd4hep::mm)/2.0 ) ;
    dd4hep::rec::Vector2D ddEndVecB(   ddLocalDirB.u(),  (stripLength * dd4hep::mm)/2.0 ) ;
    auto startPositionA = mmInverse * surfaceA->localToGlobal( ddStartVecA ).to<TVector3>() ;
    auto endPositionA   = mmInverse * surfaceA->localToGlobal( ddEndVecA ).to<TVector3>() ;
    auto startPositionB = mmInverse * surfaceB->localToGlobal( ddStartVecB ).to<TVector3>() ;
    auto endPositionB   = mmInverse * surfaceB->localToGlobal( ddEndVecB ).to<TVector3>() ;

    TVector3 point(0, 0, 0), vertex(0, 0, 0) ;
    const auto validIntersection = calculatePointBetweenTwoLinesUsingVertex( startPositionA, endPositionA, startPositionB, endPositionB, vertex, point ) ;

    if ( validIntersection != 0 ) {
      log<DEBUG3>() << "\tNo valid intersection for lines" << std::endl ;
      return nullptr ;
    }
    log<DEBUG3>() << "\tVertex: Position of space point (global) : ( " << point.x() << " " << point.y() << " " << point.z() << " )\n" ;
    // using dd4hep to check if hit within boundaries
    dd4hep::rec::Vector3D ddPoint( point.x() * dd4hep::mm, point.y() * dd4hep::mm, point.z() * dd4hep::mm ) ;
    if ( ! surfaceA->insideBounds( ddPoint ) ) {
      statistics._nOutOfBoundary++ ;
      log<DEBUG3>() << " SpacePoint position lies outside the boundary of the layer " << std::endl ;
      return nullptr ;
    }
    // create the new TrackerHit
    auto spacePoint = new IMPL::TrackerHitImpl() ;
    double pos[3] = {point.x(), point.y(), point.z() } ;
    spacePoint->setPosition(  pos  ) ;
    // set error treating the strips as stereo with equal and opposite rotation -- for reference see Karimaki NIM A 374 p367-370
    // first calculate the covariance matrix in the cartisian coordinate system defined by the sensor 
    // here we assume that du is the same for both sides
    if( std::fabs( a->getdU() - b->getdU() ) > 1.0e-06 ) {
      log<ERROR>() << "\tThe measurement errors of the two 1D hits must be equal \n\n" ;
      // measurement errors are not equal don't create a spacepoint
      return nullptr ; 
    }
    const double du2 = a->getdU() * a->getdU() ;
  
    // rotate the strip system back to double-layer wafer system
    auto uSensor = uA + uB ;
    auto vSensor = vA + vB ;
    auto wSensor = normalA + normalB ;
    
    TRotation rotationSensor ;
    rotationSensor.RotateAxes( uSensor, vSensor, wSensor ) ;
    double rotationArray[9] = { 
      rotationSensor.XX(), rotationSensor.XY(), rotationSensor.XZ(),
      rotationSensor.YX(), rotationSensor.YY(), rotationSensor.YZ(), 
      rotationSensor.ZX(), rotationSensor.ZY(), rotationSensor.ZZ() 
    } ;
    TMatrixD rotationSensorMatrix ( 3, 3, &rotationArray[0] ) ;
    double cos2Alpha = this->cos2Theta( vA, vSensor ) ; // alpha = strip angle   
    double sin2Alpha = 1 - cos2Alpha ; 

    TMatrixDSym covariancePlane( 3 ) ; // u,v,w
    covariancePlane( 1, 1 ) = ( 0.5 * du2 ) / cos2Alpha;
    covariancePlane( 2, 2 ) = ( 0.5 * du2 ) / sin2Alpha;
    auto covarianceXYZ = covariancePlane.Similarity( rotationSensorMatrix ) ;
    
    std::vector<float> covariance( 9 ) ; 
    int icov = 0 ;
    
    for(int irow=0; irow<3; ++irow ) {
      for(int jcol=0; jcol<irow+1; ++jcol) {
        covariance[icov] = covarianceXYZ[irow][jcol] ;
        ++icov ;
      }
    }
    spacePoint->setCovMatrix( covariance ) ;  
    const auto pointTime = std::min( a->getTime(), b->getTime() ) ;
    spacePoint->setTime( pointTime ) ;
    log<DEBUG3>() << "\tHit accepted " << std::endl << std::endl ;
    return spacePoint ;    
  }
  
  //--------------------------------------------------------------------------
  
  std::vector<int> DDSpacePointBuilderProcessor::getCellID0sAtBack( int cellID0 ) const {
    std::vector<int> back {} ;  
    // find out detector, layer
    UTIL::BitField64 cellID( UTIL::LCTrackerCellID::encoding_string() );
    cellID.setValue( cellID0 ) ;
    int subdet = cellID[ UTIL::LCTrackerCellID::subdet() ] ;
    int layer  = cellID[ UTIL::LCTrackerCellID::layer() ] ;
    if ( subdet != UTIL::ILDDetID::FTD ) {
      // check if sensor is in front
      // even layers are front sensors
      if( layer%2 == 0 ) { 
        cellID[ UTIL::LCTrackerCellID::layer() ] = layer + 1 ;
        // it is assumed that the even layers are the front layers
        // and the following odd ones the back layers        
        back.push_back( cellID.lowWord() ) ;
      }
    }
    else {
      dd4hep::Detector & theDetector2 = dd4hep::Detector::getInstance() ;
      dd4hep::DetElement ftdDE = theDetector2.detector( _subDetectorName ) ;
      dd4hep::rec::ZDiskPetalsData* ft = ftdDE.extension<dd4hep::rec::ZDiskPetalsData>() ;
      int sensor = cellID[ UTIL::LCTrackerCellID::sensor() ] ;
      int Nsensors = ft->layers.at(layer).sensorsPerPetal ;
      log<DEBUG3>() << " layer " << layer << " sensors " << Nsensors << std::endl ;
      log<DEBUG3>() << " so sensor " << sensor << " is connected with sensor " << sensor + Nsensors/2 << std::endl ;
      // check if sensor is in front
      if (sensor <= Nsensors / 2 ) {
        cellID[ UTIL::LCTrackerCellID::sensor() ] = sensor + Nsensors / 2 ; 
        // it is assumed, that sensors 1 until n/2 will be on front
        // and sensor n/2 + 1 until n are at the back
        // so the sensor x, will have sensor x+n/2 at the back
        back.push_back( cellID.lowWord() ) ;
      }
    }
    return back ;
  }
  
  //--------------------------------------------------------------------------
  
  std::string DDSpacePointBuilderProcessor::getCellID0Info( int cellID0 ) const {
    std::stringstream ss ;
    //find out layer, module, sensor
    UTIL::BitField64  cellID( UTIL::LCTrackerCellID::encoding_string() ) ;
    cellID.setValue( cellID0 ) ;
    int subdet = cellID[ UTIL::LCTrackerCellID::subdet() ] ;
    int side   = cellID[ UTIL::LCTrackerCellID::side() ] ;
    int module = cellID[ UTIL::LCTrackerCellID::module() ] ;
    int sensor = cellID[ UTIL::LCTrackerCellID::sensor() ] ;
    int layer  = cellID[ UTIL::LCTrackerCellID::layer() ] ;    
    ss << "(su" << subdet << ",si" << side << ",la" << layer << ",mo" << module << ",se" << sensor << ")" ;
    return ss.str() ;
  }
  
  //--------------------------------------------------------------------------
  
  int DDSpacePointBuilderProcessor::calculateCrossingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y ) const {
    float a = (x1*ey1 - y1*ex1) - (x2*ey1 - y2*ex1) ;
    float b = ex2*ey1 - ex1*ey2 ;
    // if b==0 the two directions e1 and e2 are parallel and there is no crossing!
    if( std::fabs(b) < DDSpacePointBuilderProcessor::crossingPointEpsilon ) {
      return 1 ;
    }
    float t = a / b ;
    x = x2 + t*ex2 ;
    y = y2 + t*ey2 ;
    return 0 ;
  }
  
  //--------------------------------------------------------------------------
  
  int DDSpacePointBuilderProcessor::calculatePointBetweenTwoLines( 
    const PositionXYZ& P1, 
    const PositionXYZ& V1,
    const PositionXYZ& P2, 
    const PositionXYZ& V2, 
    PositionXYZ& point ) const {
    // Richgungsvektor normal auf die anderen beiden:
    auto n = V1.Cross( VectorXYZ(V2) );    
    // Now we want to rotate into a coordinate system, where n is parallel to the z axis
    // For this: first set phi to 0
    // then: set theta to 0 (we set phi to 0 first, so we can then rotate arount the y axis)
    RotationXYZ rot ;
    rot.SetPhi( -n.phi() ) ;
    // rot.rotateZ( -n.phi() );
    PositionXYZ nPrime = rot * n ;
    // CLHEP::Hep3Vector nPrime = rot * n; //now the phi of nPrime should be 0
    log<DEBUG0>() << "phi of n' = " << nPrime.phi() << " (it should be 0!!!)\n";
    rot.SetTheta( -n.theta() ) ;
    // rot.rotateY( -n.theta() );
    nPrime = rot * n ;
    log<DEBUG0>() << "phi of n'' = " << nPrime.phi() << " (it should be 0!!!)\n";
    log<DEBUG0>() << "theta of n'' = " << nPrime.theta() <<  " (it should be 0!!!)\n";
    // Now rotate all the vectors and points into this coordinatesystem.
    PositionXYZ P1prime = rot * P1;
    PositionXYZ V1prime = rot * V1;
    PositionXYZ P2prime = rot * P2;
    PositionXYZ V2prime = rot * V2;    
    // What is the gain of rotating into this system?
    double x(0), y(0) ;
    int res = calculateCrossingPoint( P1prime.x(), P1prime.y(), V1prime.x(), V1prime.y(), P2prime.x(), P2prime.y(), V2prime.x(), V2prime.y(), x, y );
    if ( res != 0 ) {
      return res ;
    }
    point.SetX( x );
    point.SetY( y );
    point.SetZ( (P1prime.z() + P2prime.z())/2. );
    // Now transform back to the global coordinates
    point = rot.Inverse() * point ;
    return 0;
  }
  
  //--------------------------------------------------------------------------
  
  int DDSpacePointBuilderProcessor::calculatePointBetweenTwoLinesUsingVertex( 
                                                const TVector3& pa, 
                                                const TVector3& pb, 
                                                const TVector3& pc, 
                                                const TVector3& pd,
                                                const TVector3& vertex,
                                                TVector3& point) const {
    // A general point on the line joining point PA to point PB is
    // x, where 2*x=(1+m)*PA + (1-m)*PB. Similarly for 2*y=(1+n)*PC + (1-n)*PD.
    // Suppose that v is the vertex. Requiring that the two 'general
    // points' lie on a straight through v means that the vector x-v is a 
    // multiple of y-v. This condition fixes the parameters m and n.
    // We then return the 'space-point' x, supposed to be the layer containing PA and PB. 
    // We require that -1<m<1, otherwise x lies 
    // outside the segment PA to PB; and similarly for n.
    bool ok = true ;    
    TVector3 vab( pa - pb ) ;
    TVector3 vcd( pc - pd ) ;
    TVector3 s( pa + pb - 2*vertex ) ;  // twice the vector from vertex to midpoint
    TVector3 t( pc + pd - 2*vertex ) ;  // twice the vector from vertex to midpoint
    TVector3 qs ( vab.Cross( s ) ) ;
    TVector3 rt ( vcd.Cross( t ) ) ;
    double m = ( -(s*rt) / (vab*rt) ) ; // ratio for first line
    const double limit = 1.0 ;
    
    if ( m > limit || m < -limit) {
      ok = false ;
    }
    else {      
      double n = ( - ( t * qs ) / ( vcd * qs ) ) ; // ratio for second line
  	  if ( n > limit || n < -limit) {    
        ok = false ;
      }
    }
    
    if ( ok ) {
      point = 0.5 * ( pa + pb + m * vab ) ;
    }    
    return ok ? 0 : 1 ;
  }

  //--------------------------------------------------------------------------
  
  double DDSpacePointBuilderProcessor::cos2Theta(const TVector3 &p, const TVector3 & q) const {
    double arg ;
    double ptot2 = p.Mag2() ;
    double qtot2 = q.Mag2() ;
    if ( ptot2 == 0 || qtot2 == 0 ) {
      arg = 1.0 ;
    }
    else {
      double pdq = p.Dot( q ) ;
      arg = ( pdq / ptot2 ) * ( pdq / qtot2 ) ;
      // More naive methods overflow on vectors which can be squared
      // but can't be raised to the 4th power.
      if(arg >  1.0) arg =  1.0 ;
   }
   return arg ;
  }

  //--------------------------------------------------------------------------

  // processor declaration
  DDSpacePointBuilderProcessor aDDSpacePointBuilderProcessor ;
}
