// -- marlin headers
#include <marlin/Processor.h>
#include <marlin/ProcessorApi.h>
#include <marlin/PluginManager.h>
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- lcio headers
#include <IMPL/LCFlagImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>

// -- dd4hep headers
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

// -- marlinreco mt headers
#include "MarlinRecoMT/CalorimeterHitType.h"

namespace marlinreco_mt {
  
  class BruteForceEcalGapFiller : public marlin::Processor {
  public:
    static constexpr unsigned int MAXMODULE = 10 ; 
    static constexpr unsigned int MAXSTAVE = 15 ; 
    static constexpr unsigned int MAXLAYER = 50 ;
    // don't consider differences below this distance to be a gap
    static constexpr float DISTANCELIMIT = 0.01 ;
    // flexibility, as ratio
    static constexpr float SLOPDELTA = 0.01 ;
    typedef std::vector<EVENT::CalorimeterHit*> HitMapping[MAXLAYER][MAXSTAVE][MAXMODULE] ;

  public:
    BruteForceEcalGapFiller() ;
    BruteForceEcalGapFiller(const BruteForceEcalGapFiller &) = delete ;
    BruteForceEcalGapFiller &operator=(const BruteForceEcalGapFiller &) = delete ;
    
    // from marlin processor
    void init() ;
    void processEvent( EVENT::LCEvent * evt ) ;
    
  private:
    dd4hep::rec::LayeredCalorimeterData *getGeometryData( const int ihitType ) const ;
    void fillHitMap( EVENT::LCCollection *collection, HitMapping &hitMap ) const ;
    void addIntraModuleGapHits( EVENT::LCCollection* newcol, const HitMapping &hitMap, dd4hep::rec::LayeredCalorimeterData *calodata ) const ;
    void addInterModuleGapHits( EVENT::LCCollection* newcol, const HitMapping &hitMap, dd4hep::rec::LayeredCalorimeterData *calodata ) const ;
    
  private:
    marlin::InputCollectionProperty _inputHitCollection {this, EVENT::LCIO::CALORIMETERHIT, "inputHitCollection" ,
                             "input simcalhit Collection Name" } ;

    marlin::OutputCollectionProperty _outputHitCollection {this, EVENT::LCIO::CALORIMETERHIT, "outputHitCollection",
                             "output calorimeterhit Collection Name" } ;

    marlin::Property<std::string> _cellIDLayerString {this, "CellIDLayerString" ,
                               "name of the part of the cellID that holds the layer", "layer" } ;

    marlin::Property<std::string> _cellIDModuleString {this, "CellIDModuleString" ,
                               "name of the part of the cellID that holds the module", "module" } ;

    marlin::Property<std::string> _cellIDStaveString {this, "CellIDStaveString" ,
                               "name of the part of the cellID that holds the stave", "stave" } ;

    marlin::Property<float> _interModuleDist {this, "expectedInterModuleDistance",
             "size of gap across module boundaries (from edge to edge of cells, in mm ; accuracy < cell size)", 7. } ;

    marlin::Property<float>  _interModuleNonlinearFactor {this, "interModuleNonlinearFactor",
             "nonlin factor f: E_corr = interModuleCorrectionFactor*(1/f)*log(1 + f*E_calc)", 1. } ;

    marlin::Property<float> _intraModuleNonlinearFactor {this, "intraModuleNonlinearFactor",
             "nonlin factor f: E_corr = intraModuleCorrectionFactor*(1/f)*log(1 + f*E_calc)", 1. } ;

    marlin::Property<float> _interModuleFactor {this, "interModuleCorrectionFactor",
             "factor applied to calculated energy of inter-module gap hits", 0.35 } ;

    marlin::Property<float> _intraModuleFactor {this, "intraModuleCorrectionFactor",
             "factor applied to calculated energy of intra-module gap hits", 1.0 } ;

  private:
    dd4hep::rec::LayeredCalorimeterData    *_barrelGeometry {nullptr} ;
    dd4hep::rec::LayeredCalorimeterData    *_endcapGeometry {nullptr} ;
  };
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  BruteForceEcalGapFiller::BruteForceEcalGapFiller( )  : 
    marlin::Processor( "BruteForceEcalGapFiller" ) {
    _description = "makes a collection of ECAL gap hits" ;
  }

  //--------------------------------------------------------------------------

  void BruteForceEcalGapFiller::init() {
    printParameters();

    // find detector data
    dd4hep::Detector &detector = dd4hep::Detector::getInstance() ;
    
    const std::vector< dd4hep::DetElement>& barrelDetectors = dd4hep::DetectorSelector(detector).detectors( 
      ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), 
      ( dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD ) ) ;
      
    const std::vector< dd4hep::DetElement>& endcapDetectors = dd4hep::DetectorSelector(detector).detectors( 
      ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP), 
      ( dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD ) ) ;

    if( (barrelDetectors.size() == 1) ) {
      _barrelGeometry = barrelDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>() ;
    }
    if( ( endcapDetectors.size() == 1 ) ) {
      _endcapGeometry = endcapDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();  
    }
    if( (nullptr == _barrelGeometry) and (nullptr == _endcapGeometry) ) {
      marlin::ProcessorApi::abort( this, "Couldn't find any of the ecal calorimeters (endcap and barrel) !" ) ;
    }
    if( nullptr == _barrelGeometry ) {
      log<WARNING>() << "ECal barrel calorimeter data not found !" << std::endl ;
    }
    if( nullptr == _endcapGeometry ) {
      log<WARNING>() << "ECal endcap calorimeter data not found !" << std::endl ;
    }
  }

  //--------------------------------------------------------------------------

  void BruteForceEcalGapFiller::processEvent( EVENT::LCEvent * evt ) {
    log<DEBUG3>() << "looking for collection: " << _inputHitCollection << std::endl ;
    try {
      auto col = evt->getCollection( _inputHitCollection ) ;
      int numElements = col->getNumberOfElements();
      log<DEBUG3>() << _inputHitCollection << " number of elements = " << numElements << std::endl ;
      if( numElements == 0 ) {
        return ;
      }
      // get the correct geometry data 
      auto hittype = static_cast<EVENT::CalorimeterHit*>( col->getElementAt( 0 ) )->getType() ;
      dd4hep::rec::LayeredCalorimeterData *caloData = getGeometryData( hittype ) ;
      
      // fill the hit map
      HitMapping hitMap ;
      fillHitMap( col, hitMap ) ;

      // create new collection: hits
      std::string encodingString = col->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ) ;
      IMPL::LCFlagImpl flag ;
      flag.setBit( EVENT::LCIO::CHBIT_LONG ) ;
      flag.setBit( EVENT::LCIO::RCHBIT_TIME ) ; // store timing on output hits.
      auto newcol = std::make_unique<IMPL::LCCollectionVec>( EVENT::LCIO::CALORIMETERHIT ) ;
      newcol->parameters().setValue( EVENT::LCIO::CellIDEncoding, encodingString ) ;
      newcol->setFlag( flag.getFlag() ) ;

      // now make the gap hits
      addIntraModuleGapHits( newcol.get(), hitMap, caloData ) ; // gaps within a module
      addInterModuleGapHits( newcol.get(), hitMap, caloData ) ; // gaps between modules

      evt->addCollection( newcol.release(), _outputHitCollection ) ;
    } 
    catch( EVENT::DataNotAvailableException &e ) {
      log<DEBUG3>() << "could not find input collection " << _inputHitCollection << std::endl ;
    }
  }
  
  //--------------------------------------------------------------------------

  dd4hep::rec::LayeredCalorimeterData *BruteForceEcalGapFiller::getGeometryData( const int ihitType ) const {
    // get information about geometry
    // calorimeter hit type used to decide if it's in barrel or endcap
    CHT calHitType( ihitType ) ;
    dd4hep::rec::LayeredCalorimeterData *caloData {nullptr} ;
    if ( calHitType.is( CHT::barrel ) ) {
      caloData = _barrelGeometry ;
    }
    else if ( calHitType.is( CHT::endcap ) ) {
      caloData = _endcapGeometry;
    }
    else {
      log<WARNING>() << "Input ecal hit collection is neither barrel nor endcap" << std::endl ;
      marlin::ProcessorApi::skipCurrentEvent( this ) ;
    }
    if( nullptr == caloData ) {
      log<WARNING>() << "No calorimeter data found for the ecal input hit collection! Please check your settings!" << std::endl ;
      marlin::ProcessorApi::skipCurrentEvent( this ) ;
    }
    return caloData ;
  }
  
  //--------------------------------------------------------------------------
  
  void BruteForceEcalGapFiller::fillHitMap( EVENT::LCCollection *collection, HitMapping &hitMap ) const {
    UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder( collection ) ;
    auto numElements = collection->getNumberOfElements() ;
    // loop over input hits
    for (int j(0) ; j < numElements ; ++j ) {
      auto hit = static_cast<EVENT::CalorimeterHit*>( collection->getElementAt( j ) ) ;
      unsigned int layer  = idDecoder( hit )[ _cellIDLayerString ] ;
      unsigned int stave  = idDecoder( hit )[ _cellIDStaveString ] ;
      unsigned int module = idDecoder( hit )[ _cellIDModuleString ] ;
      if( (layer >= MAXLAYER) or (stave >= MAXSTAVE) or (module >= MAXMODULE) ) {
        marlin::ProcessorApi::abort( this, "Hit with incorrect layer, module or stave number!" ) ;
      }
      hitMap [ layer ][ stave ][ module ].push_back( hit ) ;
    }
  }
  
  //--------------------------------------------------------------------------

  void BruteForceEcalGapFiller::addIntraModuleGapHits( EVENT::LCCollection* newcol, const HitMapping &hitMap, dd4hep::rec::LayeredCalorimeterData *calodata ) const {
    // look for gaps within modules
    // i.e. between wafers, between towers
    log<DEBUG3>() << " starting addIntraModuleGapHits" << std::endl ;
    for (unsigned int il=0; il<MAXLAYER; il++) {
      // we have to get the cell sizes here
      const float cellsizeA = calodata->layers[il].cellSize0 / dd4hep::mm ;
      const float cellsizeB = calodata->layers[il].cellSize1 / dd4hep::mm ;
      log<DEBUG0>() << "cell sizes in layer " << il << " = " << cellsizeA << " " << cellsizeB << " mm" << std::endl ;
      for (unsigned int is=0; is<MAXSTAVE; is++) {
        for (unsigned int im=0; im<MAXMODULE; im++) {
        	auto &theseHits = hitMap [ il ][ is ][ im ] ;
        	if ( not theseHits.empty() ) {
        	  bool gap(false);
        	  float enFrac(0);
        	  for ( unsigned int ih=0; ih<theseHits.size()-1; ih++) {
        	    for ( unsigned int jh=ih+1; jh<theseHits.size(); jh++) {
        	      float dist1d[3] = {0} ;
        	      for (int i=0; i<3; i++) {
        		      dist1d[i] = std::fabs( theseHits[ih]->getPosition()[i] - theseHits[jh]->getPosition()[i] );
                }
        	      float distXY = std::sqrt( dist1d[0]*dist1d[0] + dist1d[1]*dist1d[1] ) ;
        	      gap = false ;
        	      if (calodata == _barrelGeometry) {
              		if ( dist1d[2]<DISTANCELIMIT && // same z coord
              		     distXY>(1.+SLOPDELTA)*cellsizeA && // bigger than one cell period, smaller than two
              		     distXY<(2.-SLOPDELTA)*cellsizeA ) {
              		  gap = true;
              		  enFrac = (distXY-cellsizeA)/cellsizeA;
              		} 
                  else if (distXY<DISTANCELIMIT && // same x-y coord 
              			   dist1d[2]>(1.+SLOPDELTA)*cellsizeB && 
              			   dist1d[2]<(2.-SLOPDELTA)*cellsizeB ) {
              		  gap = true;
              		  enFrac = (dist1d[2]-cellsizeB)/cellsizeB;
              		}
        	      } 
                else { // endcap
              		if ( dist1d[1]<DISTANCELIMIT &&
              		     dist1d[0]>(1.+SLOPDELTA)*cellsizeA &&
              		     dist1d[0]<(2.-SLOPDELTA)*cellsizeA ) { // be careful, if different size in x,y may have to worry about stave
              		  gap = true;
              		  enFrac = (dist1d[0]-cellsizeA)/cellsizeA;
              		} else if ( dist1d[0]<DISTANCELIMIT &&
              			    dist1d[1]>(1.+SLOPDELTA)*cellsizeB &&
              			    dist1d[1]<(2.-SLOPDELTA)*cellsizeB ) { // be careful, if different size in x,y may have to worry about stave
              		  gap = true;
              		  enFrac = (dist1d[1]-cellsizeB)/cellsizeB;
              		}
        	      }
        	      if ( gap ) {
              		log<DEBUG0>() << " GOT A GAP " << std::endl ;
              		float position[3]={0.};
              		for (int k=0; k<3; k++) {
              		  position[k] = 0.5*(theseHits[ih]->getPosition()[k] + theseHits[jh]->getPosition()[k]);
                  }
              		float extraEnergy = enFrac*(theseHits[ih]->getEnergy() + theseHits[jh]->getEnergy())/2.;
              		float mintime = std::min( theseHits[ih]->getTime(), theseHits[jh]->getTime() );
              		CHT::CaloType cht_type = CHT::em;
              		CHT::CaloID   cht_id   = CHT::ecal;
              		CHT::Layout   cht_lay  = (calodata == _barrelGeometry) ? CHT::barrel : CHT::endcap ;
              		auto newGapHit = new IMPL::CalorimeterHitImpl() ;
              		newGapHit->setEnergy( _intraModuleFactor* std::log ( 1 + _intraModuleNonlinearFactor*extraEnergy )/_intraModuleNonlinearFactor );
              		newGapHit->setPosition( position );
              		newGapHit->setTime( mintime );
              		newGapHit->setType( CHT( cht_type , cht_id , cht_lay , il) );
              		newcol->addElement( newGapHit );
        	      } // if gap
        	    } // jh
        	  } // ih
        	} // >1 hit
        } // im
      } // is
    } // ilayer
    log<DEBUG0>() << " done addIntraModuleGapHits " << newcol->getNumberOfElements() << std::endl ;
  }
  
  //--------------------------------------------------------------------------

  void BruteForceEcalGapFiller::addInterModuleGapHits( EVENT::LCCollection* newcol, const HitMapping &hitMap, dd4hep::rec::LayeredCalorimeterData *calodata ) const {
    // look for gaps between modules
    //  compare hits in same stave, same layer
    log<DEBUG3>() << " starting addInterModuleGapHits" << std::endl ;
    for (unsigned int il=0; il<MAXLAYER; il++) {
      // we have to get the cell sizes here
      const float cellsizeA = calodata->layers[il].cellSize0 / dd4hep::mm ;
      const float cellsizeB = calodata->layers[il].cellSize1 / dd4hep::mm ;      
      for (unsigned int is=0; is<MAXSTAVE; is++) {
        for (unsigned int im=0; im<MAXMODULE; im++) {
        	auto &theseHits = hitMap [ il ][ is ][ im ] ;
        	if ( theseHits.empty() ) {
            continue ;
          }
        	// look in next module
        	if ( im+1<MAXMODULE ) {
        	  auto &nextHits = hitMap [ il ][ is ][ im+1 ] ;
        	  if ( nextHits.empty() ) {
              continue ;
            }
        	  bool gap(false);
        	  float enFrac(0);
        	  for ( unsigned int ih=0; ih<theseHits.size(); ih++) {
        	    for ( unsigned int jh=0; jh<nextHits.size(); jh++) {
        	      float dist1d[3] = {0} ;
        	      for (int i=0; i<3; i++) {
        		        dist1d[i] = std::fabs( theseHits[ih]->getPosition()[i] - nextHits[jh]->getPosition()[i] );                  
                }

        	      float distXY = std::sqrt( dist1d[0]*dist1d[0] + dist1d[1]*dist1d[1] ) ;
        	      gap = false ;
        	      if (calodata == _barrelGeometry) { // intermodule gaps only along z
              		if ( distXY<DISTANCELIMIT && // same phi coord
              		     dist1d[2] < _interModuleDist + cellsizeB*1.9 ) { // _interModuleDist is expected distance between sensor edges
              		  gap = true;
                    enFrac = dist1d[2] / cellsizeB;
              		}
        	      } 
                else { // endcap
              		if ( dist1d[1]<DISTANCELIMIT && // same y
              		     dist1d[0] < _interModuleDist + 1.9*cellsizeA ) { // be careful, if different size in x,y may have to worrk about stave
              		  gap = true;
              		  enFrac = dist1d[0]/cellsizeA;
              		} else if ( dist1d[0]<DISTANCELIMIT && // same x
              			    dist1d[1] < _interModuleDist + 1.9*cellsizeB ) { // be careful, if different size in x,y may have to worrk about stave
              		  gap = true;
              		  enFrac = dist1d[1]/cellsizeB;
              		}
        	      }
        	      if ( gap ) {
              		log<DEBUG0>() << " addInterModuleGapHits: found gap " << dist1d[0] << " " << dist1d[1] << " " << dist1d[2] << std::endl ;
              		float position[3]={0.};
              		for (int k=0; k<3; k++) {
              		  position[k] = 0.5*(theseHits[ih]->getPosition()[k] + nextHits[jh]->getPosition()[k]);
                  }
              		float extraEnergy = enFrac*(theseHits[ih]->getEnergy() + nextHits[jh]->getEnergy())/2.;
              		float mintime = std::min( theseHits[ih]->getTime(), nextHits[jh]->getTime() );
              		CHT::CaloType cht_type = CHT::em;
              		CHT::CaloID   cht_id   = CHT::ecal;
              		CHT::Layout   cht_lay  = (calodata == _barrelGeometry) ? CHT::barrel : CHT::endcap ;
              		auto newGapHit = new IMPL::CalorimeterHitImpl();
              		newGapHit->setEnergy( _interModuleFactor* std::log ( 1 + _interModuleNonlinearFactor*extraEnergy )/_interModuleNonlinearFactor );
              		newGapHit->setPosition( position );
              		newGapHit->setTime( mintime );
              		newGapHit->setType( CHT( cht_type , cht_id , cht_lay , il) );
              		newcol->addElement( newGapHit );
        	      } // if gap
        	    } // jh
        	  } // ih
        	} // >1 hit
        } // im
      } // is
    } // ilayer
    log<DEBUG3>() << " done addInterModuleGapHits " << newcol->getNumberOfElements() << std::endl ;
  }

  MARLIN_DECLARE_PROCESSOR( BruteForceEcalGapFiller )
}
