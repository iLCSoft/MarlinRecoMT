#include <MarlinRecoMT/OverlayFileHandler.h>

// -- marlin headers
#include <marlin/Logging.h>
using namespace marlin::loglevel ;

// -- lcio headers
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>

namespace marlinreco_mt {
  
  /// Set the lcio file name
  void OverlayFileHandler::setFileName(const std::string& fname) {
    _fileName = fname ;
  }

  //--------------------------------------------------------------------------

  /// Get the number of events available in the file
  unsigned int OverlayFileHandler::getNumberOfEvents() {
    openFile() ;
    return _lcReader->getNumberOfEvents() ;
  }

  //--------------------------------------------------------------------------

  /// Get the event number at the specified index (look in the event map)
  unsigned int OverlayFileHandler::getEventNumber(unsigned int index) {
    openFile() ;
    return _eventMap.at( index * 2 + 1 ) ;
  }

  //--------------------------------------------------------------------------

  /// Get the run number at the specified index (look in the event map)
  unsigned int OverlayFileHandler::getRunNumber(unsigned int index) {
    openFile() ;
    return _eventMap.at( index * 2 ) ;
  }

  //--------------------------------------------------------------------------

  /// Read the specified event, by run and event number
  std::shared_ptr<EVENT::LCEvent> OverlayFileHandler::readEvent(int runNumber, int eventNumber) {
    openFile() ;
    streamlog_out( DEBUG6 ) << "*** Reading event from file : '" << _fileName 
          << "',  event number " << eventNumber << " of run " << runNumber << "." << std::endl ;
    return _lcReader->readEvent( runNumber, eventNumber, EVENT::LCIO::UPDATE ) ;
  }

  //--------------------------------------------------------------------------

  /// Proxy method to open the LCIO file
  void OverlayFileHandler::openFile() {
    if(nullptr == _lcReader) {
      _lcReader = std::make_shared<FileReader>( MT::LCReader::directAccess ) ;
      streamlog_out( MESSAGE ) << "*** Opening file for overlay, file name:" << _fileName << std::endl ;
      _lcReader->open( _fileName ) ;
      _lcReader->getEvents( _eventMap ) ;
      streamlog_out( MESSAGE ) << "*** Opening file for overlay : number of available events: " << _lcReader->getNumberOfEvents() << std::endl ;
    }
  }

}
