#ifndef MARLINRECOMT_OVERLAYFILEHANDLER_H
#define MARLINRECOMT_OVERLAYFILEHANDLER_H 1

// -- std headers
#include <string>
#include <vector>
#include <memory>

// -- lcio headers
#include <EVENT/LCEvent.h>
#include <MT/LCReader.h>

namespace marlinreco_mt {

  /**
   *  @brief  OverlayFileHandler class
   */
  class OverlayFileHandler {
    using FileReader = MT::LCReader ;
    
  public:
    // Default constructors
    OverlayFileHandler() = default ;
    OverlayFileHandler(const OverlayFileHandler&) = default ;
    OverlayFileHandler& operator =(const OverlayFileHandler&) = default ;

    /**
     *  @brief  Set the LCIO file name
     *  
     *  @param  fname the name of the LCIO file
     */
    void setFileName(const std::string& fname) ;
    
    /**
     *  @brief  Get the number of events available in the file
     */
    unsigned int getNumberOfEvents() ;
    
    /**
     *  @brief  Get the event number at the specified index (look in the event map)
     *  
     *  @param  index the nth event to get
     */
    unsigned int getEventNumber(unsigned int index) ;
    
    /**
     *  @brief  Get the run number at the specified index (look in the event map)
     *  
     *  @param  index the nth run to get
     */
    unsigned int getRunNumber(unsigned int index) ;
    
    /**
     *  @brief  Read the specified event, by run and event number
     *  
     *  @param  runNumber the run number of the event to read
     *  @param  eventNumber the event number of the event to read
     */
    std::shared_ptr<EVENT::LCEvent> readEvent(int runNumber, int eventNumber) ;
    
  private:
    /**
     *  @brief  Proxy method to open the LCIO file
     */
    void openFile() ;

  private:
    /// The LCIO file reader
    std::shared_ptr<FileReader>                    _lcReader {nullptr} ;   
    /// The run and event number map
    std::vector<int>                               _eventMap {} ; 
    /// The LCIO file name
    std::string                                    _fileName {} ;
  };
  
  typedef std::vector<OverlayFileHandler> OverlayFileHandlerList;

}

#endif
