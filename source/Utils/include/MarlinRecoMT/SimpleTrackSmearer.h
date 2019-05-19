#ifndef MARLINRECOMT_SIMPLETRACKSMEARER_h
#define MARLINRECOMT_SIMPLETRACKSMEARER_h 1

// -- MarlinRecoMT headers
#include <MarlinRecoMT/IFourVectorSmearer.h>

// -- std headers
#include <vector>
#include <chrono>
#include <random>

namespace marlinreco_mt {

  /** Small helper class to store tracker  resolutions for a polar angle range
   *
   *  @author F. Gaede, DESY
   *  @version $Id: SimpleTrackSmearer.h,v 1.2 2005-10-11 12:56:28 gaede Exp $
   */
  struct TrackResolution {

    TrackResolution() = default ;

    TrackResolution(float dPP, float thMin, float thMax) :
      DPP(dPP) ,
      ThMin(thMin) ,
      ThMax(thMax) {}

    float DPP {0.} ;
    float ThMin {0.} ;
    float ThMax {0.} ;
  } ;


  /** Smears the four vectors of charged tracks according to
   *  dP/P = r for a given range of the polar angle.
   *  The resolutions r are given in the constructor as
   *  a vector holding triplets of:
   *  r0, th_min0, th_max0, r1, th_min1, th_max1, ...<br>
   *  Perfect electron and muon ID is assumed to set the
   *  particle energy ( mass ) - all other particles are
   *  assumed to be pions.
   */
  class SimpleTrackSmearer : public IFourVectorSmearer {
    typedef std::vector<TrackResolution> ResVec ;
    static constexpr float ELECTRON_MASS = 0.0005109989 ;
    static constexpr float MUON_MASS = 0.10565836 ;
    static constexpr float PION_MASS = 0.139570 ;

  public:
    SimpleTrackSmearer(const std::vector<float>& resVec ) ;

    /** Smears the given four vector according to the resolution for the
     *  polar angle of the track. Returns a vector with all elements 0. if
     *  no resolution is defined.
     */
    TLorentzVector smearedFourVector( const TLorentzVector& v, int pdgCode ) ;

  protected:
    ResVec             _resVec {} ;
    std::mt19937       _ranEngine { unsigned(std::chrono::system_clock::now().time_since_epoch().count()) } ;
  } ;

} // end namespace

#endif // SimpleTrackSmearer_h
