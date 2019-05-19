#ifndef MARLINRECOMT_SIMPLECLUSTERSMEARER_h
#define MARLINRECOMT_SIMPLECLUSTERSMEARER_h 1

// -- MarlinRecoMT headers
#include <MarlinRecoMT/IFourVectorSmearer.h>

// -- std headers
#include <vector>
#include <chrono>
#include <random>

namespace marlinreco_mt {

  /** Small helper class to store calorimeter resolutions for a polar angle range
   *
   *  @author F. Gaede, DESY
   *  @version $Id: SimpleClusterSmearer.h,v 1.2 2005-10-11 12:56:28 gaede Exp $
   */
  struct ClusterResolution {

    ClusterResolution() = default ;

    ClusterResolution(float a, float b, float thMin, float thMax) :
      A(a) ,
      B(b) ,
      ThMin(thMin) ,
      ThMax(thMax) {}

    float A {0.} ;
    float B {0.};
    float ThMin {0.};
    float ThMax {0.};

  } ;

  /** Smears the four vectors of (neutral) clusters according to
   *  dE/E = A "+" B / sqrt( E/GeV ) for a given range of the polar angle.
   *  The resolution parameters A and B are given in the constructor as
   *  a vector holding quadruples of:
   *  A0, B0, th_min0, th_max0, A0, B0, th_min1, th_max1, ...<br>
   *  Momenta are assigned assumig massless particles.
   */

  class SimpleClusterSmearer : public IFourVectorSmearer {
    typedef std::vector<ClusterResolution> ResVec ;

  public:
    SimpleClusterSmearer(const std::vector<float>& resVec ) ;

    /** Smears the given four vector according to the resolution for the
     *  polar angle of the cluster. Returns a vector with all elements 0. if
     *  no resolution is defined.
     */
    TLorentzVector smearedFourVector( const TLorentzVector& v, int pdgCode ) ;

  protected:
    ResVec             _resVec {} ;
    std::mt19937       _ranEngine { unsigned(std::chrono::system_clock::now().time_since_epoch().count()) } ;
  } ;

} // end namespace

#endif // SimpleClusterSmearer_h
