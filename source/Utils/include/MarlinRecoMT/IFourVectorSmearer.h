#ifndef MARLINRECOMT_IFOURVECTORSMEARER_h
#define MARLINRECOMT_IFOURVECTORSMEARER_h 1

// -- ROOT headers
#include <TLorentzVector.h>

namespace marlinreco_mt {

  /** Interface for smearing of four vectors - based on TLorentzVector
   *
   *  @author F. Gaede, DESY
   *  @version $Id: IFourVectorSmearer.h,v 1.3 2006-03-30 16:12:16 gaede Exp $
   */
  class IFourVectorSmearer {
  public:

    /** Virtual d'tor.*/
    virtual ~IFourVectorSmearer() {}

    /** Smears the given four vector
     */
    virtual TLorentzVector smearedFourVector( const TLorentzVector& v, int pdgCode ) = 0 ;
  } ;


} // end namespace

#endif // IFourVectorSmearer_h
