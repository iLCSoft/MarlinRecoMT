#ifndef MARLINRECOMT_IRECOPARTICLEFACTORY_h
#define MARLINRECOMT_IRECOPARTICLEFACTORY_h 1

// -- lcio headers
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

namespace marlinreco_mt {


  /** Interface for a factory class that creates a ReconstructedParticle
   *  from an MCParticle
   *
   *  @author F. Gaede, DESY
   *  @version $Id: IRecoParticleFactory.h,v 1.2 2005-10-11 12:56:28 gaede Exp $
   */

  class IRecoParticleFactory {

  public:

    /** Virtual d'tor.*/
    virtual ~IRecoParticleFactory() {}

    /** The actual factory method that creates a new ReconstructedParticle
     *  for the given MCParticle. NULL if no ReconstructedParticle should be created
     *  due to detector acceptance.
     */
    virtual EVENT::ReconstructedParticle* createReconstructedParticle( const EVENT::MCParticle* mcp ) = 0 ;

  } ;


} // end namespace
#endif
