#ifndef MARLINRECOMT_SIMPLEPARTICLEFACTORY_h
#define MARLINRECOMT_SIMPLEPARTICLEFACTORY_h 1

// -- MarlinRecoMT headers
#include <MarlinRecoMT/IRecoParticleFactory.h>
#include <MarlinRecoMT/IFourVectorSmearer.h>
#include <MarlinRecoMT/FastMCParticleType.h>

namespace marlinreco_mt {

/** Implementation of IRecoParticleFactory that implements the default behaviour
 *  as described in SimpleFastMCProcessor, i.e. have polar angle ranges with different resolutions
 *  for charged tracks, photons and neutral hadrons.
 *
 *  @author F. Gaede, DESY
 *  @version $Id: SimpleParticleFactory.h,v 1.3 2007-11-23 20:09:12 gaede Exp $
 */

  class SimpleParticleFactory : public IRecoParticleFactory {
  public:
    SimpleParticleFactory() = default ;

    /** The actual factory method that creates a new ReconstructedParticle
     */
    EVENT::ReconstructedParticle* createReconstructedParticle( const EVENT::MCParticle* mcp ) ;

    /** Register a particle four vector smearer for the given type.
     */
    void registerIFourVectorSmearer( IFourVectorSmearer* sm , FastMCParticleType type ) ;

    /** Returns the type of the MCParticle.
     */
    FastMCParticleType getParticleType( const EVENT::MCParticle* mcp ) ;

    /** Set the momentum cut in GeV - no particles are produced for ower momenta. Default is 0.1 eV.
     */
    void setMomentumCut( double mCut ) ;

  protected:
    std::vector<IFourVectorSmearer*> _smearingVec {NUMBER_OF_FASTMCPARTICLETYPES, NULL} ;
    double _momentumCut {0.0000000001} ;
  };

}

#endif // SimpleParticleFactory_h
