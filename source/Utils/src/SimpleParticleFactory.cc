#include <MarlinRecoMT/SimpleParticleFactory.h>

// -- std headers
#include <cstdlib>

// -- lcio headers
#include "IMPL/ReconstructedParticleImpl.h"

namespace marlinreco_mt {

  void SimpleParticleFactory::setMomentumCut( double mCut ) {
    _momentumCut = mCut ;
  }

  void SimpleParticleFactory::registerIFourVectorSmearer( IFourVectorSmearer* sm , FastMCParticleType type ) {
    _smearingVec[ type ] = sm ;
  }

  FastMCParticleType SimpleParticleFactory::getParticleType( const EVENT::MCParticle* mcp ) {
    // assumes that mcp is a stable particle !
    FastMCParticleType type( UNKNOWN )  ;
    float charge =  mcp->getCharge()  ;
    if( charge > 1e-10 || charge < -1e-10  ) {
      type = CHARGED ;
    }
    else if(  mcp->getPDG() == 22 )  { // photon
      type = PHOTON ;
    }
    else if( std::abs( mcp->getPDG() ) == 12
          || std::abs( mcp->getPDG() ) == 14
          || std::abs( mcp->getPDG() ) == 16
          || std::abs( mcp->getPDG() ) == 18 )  { // neutrinos - 18 is tau-prime
      type = NEUTRINO ;
    }
    else {  // treat everything else neutral hadron
      type = NEUTRAL_HADRON ;
    }
    return type ;
  }

  EVENT::ReconstructedParticle* SimpleParticleFactory::createReconstructedParticle( const EVENT::MCParticle* mcp ) {
    // this is where we do the fast Monte Carlo ....
    TLorentzVector mc4V( mcp->getMomentum()[0], mcp->getMomentum()[1],
			   mcp->getMomentum()[2], mcp->getEnergy() )  ;
    FastMCParticleType type = getParticleType(mcp ) ;
    IFourVectorSmearer* sm = _smearingVec[ type ] ;
    if( sm == 0 ) {
      // if we don't have a smearer registered we don't reconstruct the particle, e.g for neutrinos
      return 0 ;
    }
    TLorentzVector reco4v = sm->smearedFourVector( mc4V , mcp->getPDG() ) ;
    if( reco4v.Vect().Mag() <= _momentumCut ) {
      return nullptr ;
    }
    float p[3] ;
    float vtx[3] ;
    p[0] = reco4v.Px() ;
    p[1] = reco4v.Py() ;
    p[2] = reco4v.P() ;
    vtx[0] = mcp->getVertex()[0] ;
    vtx[1] = mcp->getVertex()[1] ;
    vtx[2] = mcp->getVertex()[2] ;
    auto rec = new IMPL::ReconstructedParticleImpl() ;
    rec->setMomentum( p ) ;
    rec->setEnergy( reco4v.E() ) ;
    rec->setMass( reco4v.M() ) ;
    rec->setCharge( mcp->getCharge() ) ;
    rec->setReferencePoint( vtx ) ;
    rec->setType( type ) ;
    return rec ;
  }


} // namespace
