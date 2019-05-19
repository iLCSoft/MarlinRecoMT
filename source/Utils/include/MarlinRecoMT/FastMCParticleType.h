#ifndef MARLINRECOMT_FASTMCPARTICLETYPE_h
#define MARLINRECOMT_FASTMCPARTICLETYPE_h 1

namespace marlinreco_mt {

  /** Enumeration that defines integer constants for various particle
   *  types used in the fast Monte Carlo.
   *
   *  @author F. Gaede, DESY
   *  @version $Id: FastMCParticleType.h,v 1.2 2005-10-11 12:56:28 gaede Exp $
   */
  enum FastMCParticleType {
    UNKNOWN = 0 ,
    CHARGED,
    PHOTON,
    NEUTRAL_HADRON,
    NEUTRINO,
    NUMBER_OF_FASTMCPARTICLETYPES
  };

} // end namespace

#endif
