#ifndef LCGeometryTypes_H
#define LCGeometryTypes_H 1

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/SMatrix.h>

namespace marlinreco_mt {
  
  /** 
   *  @file Definition of geometry types used in ILC software - currently use ROOT.
   *  @author R.Ete, DESY
   */
  // 2D, 3D and Lorentz vectors
  using LCVector2D = ROOT::Math::XYVector ;
  using LCVector3D = ROOT::Math::XYZVector ;
  using LCLorentzVector = ROOT::Math::XYZTVector ;

  // Error matrix. To be used as e.g 'LCErrorMatrix<4>'
  template <unsigned int N>
  using LCErrorMatrix = ROOT::Math::SMatrix<double, N, N, ROOT::Math::MatRepSym<double, N>> ;

  namespace vector {
    
    /**
     *  @brief  Helper function to get an orthogonal vector 
     *          Ported from CLHEP. Lacking from ROOT interface
     */
    inline LCVector3D orthogonal( const LCVector3D &vec ) {
      double xx = std::fabs(vec.x());
      double yy = std::fabs(vec.y());
      double zz = std::fabs(vec.z());
      if (xx < yy) {
        return xx < zz ? LCVector3D(0,vec.z(),-vec.y()) : LCVector3D(vec.y(),-vec.x(),0);
      }
      else{
        return yy < zz ? LCVector3D(-vec.z(),0,vec.x()) : LCVector3D(vec.y(),-vec.x(),0);
      }
    }
    
  }
    
}

#endif /* ifndef LCGeometryTypes_H */
