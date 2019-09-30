#ifndef _voxeltpc_included_
#define _voxeltpc_included_ 1

// -- std headers
#include <vector>

// -- marlinreco mt headers
#include <MarlinRecoMT/LCGeometryTypes.h>

namespace marlinreco_mt {

  class VoxelTPC {

   public:
    VoxelTPC() = default ;
    ~VoxelTPC() = default ;
    
    // the intialation in the constructor here would be preferable though I don't know how to intialise
    // the array xyz[3] here with pos[3], for the mean time the constructor will be put in the .cc file
    //  VoxelTPC(int row, int phi, int z, double pos[3]) : row_index(row), phi_index(phi), z_index(z){}
    VoxelTPC(int row, int phi, int z, double pos[3], double posRPhi[2], double edep, double rPhiRes, double zRes);
    VoxelTPC(int row, int phi, int z, const LCVector3D &coord, double edep, double rPhiRes, double zRes);

    void setAdjacent(VoxelTPC * p_voxel) { _adjacent_voxels.push_back(p_voxel);}
    void setIsClusterHit() { _isClusterHit = true;};
    void setIsMerged() { _isMerged = true;};
    bool IsClusterHit() const { return _isClusterHit;}
    bool IsMerged() const { return _isMerged;}
    int clusterFind(std::vector <VoxelTPC*>* hitList);
    

    int getRowIndex() const {return _row_index;}
    int getPhiIndex() const {return _phi_index;}
    int getZIndex() const {return _z_index;}
    VoxelTPC * getFirstAdjacent() const {return *(_adjacent_voxels.begin());}
    VoxelTPC * getAdjacent(int i) const {return _adjacent_voxels[i];}
    int getNumberOfAdjacent() const {return _adjacent_voxels.size();}
    double getX() const {return _coord.x();}
    double getY() const {return _coord.y();}
    double getZ() const {return _coord.z();}
    double getR() const {return _coord.rho();}
    double getPhi() const {return _coord.phi();}
    double getEDep() const {return _edep;}
    double getRPhiRes() const {return _rPhiRes;}
    double getZRes() const {return _zRes;}
    const LCVector3D &getHep3Vector() const {return _coord;}
    
    static comparePhi( const VoxelTPC *const a, const VoxelTPC *const b )  { 
      return ( a->getPhiIndex() < b->getPhiIndex() ) ; 
    } 
    static bool compareZ( const VoxelTPC *const a, const VoxelTPC *const b ) { 
      return ( a->getZIndex() < b->getZIndex() ) ; 
    } 

   private:
    int _row_index{}; 
    int _phi_index{};
    int _z_index{};
    std::vector <VoxelTPC *> _adjacent_voxels{};
    LCVector3D _coord{};
    double _edep{};
    double _rPhiRes{};
    double _zRes{};
    bool _isMerged{};
    bool _isClusterHit{};
  };
  
}
  
#endif
