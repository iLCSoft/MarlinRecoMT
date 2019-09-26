#include <MarlinRecoMT/VoxelTPC.h>

namespace marlinreco_mt {

  VoxelTPC::VoxelTPC(int row, int phi, int z, double pos[3], double* /*posRPhi[2]*/, double edep, double RPhiRes, double ZRes)
  {
    _row_index = row;
    _phi_index = phi;
    _z_index = z;
    _coord.SetX(pos[0]);
    _coord.SetY(pos[1]);
    _coord.SetZ(pos[2]);
    _edep = edep;
    _rPhiRes = RPhiRes;
    _zRes = ZRes;
    _isMerged = false;
    _isClusterHit = false;
  }

  VoxelTPC::VoxelTPC(int row, int phi, int z, const LCVector3D &coord, double edep, double RPhiRes, double ZRes)
  {
    _row_index = row;
    _phi_index = phi;
    _z_index = z;
    _coord=coord;
    _edep = edep;
    _rPhiRes = RPhiRes;
    _zRes = ZRes;
    _isMerged = false;
    _isClusterHit = false;
  }

  int VoxelTPC::clusterFind(std::vector <VoxelTPC*>* hitList){
    
    if(!this->IsClusterHit()){
      hitList->push_back(this);
      this->setIsClusterHit();
      for(int i=0; i<this->getNumberOfAdjacent();++i){
        getAdjacent(i)->clusterFind(hitList);
      }
    }

    return hitList->size();
  }

}

