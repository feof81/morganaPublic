#include "feOscCr3d_extern.h"


// R = 0 - NW = 1 -------------------------------------------------------------
feOscCr3d_dirGen<1>::
feOscCr3d_dirGen()
{
  listH.resize(4);
  listH(1) = point3d( 1.0,  1.0,  1.0) / sqrt(3);
  listH(2) = point3d(-1.0,  0.0,  0.0);
  listH(3) = point3d( 0.0, -1.0,  0.0);
  listH(4) = point3d( 0.0,  0.0, -1.0);
  
  listY.resize(4);
  listY(1) = point3d(1.0,1.0,1.0) / 3.0;
  listY(2) = point3d(0.0,1.0,1.0) / 3.0;
  listY(3) = point3d(1.0,0.0,1.0) / 3.0;
  listY(4) = point3d(1.0,1.0,0.0) / 3.0;
}

const sVect<point3d> &
feOscCr3d_dirGen<1>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscCr3d_dirGen<1>::
getY() const
{
  return(listY);
}
