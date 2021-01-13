#include "feOscPr3d_extern.h"


// R = 0 - NW = 1 -------------------------------------------------------------
feOscPr3d_dirGen<0,1>::
feOscPr3d_dirGen()
{
  listH.resize(6);
  listH(1) = point3d( 1.0,  0.0,  0.0);
  listH(2) = point3d(-1.0,  0.0,  0.0);
  listH(3) = point3d( 0.0,  1.0,  0.0);
  listH(4) = point3d( 0.0, -1.0,  0.0);
  listH(5) = point3d( 0.0,  0.0,  1.0);
  listH(6) = point3d( 0.0,  0.0, -1.0 );
  
  listY.resize(6);
  listY(1) = point3d(0.25, 0.25, 0.25);
  listY(2) = point3d(0.25, 0.25, 0.25);
  listY(3) = point3d(0.25, 0.25, 0.25);
  listY(4) = point3d(0.25, 0.25, 0.25);
  listY(5) = point3d(0.25, 0.25, 0.25);
  listY(6) = point3d(0.25, 0.25, 0.25);
}

const sVect<point3d> &
feOscPr3d_dirGen<0,1>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscPr3d_dirGen<0,1>::
getY() const
{
  return(listY);
}


// R = 1 - NW = 1 -------------------------------------------------------------
feOscPr3d_dirGen<1,1>::
feOscPr3d_dirGen()
{
  listH.resize(24);
  listH(1)  = point3d( 1.0,  0.0,  0.0);
  listH(2)  = point3d(-1.0,  0.0,  0.0);
  listH(3)  = point3d( 0.0,  1.0,  0.0);
  listH(4)  = point3d( 0.0, -1.0,  0.0);
  listH(5)  = point3d( 0.0,  0.0,  1.0);
  listH(6)  = point3d( 0.0,  0.0, -1.0);
  
  listH(7)  = point3d( 1.0,  0.0,  0.0);
  listH(8)  = point3d(-1.0,  0.0,  0.0);
  listH(9)  = point3d( 0.0,  1.0,  0.0);
  listH(10) = point3d( 0.0, -1.0,  0.0);
  listH(11) = point3d( 0.0,  0.0,  1.0);
  listH(12) = point3d( 0.0,  0.0, -1.0);
  
  listH(13) = point3d( 1.0,  0.0,  0.0);
  listH(14) = point3d(-1.0,  0.0,  0.0);
  listH(15) = point3d( 0.0,  1.0,  0.0);
  listH(16) = point3d( 0.0, -1.0,  0.0);
  listH(17) = point3d( 0.0,  0.0,  1.0);
  listH(18) = point3d( 0.0,  0.0, -1.0);
  
  listH(19) = point3d( 1.0,  0.0,  0.0);
  listH(20) = point3d(-1.0,  0.0,  0.0);
  listH(21) = point3d( 0.0,  1.0,  0.0);
  listH(22) = point3d( 0.0, -1.0,  0.0);
  listH(23) = point3d( 0.0,  0.0,  1.0);
  listH(24) = point3d( 0.0,  0.0, -1.0);
  
  listY.resize(24);
  listY(1)  = point3d(0.0, 0.0, 0.0);
  listY(2)  = point3d(0.0, 0.0, 0.0);
  listY(3)  = point3d(0.0, 0.0, 0.0);
  listY(4)  = point3d(0.0, 0.0, 0.0);
  listY(5)  = point3d(0.0, 0.0, 0.0);
  listY(6)  = point3d(0.0, 0.0, 0.0);
  
  listY(7)  = point3d(1.0, 0.0, 0.0);
  listY(8)  = point3d(1.0, 0.0, 0.0);
  listY(9)  = point3d(1.0, 0.0, 0.0);
  listY(10) = point3d(1.0, 0.0, 0.0);
  listY(11) = point3d(1.0, 0.0, 0.0);
  listY(12) = point3d(1.0, 0.0, 0.0);
  
  listY(13) = point3d(0.0, 1.0, 0.0);
  listY(14) = point3d(0.0, 1.0, 0.0);
  listY(15) = point3d(0.0, 1.0, 0.0);
  listY(16) = point3d(0.0, 1.0, 0.0);
  listY(17) = point3d(0.0, 1.0, 0.0);
  listY(18) = point3d(0.0, 1.0, 0.0);
  
  listY(19) = point3d(0.0, 0.0, 1.0);
  listY(20) = point3d(0.0, 0.0, 1.0);
  listY(21) = point3d(0.0, 0.0, 1.0);
  listY(22) = point3d(0.0, 0.0, 1.0);
  listY(23) = point3d(0.0, 0.0, 1.0);
  listY(24) = point3d(0.0, 0.0, 1.0);
}

const sVect<point3d> &
feOscPr3d_dirGen<1,1>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscPr3d_dirGen<1,1>::
getY() const
{
  return(listY);
}
