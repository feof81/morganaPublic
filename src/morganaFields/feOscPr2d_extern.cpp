#include "feOscPr2d_extern.h"

// R = 0 - NW = 0 -------------------------------------------------------------
feOscPr2d_dirGen<0,0>::
feOscPr2d_dirGen()
{
  listH.resize(1);
  listH(1) = point3d(0.0, 0.0, 0.0);
  
  listY.resize(1);
  listY(1) = point3d(1.0/3.0, 1.0/3.0, 0.0);
}

const sVect<point3d> &
feOscPr2d_dirGen<0,0>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscPr2d_dirGen<0,0>::
getY() const
{
  return(listY);
}


// R = 1 - NW = 0 -------------------------------------------------------------
feOscPr2d_dirGen<1,0>::
feOscPr2d_dirGen()
{
  listH.resize(3);
  listH(1) = point3d(0.0, 0.0, 0.0);
  listH(2) = point3d(0.0, 0.0, 0.0);
  listH(3) = point3d(0.0, 0.0, 0.0);
  
  listY.resize(3);
  listY(1) = point3d(0.0, 0.0, 0.0);
  listY(2) = point3d(1.0, 0.0, 0.0);
  listY(3) = point3d(0.0, 1.0, 0.0);
}

const sVect<point3d> &
feOscPr2d_dirGen<1,0>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscPr2d_dirGen<1,0>::
getY() const
{
  return(listY);
}

// R = 0 - NW = 0 -------------------------------------------------------------
feOscPr2d_dirGen<0,1>::
feOscPr2d_dirGen()
{
  listH.resize(4);
  listH(1) = point3d( 1.0,  0.0, 0.0);
  listH(2) = point3d(-1.0,  0.0, 0.0);
  listH(3) = point3d( 0.0,  1.0, 0.0);
  listH(4) = point3d( 0.0, -1.0, 0.0);
  
  listY.resize(4);
  listY(1) = point3d(1.0/3.0, 1.0/3.0, 0.0);
  listY(2) = point3d(1.0/3.0, 1.0/3.0, 0.0);
  listY(3) = point3d(1.0/3.0, 1.0/3.0, 0.0);
  listY(4) = point3d(1.0/3.0, 1.0/3.0, 0.0);
}

const sVect<point3d> &
feOscPr2d_dirGen<0,1>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscPr2d_dirGen<0,1>::
getY() const
{
  return(listY);
}

// R = 1 - NW = 1 -------------------------------------------------------------
feOscPr2d_dirGen<1,1>::
feOscPr2d_dirGen()
{
  listH.resize(12);
  listH(1) = point3d( 1.0,  0.0, 0.0);
  listH(2) = point3d(-1.0,  0.0, 0.0);
  listH(3) = point3d( 0.0,  1.0, 0.0);
  listH(4) = point3d( 0.0, -1.0, 0.0);
  
  listH(5) = point3d( 1.0,  0.0, 0.0);
  listH(6) = point3d(-1.0,  0.0, 0.0);
  listH(7) = point3d( 0.0,  1.0, 0.0);
  listH(8) = point3d( 0.0, -1.0, 0.0);
  
  listH(9)  = point3d( 1.0,  0.0, 0.0);
  listH(10) = point3d(-1.0,  0.0, 0.0);
  listH(11) = point3d( 0.0,  1.0, 0.0);
  listH(12) = point3d( 0.0, -1.0, 0.0);
  
  listY.resize(12);
  listY(1) = point3d(0.0, 0.0, 0.0);
  listY(2) = point3d(0.0, 0.0, 0.0);
  listY(3) = point3d(0.0, 0.0, 0.0);
  listY(4) = point3d(0.0, 0.0, 0.0);
  
  listY(5) = point3d(1.0, 0.0, 0.0);
  listY(6) = point3d(1.0, 0.0, 0.0);
  listY(7) = point3d(1.0, 0.0, 0.0);
  listY(8) = point3d(1.0, 0.0, 0.0);
  
  listY(9)  = point3d(0.0, 1.0, 0.0);
  listY(10) = point3d(0.0, 1.0, 0.0);
  listY(11) = point3d(0.0, 1.0, 0.0);
  listY(12) = point3d(0.0, 1.0, 0.0);
}

const sVect<point3d> &
feOscPr2d_dirGen<1,1>::
getH() const
{
  return(listH);
}

const sVect<point3d> &
feOscPr2d_dirGen<1,1>::
getY() const
{
  return(listY);
}
