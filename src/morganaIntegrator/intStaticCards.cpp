/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#include "intStaticCards.h"


/*! -------------------------------------------------------------------------------------------- */
/*!                                       INTEGRATION 1D                                         */
/*! -------------------------------------------------------------------------------------------- */


//! LinearLine, STANDARD, Precision 1----------------------------------------------------------------
point3d intStaticCard<linearLine,STANDARD,1>::getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {point3d(0.5, 0.0, 0.0)};
  
  return(Yn[i-1]);
}

Real intStaticCard<linearLine,STANDARD,1>::getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[1] = {1.0};
  
  return(Wn[i-1]);
}


//! LinearLine, STANDARD, Precision 2---------------------------------------------------------------
point3d intStaticCard<linearLine,STANDARD,2>::getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q2ptx1 = ( 1 - sqrt( 1. / 3. ) ) / 2., q2ptx2 = ( 1 + sqrt( 1. / 3. ) ) / 2.;
  static const point3d Yn[2] = {
    point3d(q2ptx1, 0.0, 0.0),
    point3d(q2ptx2, 0.0, 0.0) };
  
  return(Yn[i-1]);
}

Real intStaticCard<linearLine,STANDARD,2>::getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static Real q2ptw1 = 0.5, q2ptw2 = 0.5;
  static const Real Wn[2] = {q2ptw1, q2ptw2};
  
  return(Wn[i-1]);
}


//! LinearLine, STANDARD, Precision 3---------------------------------------------------------------
point3d intStaticCard<linearLine,STANDARD,3>::getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q3ptx1 = 0.5, q3ptx2 = ( 1 - sqrt( 3. / 5. ) ) / 2., q3ptx3 = ( 1 + sqrt( 3. / 5. ) ) / 2.;
  static const point3d Yn[3] = {
     point3d(q3ptx1, 0.0, 0.0),
     point3d(q3ptx2, 0.0, 0.0),
     point3d(q3ptx3, 0.0, 0.0) };
     
  return(Yn[i-1]);
}

Real intStaticCard<linearLine,STANDARD,3>::getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q3ptw1 = 8. / 18., q3ptw2 = 5. / 18., q3ptw3 = 5. / 18.;
  static const Real Wn[3] = {q3ptw1, q3ptw2, q3ptw3};

  return(Wn[i-1]);
}


//! LinearLine, HYBRID, Precision 1----------------------------------------------------------------
point3d intStaticCard<linearLine,HYBRID,1>::getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {point3d(0.5, 0.0, 0.0)};
  
  return(Yn[i-1]);
}

Real intStaticCard<linearLine,HYBRID,1>::getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[1] = {1.0};
  
  return(Wn[i-1]);
}


//! LinearLine, HYBRID, Precision 2---------------------------------------------------------------
point3d intStaticCard<linearLine,HYBRID,2>::getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q2ptx1 = ( 1 - sqrt( 1. / 3. ) ) / 2., q2ptx2 = ( 1 + sqrt( 1. / 3. ) ) / 2.;
  static const point3d Yn[2] = {
    point3d(q2ptx1, 0.0, 0.0),
    point3d(q2ptx2, 0.0, 0.0) };
  
  return(Yn[i-1]);
}

Real intStaticCard<linearLine,HYBRID,2>::getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static Real q2ptw1 = 0.5, q2ptw2 = 0.5;
  static const Real Wn[2] = {q2ptw1, q2ptw2};
  
  return(Wn[i-1]);
}


//! LinearLine, HYBRID, Precision 3---------------------------------------------------------------
point3d intStaticCard<linearLine,HYBRID,3>::getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q3ptx1 = 0.5, q3ptx2 = ( 1 - sqrt( 3. / 5. ) ) / 2., q3ptx3 = ( 1 + sqrt( 3. / 5. ) ) / 2.;
  static const point3d Yn[3] = {
     point3d(q3ptx1, 0.0, 0.0),
     point3d(q3ptx2, 0.0, 0.0),
     point3d(q3ptx3, 0.0, 0.0) };
     
  return(Yn[i-1]);
}

Real intStaticCard<linearLine,HYBRID,3>::getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q3ptw1 = 8. / 18., q3ptw2 = 5. / 18., q3ptw3 = 5. / 18.;
  static const Real Wn[3] = {q3ptw1, q3ptw2, q3ptw3};

  return(Wn[i-1]);
}



/*! -------------------------------------------------------------------------------------------- */
/*!                                       INTEGRATION 2D                                         */
/*! -------------------------------------------------------------------------------------------- */

//! LinearTriangle, STANDARD, Precision 1------------------------------------------------------------
point3d
intStaticCard<linearTriangle,STANDARD,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {point3d(1.0/3.0, 1.0/3.0, 0.0)};
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,STANDARD,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[1] = {0.5};
  return(Wn[i-1]);
}


//! LinearTriangle, STANDARD, Precision 2------------------------------------------------------------
point3d
intStaticCard<linearTriangle,STANDARD,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[3] = {
    point3d(0.5, 0.0, 0.0),
    point3d(0.0, 0.5, 0.0),
    point3d(0.5, 0.5, 0.0) };
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,STANDARD,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[3] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
  return(Wn[i-1]);
}


//! LinearTriangle, STANDARD, Precision 3------------------------------------------------------------
point3d
intStaticCard<linearTriangle,STANDARD,3>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t4pt_xb1 =  3.0 / 5.0;
  static const Real t4pt_xb2 =  1.0 / 5.0;
  static const Real t4pt_a   =  1.0 / 3.0;
  
  static const point3d Yn[4] = {
    point3d(t4pt_xb1, t4pt_xb2, 0.0),
    point3d(t4pt_xb2, t4pt_xb1, 0.0),
    point3d(t4pt_xb2, t4pt_xb2, 0.0),
    point3d(t4pt_a, t4pt_a, 0.0) };
     
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,STANDARD,3>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t4pt_w1  = 25.0 / 96.0;
  static const Real t4pt_w2  = -9.0 / 32.0;
  static const Real Wn[4] = {t4pt_w1, t4pt_w1, t4pt_w1, t4pt_w2};
  
  return(Wn[i-1]);
}


//! LinearTriangle, P0DUAL, Precision 4------------------------------------------------------------
point3d
intStaticCard<linearTriangle,STANDARD,4>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t6pt_x1 = 0.091576213509770743;
  static const Real t6pt_x2 = 0.44594849091596488;
  static const point3d Yn[6] = {
    point3d(    t6pt_x1,     t6pt_x1, 0.0),
    point3d(    t6pt_x1, 1-2*t6pt_x1, 0.0),
    point3d(1-2*t6pt_x1,     t6pt_x1, 0.0),
    point3d(    t6pt_x2,     t6pt_x2, 0.0),
    point3d(    t6pt_x2, 1-2*t6pt_x2, 0.0),
    point3d(1-2*t6pt_x2,     t6pt_x2, 0.0)
  };
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,STANDARD,4>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t6pt_w1 = 0.054975871827660933;
  static const Real t6pt_w2 = 0.11169079483900573;
  static const Real Wn[6] = {t6pt_w1, t6pt_w1, t6pt_w1, t6pt_w2, t6pt_w2, t6pt_w2};
  
  return(Wn[i-1]);
}


//! LinearTriangle, P0DUAL, Precision 1------------------------------------------------------------
point3d
intStaticCard<linearTriangle,MINIELEMENT,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d P1(0.0,0.0,0.0);
  static const point3d P2(1.0,0.0,0.0);
  static const point3d P3(0.0,1.0,0.0);
  static const point3d PM = (P1+P2+P3)/3.0;
  
  static const point3d Yn[3] = {
    (P1+P2+PM) / 3.0,
    (P2+P3+PM) / 3.0,
    (P3+P1+PM) / 3.0 };
    
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,MINIELEMENT,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[3] = {1.0/ 6.0, 1.0/ 6.0, 1.0/6.0};
  
  return(Wn[i-1]);
}


//! LinearTriangle, P0DUAL, Precision 1------------------------------------------------------------
point3d
intStaticCard<linearTriangle,P0DUAL,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d P1(0.0,0.0,0.0);
  static const point3d P2(1.0,0.0,0.0);
  static const point3d P3(0.0,1.0,0.0);
  static const point3d PM = (P1+P2+P3)/3.0;
  
  static const point3d Yn[3] = {
    (P1+PM) / 3.0,
    (P2+PM) / 3.0,
    (P3+PM) / 3.0 };
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearTriangle,P0DUAL,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[3] = {1.0/ 6.0, 1.0/ 6.0, 1.0/6.0};
  
  return(Wn[i-1]);
}


//! LinearTriangle, HYBRID, Precision 1------------------------------------------------------------
point3d
intStaticCard<linearTriangle,HYBRID,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {point3d(1.0/3.0, 1.0/3.0, 0.0)};
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,HYBRID,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[1] = {0.5};
  return(Wn[i-1]);
}


//! LinearTriangle, HYBRID, Precision 2------------------------------------------------------------
point3d
intStaticCard<linearTriangle,HYBRID,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t4pt_xb1 =  3.0 / 5.0;
  static const Real t4pt_xb2 =  1.0 / 5.0;
  static const Real t4pt_a   =  1.0 / 3.0;
  
  static const point3d Yn[4] = {
    point3d(t4pt_xb1, t4pt_xb2, 0.0),
    point3d(t4pt_xb2, t4pt_xb1, 0.0),
    point3d(t4pt_xb2, t4pt_xb2, 0.0),
    point3d(t4pt_a, t4pt_a, 0.0) };
     
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,HYBRID,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t4pt_w1  = 25.0 / 96.0;
  static const Real t4pt_w2  = -9.0 / 32.0;
  static const Real Wn[4] = {t4pt_w1, t4pt_w1, t4pt_w1, t4pt_w2};
  
  return(Wn[i-1]);
}


//! LinearTriangle, HYBRID, Precision 3------------------------------------------------------------
point3d
intStaticCard<linearTriangle,HYBRID,3>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t6pt_x1 = 0.091576213509770743;
  static const Real t6pt_x2 = 0.44594849091596488;
  static const point3d Yn[6] = {
    point3d(    t6pt_x1,     t6pt_x1, 0.0),
    point3d(    t6pt_x1, 1-2*t6pt_x1, 0.0),
    point3d(1-2*t6pt_x1,     t6pt_x1, 0.0),
    point3d(    t6pt_x2,     t6pt_x2, 0.0),
    point3d(    t6pt_x2, 1-2*t6pt_x2, 0.0),
    point3d(1-2*t6pt_x2,     t6pt_x2, 0.0)
  };
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTriangle,HYBRID,3>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t6pt_w1 = 0.054975871827660933;
  static const Real t6pt_w2 = 0.11169079483900573;
  static const Real Wn[6] = {t6pt_w1, t6pt_w1, t6pt_w1, t6pt_w2, t6pt_w2, t6pt_w2};
  
  return(Wn[i-1]);
}


//! LinearQuad, STANDARD, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearQuad,STANDARD,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = 0.0;
  static const Real p2 = 0.5;
  static const Real p3 = 1.0;
  
  static const point3d Yn[9] = {
  point3d(p1,p1,0.0), point3d(p1,p2,0.0), point3d(p1,p3,0.0),
  point3d(p2,p1,0.0), point3d(p2,p2,0.0), point3d(p2,p3,0.0),
  point3d(p3,p1,0.0), point3d(p3,p2,0.0), point3d(p3,p3,0.0)};
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearQuad,STANDARD,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real w1 = 1.0/6.0;
  static const Real w2 = 2.0/3.0;
  static const Real w3 = 1.0/6.0;
  
  static const Real Wn[9] = {
  w1 * w1,  w1 * w2,  w1 * w3,
  w2 * w1,  w2 * w2,  w2 * w3,
  w3 * w1,  w3 * w2,  w3 * w3};
  
  return(Wn[i-1]);
}


//! LinearQuad, STANDARD, Precision 2---------------------------------------------------------------
point3d
intStaticCard<linearQuad,STANDARD,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = 0.0;
  static const Real p2 = 1.0/2.0 - sqrt(1.0/5.0);
  static const Real p3 = 1.0/2.0 + sqrt(1.0/5.0);
  static const Real p4 = 1.0;
  
  static const point3d Yn[16] = {
  point3d(p1,p1,0.0), point3d(p1,p2,0.0), point3d(p1,p3,0.0), point3d(p1,p4,0.0),
  point3d(p2,p1,0.0), point3d(p2,p2,0.0), point3d(p2,p3,0.0), point3d(p2,p4,0.0),
  point3d(p3,p1,0.0), point3d(p3,p2,0.0), point3d(p3,p3,0.0), point3d(p3,p4,0.0),
  point3d(p4,p1,0.0), point3d(p4,p2,0.0), point3d(p4,p3,0.0), point3d(p4,p4,0.0)};
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearQuad,STANDARD,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real w1 = 1.0/12.0;
  static const Real w2 = 5.0/12.0;
  static const Real w3 = 5.0/12.0;
  static const Real w4 = 1.0/12.0;
  
  static const Real Wn[16] = {
  w1 * w1,  w1 * w2,  w1 * w3,  w1 * w4,
  w2 * w1,  w2 * w2,  w2 * w3,  w2 * w4,
  w3 * w1,  w3 * w2,  w3 * w3,  w3 * w4,
  w4 * w1,  w4 * w2,  w4 * w3,  w4 * w4};
  
  return(Wn[i-1]);
}


//! LinearQuad, STANDARD, Precision 3---------------------------------------------------------------
point3d
intStaticCard<linearQuad,STANDARD,3>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = 0.0;
  static const Real p2 = 1.0/2.0 - sqrt(3.0/7.0);
  static const Real p3 = 0.5;
  static const Real p4 = 1.0/2.0 + sqrt(3.0/7.0);
  static const Real p5 = 1.0;
  
  static const point3d Yn[25] = {
  point3d(p1,p1,0.0), point3d(p1,p2,0.0), point3d(p1,p3,0.0), point3d(p1,p4,0.0), point3d(p1,p5,0.0),
  point3d(p2,p1,0.0), point3d(p2,p2,0.0), point3d(p2,p3,0.0), point3d(p2,p4,0.0), point3d(p2,p5,0.0),
  point3d(p3,p1,0.0), point3d(p3,p2,0.0), point3d(p3,p3,0.0), point3d(p3,p4,0.0), point3d(p3,p5,0.0),
  point3d(p4,p1,0.0), point3d(p4,p2,0.0), point3d(p4,p3,0.0), point3d(p4,p4,0.0), point3d(p4,p5,0.0),
  point3d(p5,p1,0.0), point3d(p5,p2,0.0), point3d(p5,p3,0.0), point3d(p5,p4,0.0), point3d(p5,p5,0.0),
  };
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearQuad,STANDARD,3>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real w1 =  1.0 /  20.0;
  static const Real w2 = 49.0 / 180.0;
  static const Real w3 = 32.0 /  90.0;
  static const Real w4 = 49.0 / 180.0;
  static const Real w5 =  1.0 /  20.0;
  
  static const Real Wn[25] = {
  w1 * w1,  w1 * w2,  w1 * w3,  w1 * w4,  w1 * w5,
  w2 * w1,  w2 * w2,  w2 * w3,  w2 * w4,  w2 * w5,
  w3 * w1,  w3 * w2,  w3 * w3,  w3 * w4,  w3 * w5,
  w4 * w1,  w4 * w2,  w4 * w3,  w4 * w4,  w4 * w5,
  w5 * w1,  w5 * w2,  w5 * w3,  w5 * w4,  w5 * w5, };
  
  return(Wn[i-1]);
}


//! LinearQuad, HYBRID, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearQuad,HYBRID,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {
  point3d(0.5, 0.5, 0.0)};
    
  return(Yn[i-1]);
}

Real
intStaticCard<linearQuad,HYBRID,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);

  static const Real Wn[1] = { 1 };
  
  return(Wn[i-1]);
}


//! LinearQuad, HYBRID, Precision 2---------------------------------------------------------------
point3d
intStaticCard<linearQuad,HYBRID,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = ( (-1.0 / sqrt(3.0)) + 1.0  ) / 2.0;
  static const Real p2 = ( ( 1.0 / sqrt(3.0)) + 1.0  ) / 2.0;
  
  static const point3d Yn[4] = {
  point3d(p1, p1, 0.0), point3d(p1, p2, 0.0),
  point3d(p2, p1, 0.0), point3d(p2, p2, 0.0)
  };
    
  return(Yn[i-1]);
}

Real
intStaticCard<linearQuad,HYBRID,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);

  static const Real Wn[4] = { 0.25, 0.25, 0.25, 0.25 };
  
  return(Wn[i-1]);
}


//! LinearQuad, HYBRID, Precision 3---------------------------------------------------------------
point3d
intStaticCard<linearQuad,HYBRID,3>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = ( - sqrt(3.0 / 5.0) + 1.0  ) / 2.0;
  static const Real p2 = 1.0 / 2.0;
  static const Real p3 = (   sqrt(3.0 / 5.0) + 1.0  ) / 2.0;
  
  static const point3d Yn[9] = {
  point3d(p1, p1, 0.0), point3d(p1, p2, 0.0), point3d(p1, p1, 0.0),
  point3d(p2, p1, 0.0), point3d(p2, p2, 0.0), point3d(p2, p1, 0.0),
  point3d(p2, p1, 0.0), point3d(p3, p2, 0.0), point3d(p3, p1, 0.0)
  };
    
  return(Yn[i-1]);
}

Real
intStaticCard<linearQuad,HYBRID,3>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);

  static const Real w1 = 5.0 / 18.0;
  static const Real w2 = 8.0 / 18.0;
  static const Real w3 = 5.0 / 18.0;
  
  static const Real Wn[9] = {
    w1 * w1,  w1 * w2,  w1 * w3,
    w2 * w1,  w2 * w2,  w2 * w3,
    w3 * w1,  w3 * w2,  w3 * w3
  };
  
  return(Wn[i-1]);
}



/*! -------------------------------------------------------------------------------------------- */
/*!                                       INTEGRATION 3D                                         */
/*! -------------------------------------------------------------------------------------------- */

//! LinearTetra, STANDARD, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearTetra,STANDARD,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {point3d(1.0/4.0, 1.0/4.0, 1.0/4.0)};
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTetra,STANDARD,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[1] = {1.0/ 6.0};
  
  return(Wn[i-1]);
}


//! LinearTetra, STANDARD, Precision 2---------------------------------------------------------------
point3d
intStaticCard<linearTetra,STANDARD,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real tet4ptx1 = ( 5.0 -   sqrt( 5.0 ) ) / 20.0;
  static const Real tet4ptx2 = ( 5.0 + 3*sqrt( 5.0 ) ) / 20.0;
  
  static const point3d Yn[4] = {
    point3d(tet4ptx1, tet4ptx1, tet4ptx1),
    point3d(tet4ptx1, tet4ptx1, tet4ptx2),
    point3d(tet4ptx1, tet4ptx2, tet4ptx1),
    point3d(tet4ptx2, tet4ptx1, tet4ptx1)};
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTetra,STANDARD,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[4] = {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0};
  
  return(Wn[i-1]);
}


//! LinearTetra, STANDARD, Precision 3---------------------------------------------------------------
point3d
intStaticCard<linearTetra,STANDARD,3>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real r5      = 0.25; // s
  static const Real s5[ 4 ] = {0.09197107805272303, 0.3197936278296299	}; // (7 \mp \sqrt(15))/34
  static const Real t5[ 4 ] = {0.7240867658418310, 0.04061911651111023 }; // (13 \pm 3*sqrt(15))/34
  static const Real u5      = 0.05635083268962915; // (10-2*sqrt(15))/40
  static const Real v5      = 0.4436491673103708; // (10+2*sqrt(15))/40
  
  static const point3d Yn[15] = {
    point3d( r5, r5, r5),
    point3d( s5[ 0 ], s5[ 0 ], s5[ 0 ] ),
    point3d( t5[ 0 ], s5[ 0 ], s5[ 0 ] ),
    point3d( s5[ 0 ], t5[ 0 ], s5[ 0 ] ),
    point3d( s5[ 0 ], s5[ 0 ], t5[ 0 ] ),
    point3d( s5[ 1 ], s5[ 1 ], s5[ 1 ] ),
    point3d( t5[ 1 ], s5[ 1 ], s5[ 1 ] ),
    point3d( s5[ 1 ], t5[ 1 ], s5[ 1 ] ),
    point3d( s5[ 1 ], s5[ 1 ], t5[ 1 ] ),
    point3d( u5, u5, v5 ),
    point3d( u5, v5, u5 ),
    point3d( v5, u5, u5 ),
    point3d( v5, v5, u5 ),
    point3d( v5, u5, v5 ),
    point3d( u5, v5, v5 ) };

    return(Yn[i-1]);
}

Real
intStaticCard<linearTetra,STANDARD,3>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real A5      = 0.01975308641975309; // 16/135*1/6
  static const Real B5[ 2 ] = {0.01198951396316977, 0.01151136787104540}; // 1/6*(2665 \pm 14*sqrt(15))/37800
  static const Real C5      = 0.008818342151675485; // 20/378*1/6
  
  static const Real Wn[15] = {
    A5,
    B5[ 0 ],
    B5[ 0 ],
    B5[ 0 ],
    B5[ 0 ],
    B5[ 1 ],
    B5[ 1 ],
    B5[ 1 ],
    B5[ 1 ],
    C5,
    C5,
    C5,
    C5,
    C5,
    C5 };
    
  return(Wn[i-1]); 
}


//! LinearTetra, STANDARD, Precision 4---------------------------------------------------------------
point3d
intStaticCard<linearTetra,STANDARD,4>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real t[ 4 ] = {0.0485005494, 0.2386007376, 0.5170472951, 0.7958514179 };
  static const Real s[ 4 ] = {0.0571041961, 0.2768430136, 0.5835904324, 0.8602401357 };
  static const Real r[ 4 ] = {0.0694318422, 0.3300094782, 0.6699905218, 0.9305681558 };
  
  static const point3d Yn[64] = {
  point3d( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] )), 
  point3d( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] )),
  point3d( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] )),      
  point3d( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] )),
  point3d( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] )),
  point3d( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] )),
  point3d( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] )) };
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTetra,STANDARD,4>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real A[ 4 ] = {0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226 };
  static const Real B[ 4 ] = {0.1355069134, 0.2034645680, 0.1298475476, 0.0311809709 };
  static const Real C[ 4 ] = {0.1108884156, 0.1434587898, 0.0686338872, 0.0103522407 };
  
  static const Real Wn[64] = {
  A[ 0 ] * B[ 0 ] * C[ 0 ],
  A[ 0 ] * B[ 0 ] * C[ 1 ],
  A[ 0 ] * B[ 0 ] * C[ 2 ],
  A[ 0 ] * B[ 0 ] * C[ 3 ],
  A[ 0 ] * B[ 1 ] * C[ 0 ],
  A[ 0 ] * B[ 1 ] * C[ 1 ],
  A[ 0 ] * B[ 1 ] * C[ 2 ],
  A[ 0 ] * B[ 1 ] * C[ 3 ],
  A[ 0 ] * B[ 2 ] * C[ 0 ],
  A[ 0 ] * B[ 2 ] * C[ 1 ],
  A[ 0 ] * B[ 2 ] * C[ 2 ],
  A[ 0 ] * B[ 2 ] * C[ 3 ],
  A[ 0 ] * B[ 3 ] * C[ 0 ],
  A[ 0 ] * B[ 3 ] * C[ 1 ],
  A[ 0 ] * B[ 3 ] * C[ 2 ],
  A[ 0 ] * B[ 3 ] * C[ 3 ],
  A[ 1 ] * B[ 0 ] * C[ 0 ],
  A[ 1 ] * B[ 0 ] * C[ 1 ],
  A[ 1 ] * B[ 0 ] * C[ 2 ],
  A[ 1 ] * B[ 0 ] * C[ 3 ],
  A[ 1 ] * B[ 1 ] * C[ 0 ],
  A[ 1 ] * B[ 1 ] * C[ 1 ],
  A[ 1 ] * B[ 1 ] * C[ 2 ],
  A[ 1 ] * B[ 1 ] * C[ 3 ],
  A[ 1 ] * B[ 2 ] * C[ 0 ],
  A[ 1 ] * B[ 2 ] * C[ 1 ],
  A[ 1 ] * B[ 2 ] * C[ 2 ],
  A[ 1 ] * B[ 2 ] * C[ 3 ],
  A[ 1 ] * B[ 3 ] * C[ 0 ],
  A[ 1 ] * B[ 3 ] * C[ 1 ],
  A[ 1 ] * B[ 3 ] * C[ 2 ],
  A[ 1 ] * B[ 3 ] * C[ 3 ],
  A[ 2 ] * B[ 0 ] * C[ 0 ],
  A[ 2 ] * B[ 0 ] * C[ 1 ],
  A[ 2 ] * B[ 0 ] * C[ 2 ],
  A[ 2 ] * B[ 0 ] * C[ 3 ],
  A[ 2 ] * B[ 1 ] * C[ 0 ],
  A[ 2 ] * B[ 1 ] * C[ 1 ],
  A[ 2 ] * B[ 1 ] * C[ 2 ],
  A[ 2 ] * B[ 1 ] * C[ 3 ],
  A[ 2 ] * B[ 2 ] * C[ 0 ],
  A[ 2 ] * B[ 2 ] * C[ 1 ],
  A[ 2 ] * B[ 2 ] * C[ 2 ],
  A[ 2 ] * B[ 2 ] * C[ 3 ],
  A[ 2 ] * B[ 3 ] * C[ 0 ],
  A[ 2 ] * B[ 3 ] * C[ 1 ],
  A[ 2 ] * B[ 3 ] * C[ 2 ],
  A[ 2 ] * B[ 3 ] * C[ 3 ],
  A[ 3 ] * B[ 0 ] * C[ 0 ],
  A[ 3 ] * B[ 0 ] * C[ 1 ],
  A[ 3 ] * B[ 0 ] * C[ 2 ],
  A[ 3 ] * B[ 0 ] * C[ 3 ],
  A[ 3 ] * B[ 1 ] * C[ 0 ],
  A[ 3 ] * B[ 1 ] * C[ 1 ],
  A[ 3 ] * B[ 1 ] * C[ 2 ],
  A[ 3 ] * B[ 1 ] * C[ 3 ],
  A[ 3 ] * B[ 2 ] * C[ 0 ],
  A[ 3 ] * B[ 2 ] * C[ 1 ],
  A[ 3 ] * B[ 2 ] * C[ 2 ],
  A[ 3 ] * B[ 2 ] * C[ 3 ],
  A[ 3 ] * B[ 3 ] * C[ 0 ],
  A[ 3 ] * B[ 3 ] * C[ 1 ],
  A[ 3 ] * B[ 3 ] * C[ 2 ],
  A[ 3 ] * B[ 3 ] * C[ 3 ] };
  
  return(Wn[i-1]); 
}


//! LinearTetra, MINIELEMENT, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearTetra,MINIELEMENT,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d P1(0.0,0.0,0.0);
  static const point3d P2(1.0,0.0,0.0);
  static const point3d P3(0.0,1.0,0.0);
  static const point3d P4(0.0,0.0,1.0);
  static const point3d PM = (P1+P2+P3+P4)/4.0;
  
  static const point3d Yn[4] = {
  (P2+P3+P4+PM)/4.0,
  (P1+P3+P4+PM)/4.0,
  (P1+P2+P4+PM)/4.0,
  (P1+P2+P3+PM)/4.0 };
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearTetra,MINIELEMENT,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[4] = {
    1.0/24.0,
    1.0/24.0,
    1.0/24.0,
    1.0/24.0 };

  return(Wn[i-1]); 
}


//! LinearTetra, P0DUAL, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearTetra,P0DUAL,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d P1(0.0,0.0,0.0);
  static const point3d P2(1.0,0.0,0.0);
  static const point3d P3(0.0,1.0,0.0);
  static const point3d P4(0.0,0.0,1.0);
  static const point3d PM = (P1+P2+P3+P4)/4.0;
  
  static const point3d Yn[4] = {
    (P1+PM)/2.0,
    (P2+PM)/2.0,
    (P3+PM)/2.0,
    (P4+PM)/2.0 };
    
  return(Yn[i-1]); 
}

Real
intStaticCard<linearTetra,P0DUAL,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[4] = {
    1.0/24.0,
    1.0/24.0,
    1.0/24.0,
    1.0/24.0 };

  return(Wn[i-1]); 
}


//! LinearHexa, STANDARD, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearHexa,STANDARD,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = 0.0;
  static const Real p2 = 0.5;
  static const Real p3 = 1.0;
  
  static const point3d Yn[27] = {
  point3d(p1,p1,p1), point3d(p1,p2,p1), point3d(p1,p3,p1),
  point3d(p2,p1,p1), point3d(p2,p2,p1), point3d(p2,p3,p1),
  point3d(p3,p1,p1), point3d(p3,p2,p1), point3d(p3,p3,p1),
  
  point3d(p1,p1,p2), point3d(p1,p2,p2), point3d(p1,p3,p2),
  point3d(p2,p1,p2), point3d(p2,p2,p2), point3d(p2,p3,p2),
  point3d(p3,p1,p2), point3d(p3,p2,p2), point3d(p3,p3,p2),
  
  point3d(p1,p1,p3), point3d(p1,p2,p3), point3d(p1,p3,p3),
  point3d(p2,p1,p3), point3d(p2,p2,p3), point3d(p2,p3,p3),
  point3d(p3,p1,p3), point3d(p3,p2,p3), point3d(p3,p3,p3)};
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearHexa,STANDARD,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real w1 = 1.0/6.0;
  static const Real w2 = 2.0/3.0;
  static const Real w3 = 1.0/6.0;
  
  static const Real Wn[27] = {
  w1 * w1 * w1,  w1 * w2 * w1,  w1 * w3 * w1,
  w2 * w1 * w1,  w2 * w2 * w1,  w2 * w3 * w1,
  w3 * w1 * w1,  w3 * w2 * w1,  w3 * w3 * w1,
  
  w1 * w1 * w2,  w1 * w2 * w2,  w1 * w3 * w2,
  w2 * w1 * w2,  w2 * w2 * w2,  w2 * w3 * w2,
  w3 * w1 * w2,  w3 * w2 * w2,  w3 * w3 * w2,
  
  w1 * w1 * w3,  w1 * w2 * w3,  w1 * w3 * w3,
  w2 * w1 * w3,  w2 * w2 * w3,  w2 * w3 * w3,
  w3 * w1 * w3,  w3 * w2 * w3,  w3 * w3 * w3};
  
  return(Wn[i-1]);
}


//! LinearQuad, STANDARD, Precision 2---------------------------------------------------------------
point3d
intStaticCard<linearHexa,STANDARD,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = 0.0;
  static const Real p2 = 1.0/2.0 - sqrt(1.0/5.0);
  static const Real p3 = 1.0/2.0 + sqrt(1.0/5.0);
  static const Real p4 = 1.0;
  
  static const point3d Yn[64] = {
  point3d(p1,p1,p1), point3d(p1,p2,p1), point3d(p1,p3,p1), point3d(p1,p4,p1),
  point3d(p2,p1,p1), point3d(p2,p2,p1), point3d(p2,p3,p1), point3d(p2,p4,p1),
  point3d(p3,p1,p1), point3d(p3,p2,p1), point3d(p3,p3,p1), point3d(p3,p4,p1),
  point3d(p4,p1,p1), point3d(p4,p2,p1), point3d(p4,p3,p1), point3d(p4,p4,p1),
  
  point3d(p1,p1,p2), point3d(p1,p2,p2), point3d(p1,p3,p2), point3d(p1,p4,p2),
  point3d(p2,p1,p2), point3d(p2,p2,p2), point3d(p2,p3,p2), point3d(p2,p4,p2),
  point3d(p3,p1,p2), point3d(p3,p2,p2), point3d(p3,p3,p2), point3d(p3,p4,p2),
  point3d(p4,p1,p2), point3d(p4,p2,p2), point3d(p4,p3,p2), point3d(p4,p4,p2),
  
  point3d(p1,p1,p3), point3d(p1,p2,p3), point3d(p1,p3,p3), point3d(p1,p4,p3),
  point3d(p2,p1,p3), point3d(p2,p2,p3), point3d(p2,p3,p3), point3d(p2,p4,p3),
  point3d(p3,p1,p3), point3d(p3,p2,p3), point3d(p3,p3,p3), point3d(p3,p4,p3),
  point3d(p4,p1,p3), point3d(p4,p2,p3), point3d(p4,p3,p3), point3d(p4,p4,p3),
  
  point3d(p1,p1,p4), point3d(p1,p2,p4), point3d(p1,p3,p4), point3d(p1,p4,p4),
  point3d(p2,p1,p4), point3d(p2,p2,p4), point3d(p2,p3,p4), point3d(p2,p4,p4),
  point3d(p3,p1,p4), point3d(p3,p2,p4), point3d(p3,p3,p4), point3d(p3,p4,p4),
  point3d(p4,p1,p4), point3d(p4,p2,p4), point3d(p4,p3,p4), point3d(p4,p4,p4)};
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearHexa,STANDARD,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real w1 = 1.0/12.0;
  static const Real w2 = 5.0/12.0;
  static const Real w3 = 5.0/12.0;
  static const Real w4 = 1.0/12.0;
  
  static const Real Wn[64] = {
  w1 * w1 * w1,  w1 * w2 * w1,  w1 * w3 * w1,  w1 * w4 * w1,
  w2 * w1 * w1,  w2 * w2 * w1,  w2 * w3 * w1,  w2 * w4 * w1,
  w3 * w1 * w1,  w3 * w2 * w1,  w3 * w3 * w1,  w3 * w4 * w1,
  w4 * w1 * w1,  w4 * w2 * w1,  w4 * w3 * w1,  w4 * w4 * w1,
  
  w1 * w1 * w2,  w1 * w2 * w2,  w1 * w3 * w2,  w1 * w4 * w2,
  w2 * w1 * w2,  w2 * w2 * w2,  w2 * w3 * w2,  w2 * w4 * w2,
  w3 * w1 * w2,  w3 * w2 * w2,  w3 * w3 * w2,  w3 * w4 * w2,
  w4 * w1 * w2,  w4 * w2 * w2,  w4 * w3 * w2,  w4 * w4 * w2,
  
  w1 * w1 * w3,  w1 * w2 * w3,  w1 * w3 * w3,  w1 * w4 * w3,
  w2 * w1 * w3,  w2 * w2 * w3,  w2 * w3 * w3,  w2 * w4 * w3,
  w3 * w1 * w3,  w3 * w2 * w3,  w3 * w3 * w3,  w3 * w4 * w3,
  w4 * w1 * w3,  w4 * w2 * w3,  w4 * w3 * w3,  w4 * w4 * w3,
  
  w1 * w1 * w4,  w1 * w2 * w4,  w1 * w3 * w4,  w1 * w4 * w4,
  w2 * w1 * w4,  w2 * w2 * w4,  w2 * w3 * w4,  w2 * w4 * w4,
  w3 * w1 * w4,  w3 * w2 * w4,  w3 * w3 * w4,  w3 * w4 * w4,
  w4 * w1 * w4,  w4 * w2 * w4,  w4 * w3 * w4,  w4 * w4 * w4 };
  
  return(Wn[i-1]);
}


//! LinearQuad, STANDARD, Precision 3---------------------------------------------------------------
point3d
intStaticCard<linearHexa,STANDARD,3>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real p1 = 0.0;
  static const Real p2 = 1.0/2.0 - sqrt(3.0/7.0);
  static const Real p3 = 0.5;
  static const Real p4 = 1.0/2.0 + sqrt(3.0/7.0);
  static const Real p5 = 1.0;
  
  static const point3d Yn[125] = {
  point3d(p1,p1,p1), point3d(p1,p2,p1), point3d(p1,p3,p1), point3d(p1,p4,p1), point3d(p1,p5,p1),
  point3d(p2,p1,p1), point3d(p2,p2,p1), point3d(p2,p3,p1), point3d(p2,p4,p1), point3d(p2,p5,p1),
  point3d(p3,p1,p1), point3d(p3,p2,p1), point3d(p3,p3,p1), point3d(p3,p4,p1), point3d(p3,p5,p1),
  point3d(p4,p1,p1), point3d(p4,p2,p1), point3d(p4,p3,p1), point3d(p4,p4,p1), point3d(p4,p5,p1),
  point3d(p5,p1,p1), point3d(p5,p2,p1), point3d(p5,p3,p1), point3d(p5,p4,p1), point3d(p5,p5,p1),
  
  point3d(p1,p1,p2), point3d(p1,p2,p2), point3d(p1,p3,p2), point3d(p1,p4,p2), point3d(p1,p5,p2),
  point3d(p2,p1,p2), point3d(p2,p2,p2), point3d(p2,p3,p2), point3d(p2,p4,p2), point3d(p2,p5,p2),
  point3d(p3,p1,p2), point3d(p3,p2,p2), point3d(p3,p3,p2), point3d(p3,p4,p2), point3d(p3,p5,p2),
  point3d(p4,p1,p2), point3d(p4,p2,p2), point3d(p4,p3,p2), point3d(p4,p4,p2), point3d(p4,p5,p2),
  point3d(p5,p1,p2), point3d(p5,p2,p2), point3d(p5,p3,p2), point3d(p5,p4,p2), point3d(p5,p5,p2),
  
  point3d(p1,p1,p3), point3d(p1,p2,p3), point3d(p1,p3,p3), point3d(p1,p4,p3), point3d(p1,p5,p3),
  point3d(p2,p1,p3), point3d(p2,p2,p3), point3d(p2,p3,p3), point3d(p2,p4,p3), point3d(p2,p5,p3),
  point3d(p3,p1,p3), point3d(p3,p2,p3), point3d(p3,p3,p3), point3d(p3,p4,p3), point3d(p3,p5,p3),
  point3d(p4,p1,p3), point3d(p4,p2,p3), point3d(p4,p3,p3), point3d(p4,p4,p3), point3d(p4,p5,p3),
  point3d(p5,p1,p3), point3d(p5,p2,p3), point3d(p5,p3,p3), point3d(p5,p4,p3), point3d(p5,p5,p3),
  
  point3d(p1,p1,p4), point3d(p1,p2,p4), point3d(p1,p3,p4), point3d(p1,p4,p4), point3d(p1,p5,p4),
  point3d(p2,p1,p4), point3d(p2,p2,p4), point3d(p2,p3,p4), point3d(p2,p4,p4), point3d(p2,p5,p4),
  point3d(p3,p1,p4), point3d(p3,p2,p4), point3d(p3,p3,p4), point3d(p3,p4,p4), point3d(p3,p5,p4),
  point3d(p4,p1,p4), point3d(p4,p2,p4), point3d(p4,p3,p4), point3d(p4,p4,p4), point3d(p4,p5,p4),
  point3d(p5,p1,p4), point3d(p5,p2,p4), point3d(p5,p3,p4), point3d(p5,p4,p4), point3d(p5,p5,p4),
  
  point3d(p1,p1,p5), point3d(p1,p2,p5), point3d(p1,p3,p5), point3d(p1,p4,p5), point3d(p1,p5,p5),
  point3d(p2,p1,p5), point3d(p2,p2,p5), point3d(p2,p3,p5), point3d(p2,p4,p5), point3d(p2,p5,p5),
  point3d(p3,p1,p5), point3d(p3,p2,p5), point3d(p3,p3,p5), point3d(p3,p4,p5), point3d(p3,p5,p5),
  point3d(p4,p1,p5), point3d(p4,p2,p5), point3d(p4,p3,p5), point3d(p4,p4,p5), point3d(p4,p5,p5),
  point3d(p5,p1,p5), point3d(p5,p2,p5), point3d(p5,p3,p5), point3d(p5,p4,p5), point3d(p5,p5,p5) };
    
  return(Yn[i-1]);  
}

Real
intStaticCard<linearHexa,STANDARD,3>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real w1 =  1.0 /  20.0;
  static const Real w2 = 49.0 / 180.0;
  static const Real w3 = 32.0 /  90.0;
  static const Real w4 = 49.0 / 180.0;
  static const Real w5 =  1.0 /  20.0;
  
  static const Real Wn[125] = {
  w1 * w1 * w1,  w1 * w2 * w1,  w1 * w3 * w1,  w1 * w4 * w1,  w1 * w5 * w1,
  w2 * w1 * w1,  w2 * w2 * w1,  w2 * w3 * w1,  w2 * w4 * w1,  w2 * w5 * w1,
  w3 * w1 * w1,  w3 * w2 * w1,  w3 * w3 * w1,  w3 * w4 * w1,  w3 * w5 * w1,
  w4 * w1 * w1,  w4 * w2 * w1,  w4 * w3 * w1,  w4 * w4 * w1,  w4 * w5 * w1,
  w5 * w1 * w1,  w5 * w2 * w1,  w5 * w3 * w1,  w5 * w4 * w1,  w5 * w5 * w1,
  
  w1 * w1 * w2,  w1 * w2 * w2,  w1 * w3 * w2,  w1 * w4 * w2,  w1 * w5 * w2,
  w2 * w1 * w2,  w2 * w2 * w2,  w2 * w3 * w2,  w2 * w4 * w2,  w2 * w5 * w2,
  w3 * w1 * w2,  w3 * w2 * w2,  w3 * w3 * w2,  w3 * w4 * w2,  w3 * w5 * w2,
  w4 * w1 * w2,  w4 * w2 * w2,  w4 * w3 * w2,  w4 * w4 * w2,  w4 * w5 * w2,
  w5 * w1 * w2,  w5 * w2 * w2,  w5 * w3 * w2,  w5 * w4 * w2,  w5 * w5 * w2,
  
  w1 * w1 * w3,  w1 * w2 * w3,  w1 * w3 * w3,  w1 * w4 * w3,  w1 * w5 * w3,
  w2 * w1 * w3,  w2 * w2 * w3,  w2 * w3 * w3,  w2 * w4 * w3,  w2 * w5 * w3,
  w3 * w1 * w3,  w3 * w2 * w3,  w3 * w3 * w3,  w3 * w4 * w3,  w3 * w5 * w3,
  w4 * w1 * w3,  w4 * w2 * w3,  w4 * w3 * w3,  w4 * w4 * w3,  w4 * w5 * w3,
  w5 * w1 * w3,  w5 * w2 * w3,  w5 * w3 * w3,  w5 * w4 * w3,  w5 * w5 * w3,
  
  w1 * w1 * w4,  w1 * w2 * w4,  w1 * w3 * w4,  w1 * w4 * w4,  w1 * w5 * w4,
  w2 * w1 * w4,  w2 * w2 * w4,  w2 * w3 * w4,  w2 * w4 * w4,  w2 * w5 * w4,
  w3 * w1 * w4,  w3 * w2 * w4,  w3 * w3 * w4,  w3 * w4 * w4,  w3 * w5 * w4,
  w4 * w1 * w4,  w4 * w2 * w4,  w4 * w3 * w4,  w4 * w4 * w4,  w4 * w5 * w4,
  w5 * w1 * w4,  w5 * w2 * w4,  w5 * w3 * w4,  w5 * w4 * w4,  w5 * w5 * w4,
  
  w1 * w1 * w5,  w1 * w2 * w5,  w1 * w3 * w5,  w1 * w4 * w5,  w1 * w5 * w5,
  w2 * w1 * w5,  w2 * w2 * w5,  w2 * w3 * w5,  w2 * w4 * w5,  w2 * w5 * w5,
  w3 * w1 * w5,  w3 * w2 * w5,  w3 * w3 * w5,  w3 * w4 * w5,  w3 * w5 * w5,
  w4 * w1 * w5,  w4 * w2 * w5,  w4 * w3 * w5,  w4 * w4 * w5,  w4 * w5 * w5,
  w5 * w1 * w5,  w5 * w2 * w5,  w5 * w3 * w5,  w5 * w4 * w5,  w5 * w5 * w5 };
  
  return(Wn[i-1]);
}


//! LinearWedge, STANDARD, Precision 1---------------------------------------------------------------
point3d
intStaticCard<linearWedge,STANDARD,1>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const point3d Yn[1] = {point3d(1.0/3.0, 1.0/3.0, 0.5)};
  return(Yn[i-1]);
}

Real
intStaticCard<linearWedge,STANDARD,1>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[1] = {0.5};
  return(Wn[i-1]);
}


//! LinearWedge, STANDARD, Precision 2---------------------------------------------------------------
point3d
intStaticCard<linearWedge,STANDARD,2>::
getYn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real q2ptx1 = ( 1 - sqrt( 1. / 3. ) ) / 2., q2ptx2 = ( 1 + sqrt( 1. / 3. ) ) / 2.;
  
  static const point3d Yn[6] = {
    point3d(0.5, 0.0, q2ptx1),
    point3d(0.0, 0.5, q2ptx1),
    point3d(0.5, 0.5, q2ptx1),
    point3d(0.5, 0.0, q2ptx2),
    point3d(0.0, 0.5, q2ptx2),
    point3d(0.5, 0.5, q2ptx2)
  };
  
  return(Yn[i-1]);
}

Real
intStaticCard<linearWedge,STANDARD,2>::
getWn(const UInt & i)
{
  assert(i >= 1);
  assert(i <= N);
  
  static const Real Wn[6] = {1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0,
                             1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0 };
  return(Wn[i-1]);
}

