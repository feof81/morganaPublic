/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "intOscSymplexSupport.h"
#include "../morganaDofs/komplex.h"
#include <map>

const Real intOscSymplexSupport::ratioToll = 1.0e-10;

komplex
intOscSymplexSupport::
intS0_1d(const point3d & K,
         const point3d & P1,
         const point3d & P2)
{
  //Exit for null wave vector----------------------------------------
  if(point3d::norm2(K) == 0)
  {
    Real lung = point3d::norm2(P2-P1);
    return(komplex(1.0,0.0) * lung);
  }
  
  //Geometry---------------------------------------------------------
  point3d D = P2 - P1;
  Real   kd = point3d::dot(K,D);
  
  //Compute----------------------------------------------------------
  if(kd == 0.0) { return( komplex::iexp(point3d::dot(K,P1)) * point3d::norm2(D)); }
  else          { return(  komplex(0.0,-1.0) *
                           komplex::iexp(point3d::dot(K,P1)) *
                          (komplex::iexp(kd) - komplex(1.0,0.0)) * (point3d::norm2(D) / kd) ); }
}

komplex
intOscSymplexSupport::
intS0_2d(const point3d & K,
         const point3d & P1,
         const point3d & P2,
         const point3d & P3)
{
  //Exit for null wave vector----------------------------------------
  Real normK = point3d::norm2(K);
  
  if(normK <= 1.0e-3)
  {
    Real area = point3d::norm2((P3-P1) ^ (P2-P1)) / 2.0;
    return(komplex(1.0,0.0) * area);
  }
  
  //Geometry---------------------------------------------------------
  point3d D1 = P2 - P1;
  point3d D2 = P3 - P2;
  point3d D3 = P1 - P3;
  
  point3d Nf = (P2-P1) ^ (P3-P1);
  Real  Area = point3d::norm2(Nf) / 2.0;
  Nf /= (Area * 2.0);
  
  point3d Nd1 = D1 ^ Nf; Nd1 /= point3d::norm2(Nd1);
  point3d Nd2 = D2 ^ Nf; Nd2 /= point3d::norm2(Nd2);
  point3d Nd3 = D3 ^ Nf; Nd3 /= point3d::norm2(Nd3);
  
  point3d Kf = K - Nf * point3d::dot(K,Nf);
  
  //Compute----------------------------------------------------------
  Real kd1 = point3d::dot(K,D1);
  Real kd2 = point3d::dot(K,D2);
  Real kd3 = point3d::dot(K,D3);
  
  komplex out;
  
  if((point3d::norm2(Kf) <= (ratioToll * normK)) || ((kd1 == 0.0) && (kd2 == 0.0) && (kd3 == 0.0)))
  {
    out += komplex::iexp(point3d::dot(K,P1)) * Area;
  }
  else
  {
    //Edge 1
    if(kd1 == 0.0) 
    { out -= komplex::iexp(point3d::dot(K,P1)) *
	     komplex(0.0,1.0) *
	    (point3d::dot(K,Nd1) / point3d::dot(Kf,Kf)) *
	     point3d::norm2(D1);
    }
    else
    {
      out -= komplex::iexp(point3d::dot(K,P1)) *
	   ((komplex::iexp(kd1) - komplex(1.0,0.0)) / kd1) *
	    (point3d::dot(K,Nd1) / point3d::dot(Kf,Kf)) *
	     point3d::norm2(D1);
    }
    
    //Edge 2
    if(kd2 == 0.0) 
    { out -= komplex::iexp(point3d::dot(K,P2)) *
	     komplex(0.0,1.0) *
	    (point3d::dot(K,Nd2) / point3d::dot(Kf,Kf)) *
	     point3d::norm2(D2);
    }
    else
    {
      out -= komplex::iexp(point3d::dot(K,P2)) *
	   ((komplex::iexp(kd2) - komplex(1.0,0.0)) / kd2) *
	    (point3d::dot(K,Nd2) / point3d::dot(Kf,Kf)) *
	     point3d::norm2(D2);
    }
    
    //Edge 3
    if(kd3 == 0.0) 
    { out -= komplex::iexp(point3d::dot(K,P3)) *
	     komplex(0.0,1.0) *
	    (point3d::dot(K,Nd3) / point3d::dot(Kf,Kf)) *
	     point3d::norm2(D3);
    }
    else
    {
      out -= komplex::iexp(point3d::dot(K,P3)) *
	   ((komplex::iexp(kd3) - komplex(1.0,0.0)) / kd3) *
	    (point3d::dot(K,Nd3) / point3d::dot(Kf,Kf)) *
	     point3d::norm2(D3);
    }
  }
  
  return(out);
}

komplex
intOscSymplexSupport::
intS0_3d(const point3d & K,
         const point3d & P1,
         const point3d & P2,
         const point3d & P3,
         const point3d & P4)
{
  //Exit for null wave vector----------------------------------------
  Real normK = point3d::norm2(K);
  
  if(normK <= 1.0e-3)
  {
    Real vol = point3d::dot(P4-P1, (P2-P1) ^ (P3-P1)) / 6.0;
    return(komplex(1.0,0.0) * vol);
  }
  
  //Geometry---------------------------------------------------------
  point3d Nf1 = (P3-P2) ^ (P4-P2);  Real Area1 = point3d::norm2(Nf1) / 2.0;  Nf1 /= (Area1 * 2.0);
  point3d Nf2 = (P4-P1) ^ (P3-P1);  Real Area2 = point3d::norm2(Nf2) / 2.0;  Nf2 /= (Area2 * 2.0);
  point3d Nf3 = (P2-P1) ^ (P4-P1);  Real Area3 = point3d::norm2(Nf3) / 2.0;  Nf3 /= (Area3 * 2.0);
  point3d Nf4 = (P3-P1) ^ (P2-P1);  Real Area4 = point3d::norm2(Nf4) / 2.0;  Nf4 /= (Area4 * 2.0);
  
  point3d D11 = P3 - P2;  const point3d & P11 = P2;
  point3d D12 = P4 - P3;  const point3d & P12 = P3;
  point3d D13 = P2 - P4;  const point3d & P13 = P4;
  
  point3d D21 = P4 - P1;  const point3d & P21 = P1;
  point3d D22 = P3 - P4;  const point3d & P22 = P4;
  point3d D23 = P1 - P3;  const point3d & P23 = P3;
  
  point3d D31 = P4 - P2;  const point3d & P31 = P2; 
  point3d D32 = P1 - P4;  const point3d & P32 = P4;
  point3d D33 = P2 - P1;  const point3d & P33 = P1;
  
  point3d D41 = P3 - P1;  const point3d & P41 = P1;
  point3d D42 = P2 - P3;  const point3d & P42 = P3;
  point3d D43 = P1 - P2;  const point3d & P43 = P2;
  
  point3d N11 = D11 ^ Nf1;  N11 /= point3d::norm2(N11); 
  point3d N12 = D12 ^ Nf1;  N12 /= point3d::norm2(N12);
  point3d N13 = D13 ^ Nf1;  N13 /= point3d::norm2(N13);
  
  point3d N21 = D21 ^ Nf2;  N21 /= point3d::norm2(N21); 
  point3d N22 = D22 ^ Nf2;  N22 /= point3d::norm2(N22);
  point3d N23 = D23 ^ Nf2;  N23 /= point3d::norm2(N23);
  
  point3d N31 = D31 ^ Nf3;  N31 /= point3d::norm2(N31); 
  point3d N32 = D32 ^ Nf3;  N32 /= point3d::norm2(N32);
  point3d N33 = D33 ^ Nf3;  N33 /= point3d::norm2(N33);
  
  point3d N41 = D41 ^ Nf4;  N41 /= point3d::norm2(N41); 
  point3d N42 = D42 ^ Nf4;  N42 /= point3d::norm2(N42);
  point3d N43 = D43 ^ Nf4;  N43 /= point3d::norm2(N43);
  
  point3d Kf1 = K - Nf1 * point3d::dot(K,Nf1);
  point3d Kf2 = K - Nf2 * point3d::dot(K,Nf2);
  point3d Kf3 = K - Nf3 * point3d::dot(K,Nf3);
  point3d Kf4 = K - Nf4 * point3d::dot(K,Nf4);
  
  //Compute----------------------------------------------------------
  Real kd11 = point3d::dot(K,D11);
  Real kd12 = point3d::dot(K,D12);
  Real kd13 = point3d::dot(K,D13);
  
  Real kd21 = point3d::dot(K,D21);
  Real kd22 = point3d::dot(K,D22);
  Real kd23 = point3d::dot(K,D23);
  
  Real kd31 = point3d::dot(K,D31);
  Real kd32 = point3d::dot(K,D32);
  Real kd33 = point3d::dot(K,D33);
  
  Real kd41 = point3d::dot(K,D41);
  Real kd42 = point3d::dot(K,D42);
  Real kd43 = point3d::dot(K,D43);
  
  komplex out;
  komplex md, mf, me;
  
  
  //FACE 1-----------------------------------------------------------
  me = (point3d::dot(K,Nf1) / (normK * normK));
  
  if((point3d::norm2(Kf1) <= (ratioToll * normK)) || ((kd11 == 0.0) && (kd12 == 0.0) && (kd13 == 0.0)))
  {
    out += me * Area1 * komplex::iexp(point3d::dot(K,P2));
  }
  else
  {
    //Face 1 - edge 1
    if(kd11 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd11) - komplex(1.0,0.0)) / kd11; }
    
    mf   = point3d::dot(K,N11) / point3d::dot(Kf1,Kf1);
    out -= komplex::iexp(point3d::dot(K,P11)) * md * mf * me * point3d::norm2(D11);
    
    //Face 1 - edge 2
    if(kd12 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd12) - komplex(1.0,0.0)) / kd12; }
    
    mf   = point3d::dot(K,N12) / point3d::dot(Kf1,Kf1);
    out -= komplex::iexp(point3d::dot(K,P12)) * md * mf * me * point3d::norm2(D12);
    
    //Face 1 - edge 3
    if(kd13 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd13) - komplex(1.0,0.0)) / kd13; }
    
    mf   = point3d::dot(K,N13) / point3d::dot(Kf1,Kf1);
    out -= komplex::iexp(point3d::dot(K,P13)) * md * mf * me * point3d::norm2(D13);
  }

  
  //FACE 2-----------------------------------------------------------
  me = (point3d::dot(K,Nf2) / (normK * normK));
  
  if((point3d::norm2(Kf2) <= (ratioToll * normK)) || ((kd21 == 0.0) && (kd22 == 0.0) && (kd23 == 0.0)))
  {
    out += me * Area2 * komplex::iexp(point3d::dot(K,P1));
  }
  else
  {
    //Face 2 - edge 1
    if(kd21 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd21) - komplex(1.0,0.0)) / kd21; }
    
    mf   = point3d::dot(K,N21) / point3d::dot(Kf2,Kf2);
    out -= komplex::iexp(point3d::dot(K,P21)) * md * mf * me * point3d::norm2(D21);
    
    //Face 2 - edge 2
    if(kd22 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd22) - komplex(1.0,0.0)) / kd22; }
    
    mf   = point3d::dot(K,N22) / point3d::dot(Kf2,Kf2);
    out -= komplex::iexp(point3d::dot(K,P22)) * md * mf * me * point3d::norm2(D22);
    
    //Face 2 - edge 3
    if(kd23 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd23) - komplex(1.0,0.0)) / kd23; }
    
    mf   = point3d::dot(K,N23) / point3d::dot(Kf2,Kf2);
    out -= komplex::iexp(point3d::dot(K,P23)) * md * mf * me * point3d::norm2(D23);
  }
  
  
  //FACE 3-----------------------------------------------------------
  me = (point3d::dot(K,Nf3) / (normK * normK));
  
  if((point3d::norm2(Kf3) <= (ratioToll * normK)) || ((kd31 == 0.0) && (kd32 == 0.0) && (kd33 == 0.0)))
  {
    out += me * Area3 * komplex::iexp(point3d::dot(K,P1));
  }
  else
  {
    //Face 3 - edge 1
    if(kd31 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd31) - komplex(1.0,0.0)) / kd31; }
    
    mf   = point3d::dot(K,N31) / point3d::dot(Kf3,Kf3);
    out -= komplex::iexp(point3d::dot(K,P31)) * md * mf * me * point3d::norm2(D31);
    
    //Face 3 - edge 2
    if(kd32 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd32) - komplex(1.0,0.0)) / kd32; }
    
    mf   = point3d::dot(K,N32) / point3d::dot(Kf3,Kf3);
    out -= komplex::iexp(point3d::dot(K,P32)) * md * mf * me * point3d::norm2(D32);
    
    //Face 3 - edge 3
    if(kd33 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd33) - komplex(1.0,0.0)) / kd33; }
    
    mf   = point3d::dot(K,N33) / point3d::dot(Kf3,Kf3);
    out -= komplex::iexp(point3d::dot(K,P33)) * md * mf * me * point3d::norm2(D33);
  }

  
  //FACE 4-----------------------------------------------------------
  me = (point3d::dot(K,Nf4) / (normK * normK));
  
  if((point3d::norm2(Kf4) <= (ratioToll * normK)) || ((kd41 == 0.0) && (kd42 == 0.0) && (kd43 == 0.0)))
  {
    out += me * Area4 * komplex::iexp(point3d::dot(K,P1));
  }
  else
  {
    //Face 4 - edge 1
    if(kd41 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd41) - komplex(1.0,0.0)) / kd41; }
    
    mf   = point3d::dot(K,N41) / point3d::dot(Kf4,Kf4);
    out -= komplex::iexp(point3d::dot(K,P41)) * md * mf * me * point3d::norm2(D41);
    
    //Face 4 - edge 2
    if(kd42 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd42) - komplex(1.0,0.0)) / kd42; }
    
    mf   = point3d::dot(K,N42) / point3d::dot(Kf4,Kf4);
    out -= komplex::iexp(point3d::dot(K,P42)) * md * mf * me * point3d::norm2(D42);
    
    //Face 4 - edge 3
    if(kd43 == 0.0) { md =  komplex(0.0,1.0); }
    else            { md = (komplex::iexp(kd43) - komplex(1.0,0.0)) / kd43; }
    
    mf   = point3d::dot(K,N43) / point3d::dot(Kf4,Kf4);
    out -= komplex::iexp(point3d::dot(K,P43)) * md * mf * me * point3d::norm2(D43);
  }

  return(out * komplex(0.0,-1.0));
}

sVect<intOscSymplexSupport::TETRA>
intOscSymplexSupport::
refineTetra(const UInt & nRef,const TETRA & tetra)
{
  //Alloc------------------------------------------------------------
  sVect<TETRA> oldTetra(1);
  oldTetra(1) = tetra;
  
  TETRA newTri;
  
  //Refine-----------------------------------------------------------
  point3d P1, P2, P3, P4, P5, P6, P7, P8, P9, P10;
  sVect<TETRA> newTetra;
  
  for(UInt s=1; s <= nRef; s++)
  {
    newTetra.clear();
    
    for(UInt i=1; i <= oldTetra.size(); ++i)
    {
      P1 = oldTetra(i).getPoint(1);
      P2 = oldTetra(i).getPoint(2);
      P3 = oldTetra(i).getPoint(3);
      P4 = oldTetra(i).getPoint(4);
      
      P5  = (P1 + P2) / 2.0;
      P6  = (P1 + P3) / 2.0;
      P7  = (P1 + P4) / 2.0;
      P8  = (P2 + P3) / 2.0;
      P9  = (P3 + P4) / 2.0;
      P10 = (P4 + P2) / 2.0;
      
      //el1
      newTri.setPoint(1,P1);
      newTri.setPoint(2,P7);
      newTri.setPoint(3,P5);
      newTri.setPoint(4,P6);
      newTetra.push_back(newTri);
      
      //el2
      newTri.setPoint(1,P7);
      newTri.setPoint(2,P4);
      newTri.setPoint(3,P10);
      newTri.setPoint(4,P9);
      newTetra.push_back(newTri);
      
      //el3
      newTri.setPoint(1,P5);
      newTri.setPoint(2,P10);
      newTri.setPoint(3,P2);
      newTri.setPoint(4,P8);
      newTetra.push_back(newTri);
      
      //el4
      newTri.setPoint(1,P6);
      newTri.setPoint(2,P9);
      newTri.setPoint(3,P8);
      newTri.setPoint(4,P3);
      newTetra.push_back(newTri);
      
      //el5
      newTri.setPoint(1,P5);
      newTri.setPoint(2,P6);
      newTri.setPoint(3,P7);
      newTri.setPoint(4,P9);
      newTetra.push_back(newTri);
      
      //el6
      newTri.setPoint(1,P5);
      newTri.setPoint(2,P7);
      newTri.setPoint(3,P10);
      newTri.setPoint(4,P9);
      newTetra.push_back(newTri);
      
      //el7
      newTri.setPoint(1,P5);
      newTri.setPoint(2,P10);
      newTri.setPoint(3,P8);
      newTri.setPoint(4,P9);
      newTetra.push_back(newTri);
      
      //el8
      newTri.setPoint(1,P5);
      newTri.setPoint(2,P8);
      newTri.setPoint(3,P6);
      newTri.setPoint(4,P9);
      newTetra.push_back(newTri);
    }
    
    oldTetra = newTetra;
  }
  
  return(oldTetra);
}

sVect<intOscSymplexSupport::TETRA>
intOscSymplexSupport::
refineTetra2(const UInt & nRef,const TETRA & tetra)
{
  //Alloc------------------------------------------------------------
  sVect<TETRA> oldTetra(1);
  oldTetra(1) = tetra;
  
  TETRA newTri;
  
  //Refine-----------------------------------------------------------
  UInt id;
  point3d P1, P2, P3, P4;
  point3d A1, A2, A3, A4;
  point3d B1, B2, B3, B4;
  point3d D1, D2, D3, D4, D5, D6;
  sVect<TETRA> newTetra;
  std::map<Real,UInt> edgeMap;
  std::map<Real,UInt>::iterator it;
  
  for(UInt s=1; s <= nRef; s++)
  {
    newTetra.clear();
    
    for(UInt i=1; i <= oldTetra.size(); ++i)
    {
      P1 = oldTetra(i).getPoint(1);
      P2 = oldTetra(i).getPoint(2);
      P3 = oldTetra(i).getPoint(3);
      P4 = oldTetra(i).getPoint(4);
      
      D1 = (P1+P2) / 2.0;
      D2 = (P2+P3) / 2.0;
      D3 = (P3+P1) / 2.0;
      D4 = (P4+P1) / 2.0;
      D5 = (P4+P2) / 2.0;
      D6 = (P4+P3) / 2.0;
      
      edgeMap.clear();
      edgeMap.insert(std::pair<Real,UInt>(point3d::norm2(P2 - P1),1));
      edgeMap.insert(std::pair<Real,UInt>(point3d::norm2(P3 - P2),2));
      edgeMap.insert(std::pair<Real,UInt>(point3d::norm2(P1 - P3),3));
      edgeMap.insert(std::pair<Real,UInt>(point3d::norm2(P4 - P1),4));
      edgeMap.insert(std::pair<Real,UInt>(point3d::norm2(P4 - P2),5));
      edgeMap.insert(std::pair<Real,UInt>(point3d::norm2(P4 - P3),6));
      
      it = edgeMap.end();
      it--;
      id = it->second;
      
      A1 = P1 * Real(id == 1) + P1 * Real(id == 2) + P1 * Real(id == 3) + P1 * Real(id == 4) + P1 * Real(id == 5) + P1 * Real(id == 6);
      A2 = D1 * Real(id == 1) + P2 * Real(id == 2) + P2 * Real(id == 3) + P2 * Real(id == 4) + P2 * Real(id == 5) + P2 * Real(id == 6);
      A3 = P3 * Real(id == 1) + D2 * Real(id == 2) + D3 * Real(id == 3) + P3 * Real(id == 4) + P3 * Real(id == 5) + P3 * Real(id == 6);
      A4 = P4 * Real(id == 1) + P4 * Real(id == 2) + P4 * Real(id == 3) + D4 * Real(id == 4) + D5 * Real(id == 5) + D6 * Real(id == 6);
      
      B1 = D1 * Real(id == 1) + P1 * Real(id == 2) + P2 * Real(id == 3) + D4 * Real(id == 4) + P1 * Real(id == 5) + P1 * Real(id == 6);
      B2 = P2 * Real(id == 1) + D2 * Real(id == 2) + P3 * Real(id == 3) + P2 * Real(id == 4) + D5 * Real(id == 5) + P2 * Real(id == 6);
      B3 = P3 * Real(id == 1) + P3 * Real(id == 2) + D3 * Real(id == 3) + P3 * Real(id == 4) + P3 * Real(id == 5) + D6 * Real(id == 6);
      B4 = P4 * Real(id == 1) + P4 * Real(id == 2) + P4 * Real(id == 3) + P4 * Real(id == 4) + P4 * Real(id == 5) + P4 * Real(id == 6);
      
      //el1
      newTri.setPoint(1,A1);
      newTri.setPoint(2,A2);
      newTri.setPoint(3,A3);
      newTri.setPoint(4,A4);
      newTetra.push_back(newTri);
      
      //el2
      newTri.setPoint(1,B1);
      newTri.setPoint(2,B2);
      newTri.setPoint(3,B3);
      newTri.setPoint(4,B4);
      newTetra.push_back(newTri);
    }
    
    oldTetra = newTetra;
  }
  
  return(oldTetra);
}

sVect<intOscSymplexSupport::TRIANGLE>
intOscSymplexSupport::
refineTriangle(const UInt & nRef, const TRIANGLE & triangle)
{
  //Alloc------------------------------------------------------------
  sVect<TRIANGLE> oldTriangle(1);
  oldTriangle(1) = triangle;
  
  TRIANGLE newTri;
  
  //Refine-----------------------------------------------------------
  point3d P1, P2, P3, P4, P5, P6;
  sVect<TRIANGLE> newTriangle;
  
  for(UInt s=1; s <= nRef; s++)
  {
    newTriangle.clear();
    
    for(UInt i=1; i <= oldTriangle.size(); ++i)
    {
      P1 = oldTriangle(i).getPoint(1);
      P2 = oldTriangle(i).getPoint(2);
      P3 = oldTriangle(i).getPoint(3);
      
      P4  = (P1 + P2) / 2.0;
      P5  = (P2 + P3) / 2.0;
      P6  = (P3 + P1) / 2.0;
      
      //el1
      newTri.setPoint(1,P1);
      newTri.setPoint(2,P4);
      newTri.setPoint(3,P6);
      newTriangle.push_back(newTri);
      
      //el2
      newTri.setPoint(1,P4);
      newTri.setPoint(2,P2);
      newTri.setPoint(3,P5);
      newTriangle.push_back(newTri);
      
      //el3
      newTri.setPoint(1,P5);
      newTri.setPoint(2,P3);
      newTri.setPoint(3,P6);
      newTriangle.push_back(newTri);
      
      //el4
      newTri.setPoint(1,P4);
      newTri.setPoint(2,P5);
      newTri.setPoint(3,P6);
      newTriangle.push_back(newTri);
    }
    
    oldTriangle = newTriangle;
  }
  
  return(oldTriangle);
}




/*
template<typename OPERATOR>
sVect<komplex>
intOscSymplexSupport::
intS1_1d(OPERATOR & op,
         const point3d & K,
         const point3d & P1,
         const point3d & P2,
         const point3d & Y1,
         const point3d & Y2)
{
  //Eval the function------------------------------------------------
  UInt matSize = op.numIndex_field() * op.numIndex_test();
  sVect<komplex> bufMat1(matSize), bufMat2(matSize);
  
  point3d YM = (Y1 + Y2) / 2.0;
  
  op.eval(Y1,bufMat1);
  op.eval(Y2,bufMat2);
}*/

/*
sVect<komplex>
intOscSymplexSupport::
intS1_1d(const sVect<Real>    & valF1,
         const sVect<Real>    & valF2,
         const sVect<point3d> & gradFM,
         const point3d & K,
         const point3d & P1,
         const point3d & P2)
{
  //Assert-----------------------------------------------------------
  assert(valF1.size() == gradFM.size());
  assert(valF2.size() == gradFM.size());
  
  //Geometry---------------------------------------------------------
  point3d D = P2 - P1;  Real length = point3d::norm2(D);  D /= length;
  Real   kd = point3d::dot(K,D);
  
  //Compute----------------------------------------------------------
  sVect<komplex> out(valF1.size());
  
  if(kd == 0)
  {    
    for(UInt i=1; i <= valF1.size(); ++i)
    {
      out(i) = (valF1(i) + valF2(i)) * (length / 2.0) * komplex::iexp(point3d::dot(K,P1));
    }
  }
  else
  {
    for(UInt i=1; i <= valF1.size(); ++i)
    {
      out(i) = komplex(0.0,-1.0) * (komplex::iexp(point3d::dot(K,P2)) * valF2(i) - komplex::iexp(point3d::dot(K,P1)) * valF1(i))
             + komplex(0.0,point3d::dot(gradFM(i),D)) * intS0_1d(K,P1,P2);
      
      out(i) /= kd;
    }
  }
  
  return(out);
}*/

/*sVect<komplex>
intOscSymplexSupport::
intS1_2d(const sVect<Real>    & valF1,
         const sVect<Real>    & valF2,
         const sVect<Real>    & valF3,
         const sVect<point3d> & gradFM,
         const point3d & K,
         const point3d & P1,
         const point3d & P2,
         const point3d & P3)
{
  //Assert-----------------------------------------------------------
  assert(valFM.size() == gradFM.size());
  
  //Geometry---------------------------------------------------------
  point3d PM = (P1 + P2 + P3) / 3.0;
  
  point3d D1 = P2 - P1;
  point3d D2 = P3 - P2;
  point3d D3 = P1 - P3;
  
  point3d Nf = (P2-P1) ^ (P3-P1);
  Real  Area = point3d::norm2(Nf) / 2.0;
  Nf /= (Area * 2.0);
  
  point3d Nd1 = D1 ^ Nf; Nd1 /= point3d::norm2(Nd1);
  point3d Nd2 = D2 ^ Nf; Nd2 /= point3d::norm2(Nd2);
  point3d Nd3 = D3 ^ Nf; Nd3 /= point3d::norm2(Nd3);
  
  point3d Kf = K - Nf * point3d::dot(K,Nf);
  
  //Compute----------------------------------------------------------
  sVect<komplex> out(valFM.size());
  
  if(point3d::norm2(Kf) == 0)
  {
    for(UInt i=1; i <= valFM.size(); ++i)
    {
      out(i) += komplex::iexp(point3d::dot(K,PM)) * valFM(i) * Area;
    }
  }
  else
  {
    for(UInt i=1; i <= valFM.size(); ++i)
    {
      out(i) -= intS1_1d(valF1(i),valF2(i),gradFM(i),K,P1,P2) * point3d::dot(Kf,Nd1) / point3d::dot(Kf,Kf);
      out(i) -= intS1_1d(valF2(i),valF3(i),gradFM(i),K,P1,P2) * point3d::dot(Kf,Nd1) / point3d::dot(Kf,Kf);
      
    }
  }
  
  return(out * komplex(0.0,1.0));
}*/
