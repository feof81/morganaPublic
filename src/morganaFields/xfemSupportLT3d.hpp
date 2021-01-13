/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef XFEMSUPPORTLT3D_HPP
#define XFEMSUPPORTLT3D_HPP

#include "feXFEM_LT3d_card.hpp"


enum Xvolumes {Xtetra, Xwedge};


/*! Simple support for the three dimensional extended finite elements on the linear tetra.
Using the local representation of the level set function the implemented methods can decompose
the tetra into wedges and tetra. */
template<typename PMAPTYPE>
class xfemSupportLT3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feXFEM_LT3d_card            FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feDynamicDofCard3d          DOFCARD;
    //@}
    
    /*! @name Constructors */ //@{
  public:
    xfemSupportLT3d();
    //@}
    
    /*! @name Splitting functions */ //@{
  public:
    void caseB(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    void caseC1(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    void caseC2(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    void caseD1(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    void caseD2(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    void caseD3(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    void split(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    //@}
    
    /*! @name Volume functions */ //@{
  public:
    Real tetraVolume(const sVect<point3d> & points) const;
    Real wedgeVolume(const sVect<point3d> & points) const;
    sVect<Real> volumes(const FECARD & FeCard, const ELCARD & ElCard) const;
    sVect<Real> volumes(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const;
    //@}
    
    /*! @name Info functions */ //@{
  public:
    bool isSplit(const FECARD & FeCard, const ELCARD & ElCard) const;
    sVect<point3d> getBarycenters(const FECARD & FeCard, const ELCARD & ElCard) const;
    //@}
};

template<typename PMAPTYPE>
xfemSupportLT3d<PMAPTYPE>::
xfemSupportLT3d()
{ 
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
caseB(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{ 
  UInt nodeX = UInt(FeCard.getPhi(1) > 0.5) * 1
             + UInt(FeCard.getPhi(2) > 0.5) * 2
	     + UInt(FeCard.getPhi(3) > 0.5) * 3
	     + UInt(FeCard.getPhi(4) > 0.5) * 4;
	     
  assert(nodeX >= 1); assert(nodeX <= 4);
  
  UInt faceX = GEOSHAPE::faceOppositeToNode(nodeX);
  
  UInt nodeO1 = GEOSHAPE::faceToPoint(faceX,1);
  UInt nodeO2 = GEOSHAPE::faceToPoint(faceX,2);
  UInt nodeO3 = GEOSHAPE::faceToPoint(faceX,3);
  
  Real csi1 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO1) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO1) - FeCard.getPhi(nodeX)) != 0.0);
  Real csi2 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO2) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO2) - FeCard.getPhi(nodeX)) != 0.0);
  Real csi3 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO3) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO3) - FeCard.getPhi(nodeX)) != 0.0);
  
  point3d P1 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO1) - ElCard.getNode(nodeX)) * csi1;
  point3d P2 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO2) - ElCard.getNode(nodeX)) * csi2;
  point3d P3 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO3) - ElCard.getNode(nodeX)) * csi3;
  
  //Build the tetra
  sVect<point3d> tetra(4);
  tetra(1) = P1;
  tetra(2) = P3;
  tetra(3) = P2;
  tetra(4) = ElCard.getNode(nodeX);
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
  
  //Build the wedge
  sVect<point3d> wedge(6);
  wedge(1) = ElCard.getNode(nodeO1);
  wedge(2) = ElCard.getNode(nodeO3);
  wedge(3) = ElCard.getNode(nodeO2);
  wedge(4) = P1;
  wedge(5) = P3;
  wedge(6) = P2;
  
  volLabels.push_back(Xwedge);
  points.push_back(wedge);
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
caseC1(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{
  sVect<UInt> nodeX, nodeO;
  
  if(FeCard.getPhi(1) > 0.5) { nodeX.push_back(1); }
  else                       { nodeO.push_back(1); }
  
  if(FeCard.getPhi(2) > 0.5) { nodeX.push_back(2); }
  else                       { nodeO.push_back(2); }
  
  if(FeCard.getPhi(3) > 0.5) { nodeX.push_back(3); }
  else                       { nodeO.push_back(3); }
  
  if(FeCard.getPhi(4) > 0.5) { nodeX.push_back(4); }
  else                       { nodeO.push_back(4); }
  
  
  Real csi11 = (0.5 - FeCard.getPhi(nodeX(1))) / (FeCard.getPhi(nodeO(1)) - FeCard.getPhi(nodeX(1))); assert( (FeCard.getPhi(nodeO(1)) - FeCard.getPhi(nodeX(1))) != 0.0);
  Real csi12 = (0.5 - FeCard.getPhi(nodeX(1))) / (FeCard.getPhi(nodeO(2)) - FeCard.getPhi(nodeX(1))); assert( (FeCard.getPhi(nodeO(2)) - FeCard.getPhi(nodeX(1))) != 0.0);
  Real csi21 = (0.5 - FeCard.getPhi(nodeX(2))) / (FeCard.getPhi(nodeO(1)) - FeCard.getPhi(nodeX(2))); assert( (FeCard.getPhi(nodeO(1)) - FeCard.getPhi(nodeX(2))) != 0.0);
  Real csi22 = (0.5 - FeCard.getPhi(nodeX(2))) / (FeCard.getPhi(nodeO(2)) - FeCard.getPhi(nodeX(2))); assert( (FeCard.getPhi(nodeO(2)) - FeCard.getPhi(nodeX(2))) != 0.0);
  
  point3d P11 = ElCard.getNode(nodeX(1)) + (ElCard.getNode(nodeO(1)) - ElCard.getNode(nodeX(1))) * csi11;
  point3d P12 = ElCard.getNode(nodeX(1)) + (ElCard.getNode(nodeO(2)) - ElCard.getNode(nodeX(1))) * csi12;
  point3d P21 = ElCard.getNode(nodeX(2)) + (ElCard.getNode(nodeO(1)) - ElCard.getNode(nodeX(2))) * csi21;
  point3d P22 = ElCard.getNode(nodeX(2)) + (ElCard.getNode(nodeO(2)) - ElCard.getNode(nodeX(2))) * csi22;
  
  
  //Build the wedge 1
  sVect<point3d> wedge(6);
  wedge(1) = ElCard.getNode(nodeX(1));
  wedge(2) = P11;
  wedge(3) = P12;
  wedge(4) = ElCard.getNode(nodeX(2));
  wedge(5) = P21;
  wedge(6) = P22;
  
  volLabels.push_back(Xwedge);
  points.push_back(wedge);
  
  
  //Build the wedge 2
  wedge(1) = ElCard.getNode(nodeO(1));
  wedge(2) = P11;
  wedge(3) = P21;
  wedge(4) = ElCard.getNode(nodeO(2));
  wedge(5) = P12;
  wedge(6) = P22;
  
  volLabels.push_back(Xwedge);
  points.push_back(wedge);
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
caseC2(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{
  UInt nodeM = UInt(FeCard.getPhi(1) == 0.5) * 1
             + UInt(FeCard.getPhi(2) == 0.5) * 2
	     + UInt(FeCard.getPhi(3) == 0.5) * 3
	     + UInt(FeCard.getPhi(4) == 0.5) * 4;

  assert(nodeM >= 1); assert(nodeM <= 4);

  UInt nodeX = UInt(FeCard.getPhi(1) > 0.5) * 1
             + UInt(FeCard.getPhi(2) > 0.5) * 2
	     + UInt(FeCard.getPhi(3) > 0.5) * 3
	     + UInt(FeCard.getPhi(4) > 0.5) * 4;
	     
  assert(nodeX >= 1); assert(nodeX <= 4);
	     
  sVect<UInt> nodeO;
  
  if((nodeM != 1) && (nodeX != 1)) {nodeO.push_back(1);}
  if((nodeM != 2) && (nodeX != 2)) {nodeO.push_back(2);}
  if((nodeM != 3) && (nodeX != 3)) {nodeO.push_back(3);}
  if((nodeM != 4) && (nodeX != 4)) {nodeO.push_back(4);}
  
  assert(nodeO.size()  == 2);
  sVect<point3d> tetra(4);
  
  Real csi1 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO(1)) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO(1)) - FeCard.getPhi(nodeX)) != 0.0);
  Real csi2 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO(2)) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO(2)) - FeCard.getPhi(nodeX)) != 0.0);
  
  point3d D1 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO(1)) - ElCard.getNode(nodeX)) * csi1;
  point3d D2 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO(2)) - ElCard.getNode(nodeX)) * csi2;
  
  point3d X1 = ElCard.getNode(nodeX);
  point3d M1 = ElCard.getNode(nodeM);
  point3d O1 = ElCard.getNode(nodeO(1));
  point3d O2 = ElCard.getNode(nodeO(2));
  
  
  //Build the tetra
  if( point3d::dot((X1 - D2) ^ (D1 - D2),M1 - D2) > 0.0 )
  {
    tetra(1) = X1;
    tetra(2) = D1;
    tetra(3) = D2;
    tetra(4) = M1;
  }
  else
  {
    tetra(1) = X1;
    tetra(2) = D2;
    tetra(3) = D1;
    tetra(4) = M1;
  }
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
  
  //Build the tetra
  if( point3d::dot((O1 - D2) ^ (D1 - D2), M1 - D2) > 0.0 )
  {
    tetra(1) = D2;
    tetra(2) = O1;
    tetra(3) = D1;
    tetra(4) = M1;
  }
  else
  {
    tetra(1) = D2;
    tetra(2) = D1;
    tetra(3) = O1;
    tetra(4) = M1;
  }
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
  
  //Build the tetra
  if( point3d::dot((O1 - O2) ^ (D2 - O2), M1 - O1) > 0.0 )
  {
    tetra(1) = D2;
    tetra(2) = O2;
    tetra(3) = O1;
    tetra(4) = M1;
  }
  else
  {
    tetra(1) = D2;
    tetra(2) = O1;
    tetra(3) = O2;
    tetra(4) = M1;
  }
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
caseD1(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{ 
  UInt nodeX = UInt(FeCard.getPhi(1) < 0.5) * 1
             + UInt(FeCard.getPhi(2) < 0.5) * 2
	     + UInt(FeCard.getPhi(3) < 0.5) * 3
	     + UInt(FeCard.getPhi(4) < 0.5) * 4;
	     
  assert(nodeX >= 1); assert(nodeX <= 4);
  
  UInt faceX = GEOSHAPE::faceOppositeToNode(nodeX);
  
  UInt nodeO1 = GEOSHAPE::faceToPoint(faceX,1);
  UInt nodeO2 = GEOSHAPE::faceToPoint(faceX,2);
  UInt nodeO3 = GEOSHAPE::faceToPoint(faceX,3);
  
  Real csi1 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO1) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO1) - FeCard.getPhi(nodeX)) != 0.0);
  Real csi2 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO2) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO2) - FeCard.getPhi(nodeX)) != 0.0);
  Real csi3 = (0.5 - FeCard.getPhi(nodeX)) / (FeCard.getPhi(nodeO3) - FeCard.getPhi(nodeX)); assert( (FeCard.getPhi(nodeO3) - FeCard.getPhi(nodeX)) != 0.0);
  
  point3d P1 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO1) - ElCard.getNode(nodeX)) * csi1;
  point3d P2 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO2) - ElCard.getNode(nodeX)) * csi2;
  point3d P3 = ElCard.getNode(nodeX) + (ElCard.getNode(nodeO3) - ElCard.getNode(nodeX)) * csi3;
  
  //Build the tetra
  sVect<point3d> tetra(4);
  tetra(1) = P1;
  tetra(2) = P3;
  tetra(3) = P2;
  tetra(4) = ElCard.getNode(nodeX);
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
  
  //Build the wedge
  sVect<point3d> wedge(6);
  wedge(1) = ElCard.getNode(nodeO1);
  wedge(2) = ElCard.getNode(nodeO3);
  wedge(3) = ElCard.getNode(nodeO2);
  wedge(4) = P1;
  wedge(5) = P3;
  wedge(6) = P2;
  
  volLabels.push_back(Xwedge);
  points.push_back(wedge);
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
caseD2(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{
  UInt nodeM = UInt(FeCard.getPhi(1) == 0.5) * 1
             + UInt(FeCard.getPhi(2) == 0.5) * 2
	     + UInt(FeCard.getPhi(3) == 0.5) * 3
	     + UInt(FeCard.getPhi(4) == 0.5) * 4;

  assert(nodeM >= 1); assert(nodeM <= 4);

  UInt nodeO = UInt(FeCard.getPhi(1) < 0.5) * 1
             + UInt(FeCard.getPhi(2) < 0.5) * 2
	     + UInt(FeCard.getPhi(3) < 0.5) * 3
	     + UInt(FeCard.getPhi(4) < 0.5) * 4;
	     
  sVect<UInt> nodeX;
  
  if((nodeM != 1) && (nodeO != 1)) {nodeX.push_back(1);}
  if((nodeM != 2) && (nodeO != 2)) {nodeX.push_back(2);}
  if((nodeM != 3) && (nodeO != 3)) {nodeX.push_back(3);}
  if((nodeM != 4) && (nodeO != 4)) {nodeX.push_back(4);}
  
  assert(nodeX.size()  == 2);
  sVect<point3d> tetra(4);
  
  
  Real csi1 = (0.5 - FeCard.getPhi(nodeO)) / (FeCard.getPhi(nodeX(1)) - FeCard.getPhi(nodeO)); assert( (FeCard.getPhi(nodeX(1)) - FeCard.getPhi(nodeO)) != 0.0);
  Real csi2 = (0.5 - FeCard.getPhi(nodeO)) / (FeCard.getPhi(nodeX(2)) - FeCard.getPhi(nodeO)); assert( (FeCard.getPhi(nodeX(2)) - FeCard.getPhi(nodeO)) != 0.0);
  
  point3d D1 = ElCard.getNode(nodeO) + (ElCard.getNode(nodeX(1)) - ElCard.getNode(nodeO)) * csi1;
  point3d D2 = ElCard.getNode(nodeO) + (ElCard.getNode(nodeX(2)) - ElCard.getNode(nodeO)) * csi2;
  
  point3d O1 = ElCard.getNode(nodeO);
  point3d M1 = ElCard.getNode(nodeM);
  point3d X1 = ElCard.getNode(nodeX(1));
  point3d X2 = ElCard.getNode(nodeX(2));
  
  //Build the tetra
  if( point3d::dot((O1 - D2) ^ (D1 - D2), M1 - D2) > 0.0 )
  {
    tetra(1) = O1;
    tetra(2) = D1;
    tetra(3) = D2;
    tetra(4) = M1;
  }
  else
  {
    tetra(1) = O1;
    tetra(2) = D2;
    tetra(3) = D1;
    tetra(4) = M1;
  }
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
  
  //Build the tetra
  if( point3d::dot((X1 - D2) ^ (D1 - D2), M1 - D2) > 0.0 )
  {
    tetra(1) = D2;
    tetra(2) = X1;
    tetra(3) = D1;
    tetra(4) = M1;
  }
  else
  {
    tetra(1) = D2;
    tetra(2) = D1;
    tetra(3) = X1;
    tetra(4) = M1;
  }
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
  
  //Build the tetra
  if( point3d::dot((X1 - X2) ^ (D2 - X2), M1 - X1) > 0.0 )
  {
    tetra(1) = D2;
    tetra(2) = X2;
    tetra(3) = X1;
    tetra(4) = M1;
  }
  else
  {
    tetra(1) = D2;
    tetra(2) = X1;
    tetra(3) = X2;
    tetra(4) = M1;
  }
  
  volLabels.push_back(Xtetra);
  points.push_back(tetra);
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
caseD3(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{
  UInt nodeX = UInt(FeCard.getPhi(1) > 0.5) * 1
             + UInt(FeCard.getPhi(2) > 0.5) * 2
	     + UInt(FeCard.getPhi(3) > 0.5) * 3
	     + UInt(FeCard.getPhi(4) > 0.5) * 4;
	     
  assert(nodeX >= 1); assert(nodeX <= 4);
  
  UInt nodeO = UInt(FeCard.getPhi(1) < 0.5) * 1
             + UInt(FeCard.getPhi(2) < 0.5) * 2
	     + UInt(FeCard.getPhi(3) < 0.5) * 3
	     + UInt(FeCard.getPhi(4) < 0.5) * 4;
	     
  assert(nodeO >= 1); assert(nodeO <= 4);
  
  sVect<UInt> nodeM;
  
  if((nodeO != 1) && (nodeX != 1)) {nodeM.push_back(1);}
  if((nodeO != 2) && (nodeX != 2)) {nodeM.push_back(2);}
  if((nodeO != 3) && (nodeX != 3)) {nodeM.push_back(3);}
  if((nodeO != 4) && (nodeX != 4)) {nodeM.push_back(4);}
  
  assert(nodeM.size()  == 2);
  sVect<point3d> tetra(4);
  
  Real csi  = (0.5 - FeCard.getPhi(nodeO)) / (FeCard.getPhi(nodeX) - FeCard.getPhi(nodeO)); assert( (FeCard.getPhi(nodeX) - FeCard.getPhi(nodeO)) != 0.0);
  point3d D = ElCard.getNode(nodeO) + (ElCard.getNode(nodeX) - ElCard.getNode(nodeO)) * csi;
  
  point3d O1 = ElCard.getNode(nodeO);
  point3d X1 = ElCard.getNode(nodeX);
  point3d M1 = ElCard.getNode(nodeM(1));
  point3d M2 = ElCard.getNode(nodeM(2));
  
  //Build the tetra
  if( point3d::dot((M1 - D) ^ (M2 - D), O1 - D) > 0.0 )
  {
    tetra(1) = D;
    tetra(2) = M1;
    tetra(3) = M2;
    tetra(4) = O1;
    
    volLabels.push_back(Xtetra);
    points.push_back(tetra);
    
    tetra(1) = D;
    tetra(2) = M2;
    tetra(3) = M1;
    tetra(4) = X1;
    
    volLabels.push_back(Xtetra);
    points.push_back(tetra);
  }
  else
  {
    tetra(1) = D;
    tetra(2) = M2;
    tetra(3) = M1;
    tetra(4) = O1;
    
    volLabels.push_back(Xtetra);
    points.push_back(tetra);
    
    tetra(1) = D;
    tetra(2) = M1;
    tetra(3) = M2;
    tetra(4) = X1;
    
    volLabels.push_back(Xtetra);
    points.push_back(tetra);
  }
}

template<typename PMAPTYPE>
void
xfemSupportLT3d<PMAPTYPE>::
split(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & points) const
{
  //Clearing variables
  volLabels.clear();
  points.clear();
  
  //Counting
  UInt valAbove = UInt(FeCard.getPhi(1) >= 0.5)
                + UInt(FeCard.getPhi(2) >= 0.5)
	        + UInt(FeCard.getPhi(3) >= 0.5)
	        + UInt(FeCard.getPhi(4) >= 0.5);

  UInt valEqual = UInt(FeCard.getPhi(1) == 0.5)
                + UInt(FeCard.getPhi(2) == 0.5)
	        + UInt(FeCard.getPhi(3) == 0.5)
	        + UInt(FeCard.getPhi(4) == 0.5);

  switch(valAbove)
  {
    case 1 :
      caseB(FeCard,ElCard,volLabels,points);
      break;
      
    case 2 :
      if(valEqual == 0) {caseC1(FeCard,ElCard,volLabels,points);}
      if(valEqual == 1) {caseC2(FeCard,ElCard,volLabels,points);}
      break;
      
    case 3 :
      if(valEqual == 0) {caseD1(FeCard,ElCard,volLabels,points);}
      if(valEqual == 1) {caseD2(FeCard,ElCard,volLabels,points);}
      if(valEqual == 2) {caseD3(FeCard,ElCard,volLabels,points);}
      break;
  }
}

template<typename PMAPTYPE>
Real
xfemSupportLT3d<PMAPTYPE>::
tetraVolume(const sVect<point3d> & points) const
{
  assert(points.size() == 4);
  
  geoMapSupport3d<linearTetra> interface;
  interface.setPoints(points);
  
  Real vol = interface.volume<STANDARD,1>();
  
  assert(vol >= 0.0);
  return(vol);
}

template<typename PMAPTYPE>
Real
xfemSupportLT3d<PMAPTYPE>::
wedgeVolume(const sVect<point3d> & points) const
{
  assert(points.size() == 6);
  
  typedef linearWedge WEDGE;
  typedef intStaticCard<WEDGE,STANDARD,2> INTTAB;
  
  Real vol = 0.0;
  geoMapInterface<WEDGE> interface;
  
  for(UInt i=1; i <= INTTAB::N; ++i)
  {
    vol += interface.getGradientDet(points,INTTAB::getYn(i)) * INTTAB::getWn(i);
  }
  
  assert(vol >= 0.0);
  return(vol);
}

template<typename PMAPTYPE>
sVect<Real>
xfemSupportLT3d<PMAPTYPE>::
volumes(const FECARD & FeCard, const ELCARD & ElCard) const
{
  sVect<Real>             volumes;
  sVect<Xvolumes>         volLabels;
  sVect<sVect<point3d> >  subPoints;
  
  split(FeCard,ElCard,volLabels,subPoints);
  
  if(volLabels.size() == 0)
  {
    volumes.push_back(tetraVolume(ElCard.getNodes()));
  }
  else
  {
    for(UInt i=1; i<=volLabels.size(); ++i)
    {
      if(volLabels(i) == Xtetra) { volumes.push_back( tetraVolume(subPoints(i)) ); }
      else                       { volumes.push_back( wedgeVolume(subPoints(i)) ); }
    }
  }
  
  return(volumes);
}

template<typename PMAPTYPE>
sVect<Real>
xfemSupportLT3d<PMAPTYPE>::
volumes(const FECARD & FeCard, const ELCARD & ElCard, sVect<Xvolumes> & volLabels, sVect< sVect<point3d> > & subPoints) const
{
  volLabels.clear();
  subPoints.clear();
  sVect<Real> volumes;
  
  split(FeCard,ElCard,volLabels,subPoints);
  
  if(volLabels.size() == 0)
  {
    volumes.push_back(tetraVolume(ElCard.getNodes()));
  }
  else
  {
    for(UInt i=1; i<=volLabels.size(); ++i)
    {
      if(volLabels(i) == Xtetra) { volumes.push_back( tetraVolume(subPoints(i)) ); }
      else                       { volumes.push_back( wedgeVolume(subPoints(i)) ); }
    }
  }
  
  return(volumes);
}

template<typename PMAPTYPE>
bool
xfemSupportLT3d<PMAPTYPE>::
isSplit(const FECARD & FeCard, const ELCARD & ElCard) const
{
  bool easyCheck = 
  ( (FeCard.getPhi(1) > 0.5) && (FeCard.getPhi(2) > 0.5) && (FeCard.getPhi(3) > 0.5) && (FeCard.getPhi(4) > 0.5) ) || 
  ( (FeCard.getPhi(1) < 0.5) && (FeCard.getPhi(2) < 0.5) && (FeCard.getPhi(3) < 0.5) && (FeCard.getPhi(4) < 0.5) );
  
  if(easyCheck)
  {
    return(false);
  }
  else
  {
    sVect<Real> vols = volumes(FeCard,ElCard);
    
    if(vols.size() == 1) {return(false);}
    else
    {
      return(( std::min(vols(1),vols(2)) / std::max(vols(1),vols(2)) ) >= FeCard.getXtoll() );
    }
  }
}

template<typename PMAPTYPE>
sVect<point3d>
xfemSupportLT3d<PMAPTYPE>::
getBarycenters(const FECARD & FeCard, const ELCARD & ElCard) const
{
  sVect<Xvolumes>         volLabels;
  sVect< sVect<point3d> > points;
  sVect<point3d>          out;
  point3d B;
  
  split(FeCard,ElCard,volLabels,points);
  
  for(UInt i=1; i <= points.size(); ++i)
  {
    B.setX(0.0); B.setY(0.0); B.setZ(0.0);
    
    for(UInt j=1; j <= points(i).size(); ++j)
    {
      B += points(i)(j);
    }
    
    out.push_back(B / points(i).size());
  }
  
  return(out);
}


#endif

