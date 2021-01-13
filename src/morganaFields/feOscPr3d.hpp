/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEOSCPR3D_HPP
#define FEOSCPR3D_HPP

#include "geoShapes.h"
#include "morganaFields.hpp"
#include "feOsc3d_card.h"
#include "feStaticDofCard3d.h"
#include "feStaticEvalIterators.hpp"
#include "../morganaDofs/komplex.h"


//Forward declaration______________________________________________________________________________
template<Int R, Int NW, typename PMAPTYPE>
class feOscPr3d
{
};


//_________________________________________________________________________________________________
//   P0/OSC1   P0/OSC1   P0/OSC1   P0/OSC1   P0/OSC1   P0/OSC1   P0/OSC1   P0/OSC1   P0/OSC1
//-------------------------------------------------------------------------------------------------

/*! Tetrahedral nodal finite element - P0 with oscillating functions */
template<typename PMAPTYPE>
class feOscPr3d<0,1,PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef komplex                     BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feOsc3d_card                FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = komplexx;
    static const FELabel     feLabel     = FE_OS1P0_3d;
    static const FEBaseLabel feBaseLabel = BL_OSwPr_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 6;
    static const UInt numBasis = 6;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 0;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 6;
    static const bool isNodal = true;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool elCardLoaded;
    bool feCardLoaded;
    FECARD FeCard;
    ELCARD ElCard;
    //@}
    
    /*! @name Constructors and setting functions */ //@{
  public:
    feOscPr3d();
    feOscPr3d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    //@}
    
    /*! @name Eval functions */ //@{
  public:
    static DOFCARD getDofCard(const UInt & i);
    static point3d getRefNode(const UInt & i);
    void localEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    void globalEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    //@}
};

template<typename PMAPTYPE>
feOscPr3d<0,1,PMAPTYPE>::
feOscPr3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feOscPr3d<0,1,PMAPTYPE>::
feOscPr3d(const FECARD & FECard, const ELCARD & ELCard)
{
  assert(FECard.size() == 6);
  
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feOscPr3d<0,1,PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  assert(FECard.size() == 6);
  
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feOscPr3d<0,1,PMAPTYPE>::DOFCARD
feOscPr3d<0,1,PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 6);
  
  static DOFCARD list[6] =
  {
    DOFCARD(VOLUME,1,1),
    DOFCARD(VOLUME,2,1),
    DOFCARD(VOLUME,3,1),
    DOFCARD(VOLUME,4,1),
    DOFCARD(VOLUME,5,1),
    DOFCARD(VOLUME,6,1)
  };
  
  return(list[i-1]);
}

template<typename PMAPTYPE>
point3d
feOscPr3d<0,1,PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 6);
  
  return(point3d(0.25, 0.25, 0.25));
}

template<typename PMAPTYPE>
void
feOscPr3d<0,1,PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  Real theta;
  
  for(UInt i=1; i <= 6; ++i)
  {
    theta  = point3d::dot(FeCard.getH(i), Y - FeCard.getY(i));
    val(i) = komplex(cos(theta), sin(theta));
  }
}

template<typename PMAPTYPE>
void
feOscPr3d<0,1,PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  Real theta;
  komplex C;
  
  for(UInt i=1; i <= 6; ++i)
  {
    theta = point3d::dot(FeCard.getH(i), Y - FeCard.getY(i));
    C     = komplex(-sin(theta),cos(theta));
    
    gradX(i) = C * FeCard.getH(i).getX();
    gradY(i) = C * FeCard.getH(i).getY();
    gradZ(i) = C * FeCard.getH(i).getZ();
  }
}

template<typename PMAPTYPE>
void
feOscPr3d<0,1,PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  localEval(Y,val);
}

template<typename PMAPTYPE>
void
feOscPr3d<0,1,PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.computeInverse();
  
  //Local grad
  sVect<BASETYPE> locX(numBasis), locY(numBasis), locZ(numBasis);
  localGrad(Y,locX,locY,locZ);
  
  //Compute gradient
  for(UInt i=1; i <= numBasis; ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}


//_________________________________________________________________________________________________
//   P1/OSC1   P1/OSC1   P1/OSC1   P1/OSC1   P1/OSC1   P1/OSC1   P1/OSC1   P1/OSC1   P1/OSC1
//-------------------------------------------------------------------------------------------------

/*! Tetrahedral nodal finite element - P1 with oscillating functions */
template<typename PMAPTYPE>
class feOscPr3d<1,1,PMAPTYPE>
{
  /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef komplex                     BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feOsc3d_card                FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = komplexx;
    static const FELabel     feLabel     = FE_OS1P1_3d;
    static const FEBaseLabel feBaseLabel = BL_OSwPr_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 24;
    static const UInt numBasis = 24;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 6;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 0;
    static const bool isNodal = true;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool elCardLoaded;
    bool feCardLoaded;
    FECARD FeCard;
    ELCARD ElCard;
    //@}
    
    /*! @name Constructors and setting functions */ //@{
  public:
    feOscPr3d();
    feOscPr3d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    //@}
    
    /*! @name Eval functions */ //@{
  public:
    static DOFCARD getDofCard(const UInt & i);
    static point3d getRefNode(const UInt & i);
    void localEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    void globalEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    //@}
};


template<typename PMAPTYPE>
feOscPr3d<1,1,PMAPTYPE>::
feOscPr3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feOscPr3d<1,1,PMAPTYPE>::
feOscPr3d(const FECARD & FECard, const ELCARD & ELCard)
{
  assert(FECard.size() == 24);
  
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feOscPr3d<1,1,PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  assert(FECard.size() == 24);
  
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feOscPr3d<1,1,PMAPTYPE>::DOFCARD
feOscPr3d<1,1,PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 24);
  
  static DOFCARD list[24] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,2,1),
    DOFCARD(VERTEX,3,1),
    DOFCARD(VERTEX,4,1),
    DOFCARD(VERTEX,5,1),
    DOFCARD(VERTEX,6,1),
    
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,2,2),
    DOFCARD(VERTEX,3,2),
    DOFCARD(VERTEX,4,2),
    DOFCARD(VERTEX,5,2),
    DOFCARD(VERTEX,6,2),
    
    DOFCARD(VERTEX,1,3),
    DOFCARD(VERTEX,2,3),
    DOFCARD(VERTEX,3,3),
    DOFCARD(VERTEX,4,3),
    DOFCARD(VERTEX,5,3),
    DOFCARD(VERTEX,6,3),
    
    DOFCARD(VERTEX,1,4),
    DOFCARD(VERTEX,2,4),
    DOFCARD(VERTEX,3,4),
    DOFCARD(VERTEX,4,4),
    DOFCARD(VERTEX,5,4),
    DOFCARD(VERTEX,6,4)
  };
  
  return(list[i-1]);
}

template<typename PMAPTYPE>
point3d
feOscPr3d<1,1,PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >=  1);
  assert(i <= 24);
  
  point3d Y = point3d(0.0, 0.0, 0.0) * Real(i >= 1)  * Real(i <= 6)  +
              point3d(1.0, 0.0, 0.0) * Real(i >= 7)  * Real(i <= 12) +
              point3d(0.0, 1.0, 0.0) * Real(i >= 13) * Real(i <= 18) +
              point3d(0.0, 0.0, 1.0) * Real(i >= 19) * Real(i <= 24);
  
  return(Y);
}

template<typename PMAPTYPE>
void
feOscPr3d<1,1,PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation P1
  sVect<Real> val_P1(4);
  feStaticRealEval<FE_P1_3d,4>::eval(val_P1,Y); 
  
  //Compute basis 
  Real theta;
  UInt k=1;
  
  for(UInt i=1; i <= 4; ++i)
  {    
    for(UInt j=1; j <= 6; ++j)
    {
      theta  = point3d::dot(FeCard.getH(k), Y - FeCard.getY(k));
      val(k) = komplex(cos(theta), sin(theta)) * val_P1(i);
      
      k++;
    }
  }
}

template<typename PMAPTYPE>
void
feOscPr3d<1,1,PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  //Polynomial evaluation P1
  sVect<Real> val_P1(4);
  sVect<Real> gradX_P1(4);
  sVect<Real> gradY_P1(4);
  sVect<Real> gradZ_P1(4);
  
  feStaticRealEval<FE_P1_3d,4>::eval(val_P1,Y);
  feStaticRealDx<FE_P1_3d,4>::eval(gradX_P1,Y);
  feStaticRealDy<FE_P1_3d,4>::eval(gradY_P1,Y);
  feStaticRealDz<FE_P1_3d,4>::eval(gradZ_P1,Y);
  
  //Compute basis
  komplex C, G;
  Real theta;
  UInt k=1;
  
  for(UInt i=1; i <= 4; ++i)
  {    
    for(UInt j=1; j <= 6; ++j)
    {
      theta = point3d::dot(FeCard.getH(k), Y - FeCard.getY(k));
      G     = komplex(-sin(theta), cos(theta));
      C     = komplex( cos(theta), sin(theta));
      
      gradX(k) = G * FeCard.getH(j).getX() * val_P1(i) + C * gradX_P1(i);
      gradY(k) = G * FeCard.getH(j).getY() * val_P1(i) + C * gradY_P1(i);
      gradZ(k) = G * FeCard.getH(j).getZ() * val_P1(i) + C * gradZ_P1(i);
      
      k++;
    }
  }
}

template<typename PMAPTYPE>
void
feOscPr3d<1,1,PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  localEval(Y,val);
}

template<typename PMAPTYPE>
void
feOscPr3d<1,1,PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.computeInverse();
  
  //Local grad
  sVect<BASETYPE> locX(numBasis), locY(numBasis), locZ(numBasis);
  localGrad(Y,locX,locY,locZ);
  
  //Compute gradient
  for(UInt i=1; i <= numBasis; ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}

#endif
