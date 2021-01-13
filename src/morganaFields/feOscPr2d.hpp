/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FEOSCPR2D_HPP
#define FEOSCPR2D_HPP

#include "geoShapes.h"
#include "elCard2d.hpp"
#include "morganaFields.hpp"
#include "feOsc3d_card.h"
#include "feStaticDofCard2d.h"
#include "feStaticEvalIterators.hpp"
#include "../morganaDofs/komplex.h"


//Forward declaration______________________________________________________________________________
template<Int R, Int NW, typename PMAPTYPE>
class feOscPr2d
{
};


//_________________________________________________________________________________________________
//   P0/OSC0   P0/OSC0   P0/OSC0   P0/OSC0   P0/OSC0   P0/OSC0   P0/OSC0   P0/OSC0
//-------------------------------------------------------------------------------------------------

/*! Tetrahedral nodal finite element - P0 with oscillating functions */
template<typename PMAPTYPE>
class feOscPr2d<0,0,PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef komplex                     BASETYPE;
    typedef linearTriangle              GEOSHAPE;
    typedef feOsc3d_card                FECARD;
    typedef elCard2d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard2d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = komplexx;
    static const FELabel     feLabel     = FE_OS0P0_2d;
    static const FEBaseLabel feBaseLabel = BL_OSwPr_2d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 1;
    static const UInt numBasis = 1;
    static const UInt dim      = 2;
    static const UInt dofPerVertex = 0;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 1;
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
    feOscPr2d();
    feOscPr2d(const FECARD & FECard, const ELCARD & ELCard);
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
feOscPr2d<0,0,PMAPTYPE>::
feOscPr2d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feOscPr2d<0,0,PMAPTYPE>::
feOscPr2d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feOscPr2d<0,0,PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feOscPr2d<0,0,PMAPTYPE>::DOFCARD
feOscPr2d<0,0,PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i == 1);  
  return( DOFCARD(VOLUME,1,1) );
}

template<typename PMAPTYPE>
point3d
feOscPr2d<0,0,PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 1);
  
  static point3d list[1] =
  {
    point3d(1.0/3.0,  1.0/3.0,  0.0)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feOscPr2d<0,0,PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation P1
  val(1) = komplex(1.0,0.0);
}

template<typename PMAPTYPE>
void
feOscPr2d<0,0,PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  //Polynomial evaluation
  gradX(1) = komplex(0.0,0.0);
  gradY(1) = komplex(0.0,0.0);
  gradZ(1) = komplex(0.0,0.0);
}

template<typename PMAPTYPE>
void
feOscPr2d<0,0,PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  val(1) = komplex(1.0,0.0);
}

template<typename PMAPTYPE>
void
feOscPr2d<0,0,PMAPTYPE>::
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
  T.completeThirdColoumn();
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
//   P1/OSC0   P1/OSC0   P1/OSC0   P1/OSC0   P1/OSC0   P1/OSC0   P1/OSC0   P1/OSC0
//-------------------------------------------------------------------------------------------------

/*! Triangular nodal finite element - P1 */
template<typename PMAPTYPE>
class feOscPr2d<1,0,PMAPTYPE>
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef komplex                     BASETYPE;
    typedef linearTriangle              GEOSHAPE;
    typedef feOsc3d_card                FECARD;
    typedef elCard2d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard2d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = komplexx;
    static const FELabel     feLabel     = FE_OS0P0_2d;
    static const FEBaseLabel feBaseLabel = BL_OSwPr_2d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 3;
    static const UInt numBasis = 3;
    static const UInt dim      = 2;
    static const UInt dofPerVertex = 1;
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
    feOscPr2d();
    feOscPr2d(const FECARD & FECard, const ELCARD & ELCard);
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
feOscPr2d<1,0,PMAPTYPE>::
feOscPr2d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feOscPr2d<1,0,PMAPTYPE>::
feOscPr2d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feOscPr2d<1,0,PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feOscPr2d<1,0,PMAPTYPE>::DOFCARD
feOscPr2d<1,0,PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 3);
  
  static DOFCARD list[3] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,1,3) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feOscPr2d<1,0,PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 3);
  
  static point3d list[3] =
  {
    point3d(0.0,  0.0,  0.0),
    point3d(1.0,  0.0,  0.0),
    point3d(0.0,  1.0,  0.0)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feOscPr2d<1,0,PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation P1
  sVect<Real> val_P1(3);
  feStaticRealEval<FE_P1_2d,3>::eval(val_P1,Y); 
  
  //Compute basis  
  for(UInt i=1; i <= 3; ++i)
  { val(i) = komplex(val_P1(i));}
}

template<typename PMAPTYPE>
void
feOscPr2d<1,0,PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  //Polynomial evaluation P1
  sVect<Real> gradX_P1(3);
  sVect<Real> gradY_P1(3);
  sVect<Real> gradZ_P1(3);
  
  feStaticRealDx<FE_P1_2d,3>::eval(gradX_P1,Y);
  feStaticRealDy<FE_P1_2d,3>::eval(gradY_P1,Y);
  feStaticRealDz<FE_P1_2d,3>::eval(gradZ_P1,Y);
  
  //Compute basis
  for(UInt i=1; i <= 3; ++i)
  {
    gradX(i) = komplex(gradX_P1(i),0.0);
    gradY(i) = komplex(gradY_P1(i),0.0);
    gradZ(i) = komplex(gradZ_P1(i),0.0);
  }
}

template<typename PMAPTYPE>
void
feOscPr2d<1,0,PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation P1
  sVect<Real> val_P1(3);
  feStaticRealEval<FE_P1_2d,3>::eval(val_P1,Y); 
  
  //Compute basis  
  for(UInt i=1; i <= 3; ++i)
  { val(i) = komplex(val_P1(i),0.0);}
}

template<typename PMAPTYPE>
void
feOscPr2d<1,0,PMAPTYPE>::
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
  T.completeThirdColoumn();
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
