/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FED0LT2D_HPP
#define FED0LT2D_HPP

#include "morganaTypes.hpp"

#include "geoMapInterface.hpp"

#include "morganaFields.hpp"
#include "elCard2d.hpp"
#include "feStaticDofCard2d.h"
#include "feStaticEvalIterators.hpp"



//_________________________________________________________________________________________________
// FE CARD
//-------------------------------------------------------------------------------------------------

/*! Dummy card */
class feD0LT2d_card
{
  public:
    feD0LT2d_card() { };
};



//_________________________________________________________________________________________________
//   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0   D0  
//-------------------------------------------------------------------------------------------------

/*! Dual grid finite element 2d */
template<typename PMAPTYPE>
class feD0LT2d
{
  /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTriangle              GEOSHAPE;
    typedef feD0LT2d_card               FECARD;
    typedef elCard2d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard2d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = dual;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_D0LT_2d;
    static const FEBaseLabel feBaseLabel = BL_D0LT_2d;
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
    feD0LT2d();
    feD0LT2d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    bool baseIsActive(const UInt & i, const point3d & Y) const;
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
feD0LT2d<PMAPTYPE>::
feD0LT2d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feD0LT2d<PMAPTYPE>::
feD0LT2d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feD0LT2d<PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
bool
feD0LT2d<PMAPTYPE>::
baseIsActive(const UInt & i, const point3d & Y) const
{
  assert(i >= 1);
  assert(i <= 3);
  
  Real dist[3] = {
    point3d::norm2( getRefNode(1) - Y ),
    point3d::norm2( getRefNode(2) - Y ),
    point3d::norm2( getRefNode(3) - Y )
  };
  
  Real minDist = min( min(dist[0],dist[1]), dist[2] );
  
  return(minDist == dist[i-1]);
}

template<typename PMAPTYPE>
typename feD0LT2d<PMAPTYPE>::DOFCARD
feD0LT2d<PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 3);
  
  static DOFCARD list[3] =
  {
    DOFCARD(VERTEX,1,1),
    DOFCARD(VERTEX,1,2),
    DOFCARD(VERTEX,1,3)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feD0LT2d<PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 3);
  
  static point3d list[3] =
  {
    point3d(0.0, 0.0, 0.0),
    point3d(1.0, 0.0, 0.0),
    point3d(0.0, 1.0, 0.0)
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feD0LT2d<PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  sVect<Real> polyVal(numPoly);
  feStaticRealEval<feLabel,numPoly>::eval(polyVal,Y);
  
  //Active terms
  val(1) = polyVal(1) * baseIsActive(1,Y);
  val(2) = polyVal(2) * baseIsActive(2,Y);
  val(3) = polyVal(3) * baseIsActive(3,Y);
}

template<typename PMAPTYPE>
void
feD0LT2d<PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Polynomial evaluation
  sVect<Real> polyGradX(numPoly);
  sVect<Real> polyGradY(numPoly);
  sVect<Real> polyGradZ(numPoly);
  
  feStaticRealDx<feLabel,numPoly>::eval(polyGradX,Y);
  feStaticRealDy<feLabel,numPoly>::eval(polyGradY,Y);
  feStaticRealDz<feLabel,numPoly>::eval(polyGradZ,Y);
  
  //Active terms X
  gradX(1) = polyGradX(1) * baseIsActive(1,Y);
  gradX(2) = polyGradX(2) * baseIsActive(2,Y);
  gradX(3) = polyGradX(3) * baseIsActive(3,Y);
	     
  //Active terms Y	     
  gradY(1) = polyGradY(1) * baseIsActive(1,Y);
  gradY(2) = polyGradY(2) * baseIsActive(2,Y);
  gradY(3) = polyGradY(3) * baseIsActive(3,Y);
	          
  //Active terms Y	     
  gradZ(1) = polyGradZ(1) * baseIsActive(1,Y);
  gradZ(2) = polyGradZ(2) * baseIsActive(2,Y);
  gradZ(3) = polyGradZ(3) * baseIsActive(3,Y);
}

template<typename PMAPTYPE>
void
feD0LT2d<PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numBasis);
  
  //Polynomial evaluation
  localEval(Y,val);
}

template<typename PMAPTYPE>
void
feD0LT2d<PMAPTYPE>::
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
