/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FERT0LOC3D_HPP
#define FERT0LOC3D_HPP

#include "morganaTypes.hpp"

#include "geoMapInterface.hpp"

#include "morganaFields.hpp"
#include "elCard3d.hpp"
#include "feStaticDofCard3d.h"
#include "feStaticEvalIterators.hpp"

#include "feMapPiola3d.hpp"


/*! Local raviart Thomas - card */
class feRt0Loc3d_card
{ };


/*! Local Raviart-Thomas 3d zero order */
template<typename PMAPTYPE>
class feRt0Loc3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef point3d                     BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feRt0Loc3d_card             FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = vectorial;
    static const FELabel     feLabel     = FE_RT0LOC_3d;
    static const FEBaseLabel feBaseLabel = BL_RT0LOC_3d;    
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 4;
    static const UInt numBasis = 4;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 0;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 0;
    static const UInt dofPerVolume = 4;
    static const bool isNodal = false;
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
    feRt0Loc3d();
    feRt0Loc3d(const FECARD & FECard, const ELCARD & ELCard);
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
feRt0Loc3d<PMAPTYPE>::
feRt0Loc3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feRt0Loc3d<PMAPTYPE>::
feRt0Loc3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feRt0Loc3d<PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feRt0Loc3d<PMAPTYPE>::DOFCARD
feRt0Loc3d<PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  static DOFCARD list[4] =
  {
    DOFCARD(VOLUME,1,1),
    DOFCARD(VOLUME,2,1),
    DOFCARD(VOLUME,3,1),
    DOFCARD(VOLUME,4,1) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feRt0Loc3d<PMAPTYPE>::
getRefNode(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  static point3d list[4] =
  {
    point3d(1.0,1.0,0.0) / 3.0,
    point3d(1.0,0.0,1.0) / 3.0,
    point3d(1.0,1.0,1.0) / 3.0,
    point3d(0.0,1.0,1.0) / 3.0 
  };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
void
feRt0Loc3d<PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  feStaticPoint3dEval<feLabel,numPoly>::eval(val,Y);
}

template<typename PMAPTYPE>
void
feRt0Loc3d<PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  //Polynomial evaluation
  feStaticPoint3dDx<feLabel,numPoly>::eval(gradX,Y);
  feStaticPoint3dDy<feLabel,numPoly>::eval(gradY,Y);
  feStaticPoint3dDz<feLabel,numPoly>::eval(gradZ,Y);
}

template<typename PMAPTYPE>
void
feRt0Loc3d<PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  //Polynomial evaluation
  localEval(Y,val);
  
  //Mapping
  feMapPiola3d<GEOSHAPE,PMAPTYPE> piolaMap(ElCard);
  val(1) = piolaMap.eval(val(1),Y);
  val(2) = piolaMap.eval(val(2),Y);
  val(3) = piolaMap.eval(val(3),Y);
  val(4) = piolaMap.eval(val(4),Y);
}

template<typename PMAPTYPE>
void
feRt0Loc3d<PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  //Asserts
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numBasis);
  assert(gradY.size() == numBasis);
  assert(gradZ.size() == numBasis);
  
  //Alloc
  sVect<BASETYPE> val(numPoly);
  sVect<BASETYPE> locX(numBasis), locY(numBasis), locZ(numBasis);
  
  //Local gradient 
  localEval(Y,val);
  localGrad(Y,locX,locY,locZ);
  
  //Gradient mapping
  feMapPiola3d<GEOSHAPE,PMAPTYPE> piolaMap(ElCard);
  
  
  //Map -1 or 1 
  for(UInt i=1; i <= 4; ++i)
  {
    locX(i) = piolaMap.evalGradX(val(i), locX(i), Y);
    locY(i) = piolaMap.evalGradY(val(i), locY(i), Y);
    locZ(i) = piolaMap.evalGradZ(val(i), locZ(i), Y);
  }
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.computeInverse();
  
  //Compute gradient
  for(UInt i=1; i <= numBasis; ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}



#endif
