/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEFVFLUX3D_HPP
#define FEFVFLUX3D_HPP


#include "morganaTypes.hpp"
#include "geoMapInterface.hpp"

#include "morganaFields.hpp"
#include "elCard3d.hpp"
#include "feStaticDofCard3d.h"
#include "feStaticEvalIterators.hpp"


//_________________________________________________________________________________________________
// FE CARD
//-------------------------------------------------------------------------------------------------
class feFvFlux3d_card
{
  public:
    feFvFlux3d_card() { };
};



//_________________________________________________________________________________________________
//  FE-FLUX
//-------------------------------------------------------------------------------------------------

/*! Flux finite element for Finite Volume applications */
template<typename PMAPTYPE>
class feFvFlux3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feFvFlux3d_card             FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feStaticDofCard3d           DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feStatic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_FVFLUX_3d;
    static const FEBaseLabel feBaseLabel = BL_FVFLUX_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt numPoly  = 4;
    static const UInt numBasis = 4;
    static const UInt dim      = 3;
    static const UInt dofPerVertex = 0;
    static const UInt dofPerEdge   = 0;
    static const UInt dofPerFace   = 1;
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
    feFvFlux3d();
    feFvFlux3d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    //@}
    
    /*! @name Info functions */ //@{
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
feFvFlux3d<PMAPTYPE>::
feFvFlux3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feFvFlux3d<PMAPTYPE>::
feFvFlux3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feFvFlux3d<PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
typename feFvFlux3d<PMAPTYPE>::DOFCARD
feFvFlux3d<PMAPTYPE>::
getDofCard(const UInt & i)
{
  assert(i >= 1);
  assert(i <= 4);
  
  static DOFCARD list[4] =
  {
    DOFCARD(FACE,1,1),
    DOFCARD(FACE,1,2),
    DOFCARD(FACE,1,3),
    DOFCARD(FACE,1,4) };
  
  return( list[i-1] );
}

template<typename PMAPTYPE>
point3d
feFvFlux3d<PMAPTYPE>::
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
feFvFlux3d<PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  /*
  static const Real relGeoToll = 1.0e-3;
    
  val(1) = Real( (Y.getZ() <= relGeoToll)    &&  (Y.getZ() >= (-relGeoToll))      &&
                 (Y.getY() >= (-relGeoToll)) &&  (Y.getY() <= (1.0 + relGeoToll)) &&
                 (Y.getX() >= (-relGeoToll)) &&  (Y.getX() <= (1.0 + relGeoToll)) &&
                 ( (Y.getX() + Y.getY()) <= (1.0 + relGeoToll) ) );
   
  val(2) = Real( (Y.getY() <= relGeoToll)    &&  (Y.getY() >= (-relGeoToll))      &&
                 (Y.getZ() >= (-relGeoToll)) &&  (Y.getZ() <= (1.0 + relGeoToll)) &&
                 (Y.getX() >= (-relGeoToll)) &&  (Y.getX() <= (1.0 + relGeoToll)) &&
                 ( (Y.getX() + Y.getZ()) <= (1.0 + relGeoToll) ) );
  
  val(3) = Real( (Y.getX() >= (-relGeoToll)) &&  (Y.getX() <= (1.0 + relGeoToll)) &&
                 (Y.getY() >= (-relGeoToll)) &&  (Y.getY() <= (1.0 + relGeoToll)) &&
                 (Y.getZ() >= (-relGeoToll)) &&  (Y.getZ() <= (1.0 + relGeoToll)) &&
                 ( (Y.getX() + Y.getY() + Y.getZ()) <= (1.0 + relGeoToll) )    &&
                 ( (Y.getX() + Y.getY() + Y.getZ()) >= (1.0 - relGeoToll) ) );
   
  val(4) = Real( (Y.getX() <= relGeoToll)    &&  (Y.getX() >= (-relGeoToll))      &&
                 (Y.getZ() >= (-relGeoToll)) &&  (Y.getZ() <= (1.0 + relGeoToll)) &&
                 (Y.getY() >= (-relGeoToll)) &&  (Y.getY() <= (1.0 + relGeoToll)) &&
                 ( (Y.getY() + Y.getZ()) <= (1.0 + relGeoToll) ) );
  */
  
  
   Real dist[4];
   
   static const point3d N(1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0));
   static const point3d X(1.0,0.0,0.0);
    
   dist[0] = abs( Y.getZ() );
   dist[1] = abs( Y.getY() );
   dist[2] = abs( point3d::dot(Y-X, N) );
   dist[3] = abs( Y.getX() );
   
   Real minDist = std::min(std::min(dist[0],dist[1]), std::min(dist[2],dist[3]));
   
   dist[0] = Real(dist[0] == minDist);
   dist[1] = Real(dist[1] == minDist);
   dist[2] = Real(dist[2] == minDist);
   dist[3] = Real(dist[3] == minDist);
   
   Real tot = (dist[0] + dist[1] + dist[2] + dist[3]);
   
   val(1) = dist[0] / tot;
   val(2) = dist[1] / tot;
   val(3) = dist[2] / tot;
   val(4) = dist[3] / tot;
}

template<typename PMAPTYPE>
void
feFvFlux3d<PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  gradX(1) = 0.0; gradY(1) = 0.0; gradZ(1) = 0.0;
  gradX(2) = 0.0; gradY(2) = 0.0; gradZ(2) = 0.0;
  gradX(3) = 0.0; gradY(3) = 0.0; gradZ(3) = 0.0;
  gradX(4) = 0.0; gradY(4) = 0.0; gradZ(4) = 0.0;
}

template<typename PMAPTYPE>
void
feFvFlux3d<PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == numPoly);
  
  localEval(Y,val);
}

template<typename PMAPTYPE>
void
feFvFlux3d<PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == numPoly);
  assert(gradY.size() == numPoly);
  assert(gradZ.size() == numPoly);
  
  gradX(1) = 0.0; gradY(1) = 0.0; gradZ(1) = 0.0;
  gradX(2) = 0.0; gradY(2) = 0.0; gradZ(2) = 0.0;
  gradX(3) = 0.0; gradY(3) = 0.0; gradZ(3) = 0.0;
  gradX(4) = 0.0; gradY(4) = 0.0; gradZ(4) = 0.0;
}

#endif
