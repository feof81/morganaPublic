/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEHYB3D_HPP
#define FEHYB3D_HPP

#include "feHYB3d_card.h"
#include "morganaFields.hpp"
#include "morganaGeometry.hpp"
#include "morganaTypes.hpp"
#include "feDynamicDofCard3d.h"
#include "elCard3d.hpp"


/*! The hybrid interconnecting finite element for primal methods */
template<typename PMAPTYPE>
class feHYB3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feHYB3d_card                FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feDynamicDofCard3d          DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feDynamic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_HYB_3d;
    static const FEBaseLabel feBaseLabel = BL_HYB_3d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt dim     = 3;
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
    feHYB3d();
    feHYB3d(const FECARD & FECard, const ELCARD & ELCard);
    void setCards(const FECARD & FECard, const ELCARD & ELCard);
    //@}
    
    /*! @name Info functions */ //@{
  public:
    UInt           getNumBasis() const;
    sVect<DOFCARD> getDofCards() const;
    UInt           getMaxLevVertex() const;
    UInt           getMaxLevEdge() const;
    UInt           getMaxLevFace() const;
    UInt           getMaxLevVolume() const;
    //@}
    
    /*! @name Eval functions */ //@{
  public:
    void localEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    void globalEval(const point3d & Y, sVect<BASETYPE> & val) const;
    void globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const;
    //@}
};



//_________________________________________________________________________________________________
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------
template<typename PMAPTYPE>
feHYB3d<PMAPTYPE>::
feHYB3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feHYB3d<PMAPTYPE>::
feHYB3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<typename PMAPTYPE>
void
feHYB3d<PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<typename PMAPTYPE>
UInt
feHYB3d<PMAPTYPE>::
getNumBasis() const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  
  return( UInt(FeCard.getActiveFaces()(1)) +
          UInt(FeCard.getActiveFaces()(2)) +
          UInt(FeCard.getActiveFaces()(3)) + 
          UInt(FeCard.getActiveFaces()(4)) );
}

template<typename PMAPTYPE>
sVect<typename feHYB3d<PMAPTYPE>::DOFCARD>
feHYB3d<PMAPTYPE>::
getDofCards() const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  
  sVect<DOFCARD> out;
  
  for(UInt i=1; i <= 4; ++i)
  {
    if(FeCard.getActiveFaces()(i))
    { out.push_back(DOFCARD(FACE,1,i)); }
  }
  
  return(out);
}

template<typename PMAPTYPE>
UInt
feHYB3d<PMAPTYPE>::
getMaxLevVertex() const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  
  return(0);
}

template<typename PMAPTYPE>
UInt
feHYB3d<PMAPTYPE>::
getMaxLevEdge() const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  
  return(0);
}

template<typename PMAPTYPE>
UInt
feHYB3d<PMAPTYPE>::
getMaxLevFace() const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  
  return(1);
}

template<typename PMAPTYPE>
UInt
feHYB3d<PMAPTYPE>::
getMaxLevVolume() const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  
  return(0);
}

template<typename PMAPTYPE>
void
feHYB3d<PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(Y.getX() == Y.getX());
  assert(val.size() == getNumBasis());
  
  //Distance info
  sVect<bool> act(4);
  act(1) = (Y.getZ() <= geoToll)    &&  (Y.getZ() >= (-geoToll))      &&
           (Y.getY() >= (-geoToll)) &&  (Y.getY() <= (1.0 + geoToll)) &&
           (Y.getX() >= (-geoToll)) &&  (Y.getX() <= (1.0 + geoToll)) &&
           ( (Y.getX() + Y.getY()) <= (1.0 + geoToll) );
	   
  act(2) = (Y.getY() <= geoToll)    &&  (Y.getY() >= (-geoToll))      &&
           (Y.getZ() >= (-geoToll)) &&  (Y.getZ() <= (1.0 + geoToll)) &&
           (Y.getX() >= (-geoToll)) &&  (Y.getX() <= (1.0 + geoToll)) &&
           ( (Y.getX() + Y.getZ()) <= (1.0 + geoToll) );
  
  act(3) = (Y.getX() >= (-geoToll)) &&  (Y.getX() <= (1.0 + geoToll)) &&
           (Y.getY() >= (-geoToll)) &&  (Y.getY() <= (1.0 + geoToll)) &&
           (Y.getZ() >= (-geoToll)) &&  (Y.getZ() <= (1.0 + geoToll)) &&
           ( (Y.getX() + Y.getY() + Y.getZ()) <= (1.0 + geoToll) )    &&
           ( (Y.getX() + Y.getY() + Y.getZ()) >= (1.0 - geoToll) );
   
  act(4) = (Y.getX() <= geoToll)    &&  (Y.getX() >= (-geoToll))      &&
           (Y.getZ() >= (-geoToll)) &&  (Y.getZ() <= (1.0 + geoToll)) &&
           (Y.getY() >= (-geoToll)) &&  (Y.getY() <= (1.0 + geoToll)) &&
           ( (Y.getY() + Y.getZ()) <= (1.0 + geoToll) );
  
  UInt k=1;
  
  if(FeCard.getActiveFaces(1))
  {val(k) = Real(act(1)) * (Real(FeCard.getFaceOrientation(1) == true) - Real(FeCard.getFaceOrientation(1) == false)); ++k;}
  
  if(FeCard.getActiveFaces(2))
  {val(k) = Real(act(2)) * (Real(FeCard.getFaceOrientation(2) == true) - Real(FeCard.getFaceOrientation(2) == false)); ++k;}
  
  if(FeCard.getActiveFaces(3))
  {val(k) = Real(act(3)) * (Real(FeCard.getFaceOrientation(3) == true) - Real(FeCard.getFaceOrientation(3) == false)); ++k;}
  
  if(FeCard.getActiveFaces(4))
  {val(k) = Real(act(4)) * (Real(FeCard.getFaceOrientation(4) == true) - Real(FeCard.getFaceOrientation(4) == false)); ++k;}
	   
  assert((k-1) == getNumBasis());
  assert( ( UInt(act(1)) +  UInt(act(2)) +  UInt(act(3)) +  UInt(act(4)) ) == 1);
}

template<typename PMAPTYPE>
void
feHYB3d<PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(Y.getX() == Y.getX());
  assert(gradX.size() == getNumBasis());
  assert(gradY.size() == getNumBasis());
  assert(gradZ.size() == getNumBasis());
  
  for(UInt i=1; i <= getNumBasis(); ++i)
  {
    gradX(i) = 0;
    gradY(i) = 0;
    gradZ(i) = 0;
  }
}

template<typename PMAPTYPE>
void
feHYB3d<PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(Y.getX() == Y.getX());
  assert(val.size() == getNumBasis());
  
  localEval(Y,val);
}

template<typename PMAPTYPE>
void
feHYB3d<PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(Y.getX() == Y.getX());
  assert(gradX.size() == getNumBasis());
  assert(gradY.size() == getNumBasis());
  assert(gradZ.size() == getNumBasis());
  
  for(UInt i=1; i <= getNumBasis(); ++i)
  {
    gradX(i) = 0;
    gradY(i) = 0;
    gradZ(i) = 0;
  }
}

#endif
