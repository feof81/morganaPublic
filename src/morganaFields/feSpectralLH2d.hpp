/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FESPECTRALLH2D_HPP
#define FESPECTRALLH2D_HPP

#include "glBase.h"
#include "polyDynamic.h"

#include "geoShapes.h"

#include "morganaFields.hpp"
#include "feSpectralLH2d_card.h"
#include "geoShapes.h"


/*! Spectral 2d finite elements */
template<typename PMAPTYPE>
class feSpectralLH2d
{
  /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearQuad                  GEOSHAPE;
    typedef feSpectralLH2d_card         FECARD;
    typedef elCard2d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feDynamicDofCard2d          DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feDynamic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_SK_2d;
    static const FEBaseLabel feBaseLabel = BL_SK_2d;
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt dim     = 2;
    static const bool isNodal = false;
    //@}
    
    /*! @name Links and flags */ //@{
  public:
    bool elCardLoaded;
    bool feCardLoaded;
    FECARD FeCard;
    ELCARD ElCard;
    //@}
    
    /*! @name Internal data */ //@{
  public:
    sVect<polyDynamic> basis;
    //@}
    
    /*! @name Constructors and setting functions */ //@{
  public:
    feSpectralLH2d();
    feSpectralLH2d(const FECARD & FECard, const ELCARD & ELCard);
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
feSpectralLH2d<PMAPTYPE>::
feSpectralLH2d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feSpectralLH2d<PMAPTYPE>::
feSpectralLH2d(const FECARD & FECard, const ELCARD & ELCard)
{
  //Copy the data
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
  
  //Polynomial generation
  UInt k=1;
  UInt numPoly = (FeCard.getRy()+1) * (FeCard.getRx()+1);
  
  glBase spectral;
  polyDynamicCard outPoly;
  basis.resize(numPoly);
  
  for(UInt iy=1; iy <= (FeCard.getRy()+1); ++iy)
  {
    for(UInt ix=1; ix <= (FeCard.getRx()+1); ++ix)
    {
      outPoly = spectral.getPolynomial(FeCard.getRx(), FeCard.getRy(), 0, ix, iy, 1);
      basis(k).setPolyDynamicCard(outPoly);	
      ++k;
    }
  }
}    
    
template<typename PMAPTYPE>    
void
feSpectralLH2d<PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  //Copy the data
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
  
  //Polynomial generation
  UInt k=1;
  UInt numPoly = (FeCard.getRy()+1) * (FeCard.getRx()+1);
  
  glBase spectral;
  polyDynamicCard outPoly;
  basis.resize(numPoly);
   
  for(UInt iy=1; iy <= (FeCard.getRy()+1); ++iy)
  {
    for(UInt ix=1; ix <= (FeCard.getRx()+1); ++ix)
    {
      outPoly = spectral.getPolynomial(FeCard.getRx(), FeCard.getRy(), 0, ix, iy, 1);
      basis(k).setPolyDynamicCard(outPoly);
      ++k;
    }
  }
}

template<typename PMAPTYPE> 
UInt
feSpectralLH2d<PMAPTYPE>::
getNumBasis() const
{
  return( (FeCard.getRy()+1) * (FeCard.getRx()+1) );
}

template<typename PMAPTYPE> 
sVect<typename feSpectralLH2d<PMAPTYPE>::DOFCARD>
feSpectralLH2d<PMAPTYPE>::
getDofCards() const
{
  UInt lev=1;
  sVect<DOFCARD> out(getNumBasis());
  
  for(UInt iy=1; iy <= (FeCard.getRy()+1); ++iy)
  {
    for(UInt ix=1; ix <= (FeCard.getRx()+1); ++ix)
    {
      out(lev) = DOFCARD(VOLUME,lev,1);
      lev++;
    }
  }
  
  return(out);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH2d<PMAPTYPE>::
getMaxLevVertex() const
{
  return(0);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH2d<PMAPTYPE>::
getMaxLevEdge() const
{
  return(0);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH2d<PMAPTYPE>::
getMaxLevFace() const
{
  return(0);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH2d<PMAPTYPE>::
getMaxLevVolume() const
{
  return(getNumBasis());
}

template<typename PMAPTYPE>
void
feSpectralLH2d<PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(  val.size() == getNumBasis());
  assert(basis.size() == getNumBasis());
  
  for(UInt i=1; i <= getNumBasis(); ++i)
  {
    val(i) = basis(i).evaluate(Y);
  }
}

template<typename PMAPTYPE>
void
feSpectralLH2d<PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == getNumBasis());
  assert(gradY.size() == getNumBasis());
  assert(gradZ.size() == getNumBasis());
  assert(basis.size() == getNumBasis());
  
  for(UInt i=1; i <= getNumBasis(); ++i)
  {
    gradX(i) = basis(i).evaluateGradientX(Y);
    gradY(i) = basis(i).evaluateGradientY(Y);
    gradZ(i) = basis(i).evaluateGradientZ(Y);
  }
}

template<typename PMAPTYPE>
void
feSpectralLH2d<PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(  val.size() == getNumBasis());
  assert(basis.size() == getNumBasis());
  
  localEval(Y,val);
}

template<typename PMAPTYPE>
void
feSpectralLH2d<PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == getNumBasis());
  assert(gradY.size() == getNumBasis());
  assert(gradZ.size() == getNumBasis());
  assert(basis.size() == getNumBasis());
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
  T.completeThirdColoumn();
  T.computeInverse();
  
  //Local grad
  sVect<BASETYPE> locX(getNumBasis()), locY(getNumBasis()), locZ(getNumBasis());
  localGrad(Y,locX,locY,locZ);
  
  //Compute gradient
  for(UInt i=1; i <= getNumBasis(); ++i)
  {
    gradX(i) = locX(i) * T(1,1) + locY(i) * T(2,1) + locZ(i) * T(3,1);
    gradY(i) = locX(i) * T(1,2) + locY(i) * T(2,2) + locZ(i) * T(3,2);
    gradZ(i) = locX(i) * T(1,3) + locY(i) * T(2,3) + locZ(i) * T(3,3);
  }
}

#endif
