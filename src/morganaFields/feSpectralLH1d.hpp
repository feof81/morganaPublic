/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FESPECTRALLH1D_HPP
#define FESPECTRALLH1D_HPP

#include "glBase.h"
#include "polyDynamic.h"

#include "geoShapes.h"

#include "morganaFields.hpp"
#include "feSpectralLH1d_card.h"
#include "geoShapes.h"


/*! Spectral 1d finite elements */
template<typename PMAPTYPE>
class feSpectralLH1d
{
  /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearLine                  GEOSHAPE;
    typedef feSpectralLH1d_card         FECARD;
    typedef elCard1d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feDynamicDofCard1d          DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public: 
    static const FEClass     feClass     = feDynamic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_SK_1d;
    static const FEBaseLabel feBaseLabel = BL_SK_1d;    
    //@}
    
    /*! @name Static numbers */ //@{
  public:
    static const UInt dim     = 1;
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
    feSpectralLH1d();
    feSpectralLH1d(const FECARD & FECard, const ELCARD & ELCard);
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
feSpectralLH1d<PMAPTYPE>::
feSpectralLH1d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<typename PMAPTYPE>
feSpectralLH1d<PMAPTYPE>::
feSpectralLH1d(const FECARD & FECard, const ELCARD & ELCard)
{
  //Copy the data
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
  
  //Polynomial generation
  UInt k=1;
  UInt numPoly = FeCard.getRx() + 1;
  
  glBase spectral;
  polyDynamicCard outPoly;
  basis.resize(numPoly);
  
  for(UInt ix=1; ix <= (FeCard.getRx() + 1); ++ix)
  {
    outPoly = spectral.getPolynomial(FeCard.getRx(), 0, 0, ix, 1, 1);
    basis(k).setPolyDynamicCard(outPoly);	
    ++k;
  }
}    
    
template<typename PMAPTYPE>    
void
feSpectralLH1d<PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  //Copy the data
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
  
  //Polynomial generation
  UInt k=1;
  UInt numPoly = FeCard.getRx() + 1;
  
  glBase spectral;
  polyDynamicCard outPoly;
  basis.resize(numPoly);
   
  for(UInt ix=1; ix <= (FeCard.getRx() + 1); ++ix)
  {
    outPoly = spectral.getPolynomial(FeCard.getRx(), 0, 0, ix, 1, 1);
    basis(k).setPolyDynamicCard(outPoly);
    ++k;
  }
}

template<typename PMAPTYPE> 
UInt
feSpectralLH1d<PMAPTYPE>::
getNumBasis() const
{
  return( FeCard.getRx() + 1 );
}

template<typename PMAPTYPE> 
sVect<typename feSpectralLH1d<PMAPTYPE>::DOFCARD>
feSpectralLH1d<PMAPTYPE>::
getDofCards() const
{
  UInt lev=1;
  sVect<DOFCARD> out(getNumBasis());
  
  for(UInt ix=1; ix <= (FeCard.getRx() + 1); ++ix)
  {
    out(lev) = DOFCARD(VOLUME,lev,1);
    lev++;
  }
  
  return(out);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH1d<PMAPTYPE>::
getMaxLevVertex() const
{
  return(0);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH1d<PMAPTYPE>::
getMaxLevEdge() const
{
  return(0);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH1d<PMAPTYPE>::
getMaxLevFace() const
{
  return(0);
}

template<typename PMAPTYPE> 
UInt
feSpectralLH1d<PMAPTYPE>::
getMaxLevVolume() const
{
  return(getNumBasis());
}

template<typename PMAPTYPE>
void
feSpectralLH1d<PMAPTYPE>::
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
feSpectralLH1d<PMAPTYPE>::
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
feSpectralLH1d<PMAPTYPE>::
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
feSpectralLH1d<PMAPTYPE>::
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
  T.completeSecondColoumn();
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
