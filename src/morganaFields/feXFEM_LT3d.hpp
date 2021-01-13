/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This file is part of Morgana.
Author: Andrea Villa, andrea.villa81@fastwebnet.it

Morgana is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Morgana is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Morgana. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef FEXFEM_LT3D_HPP
#define FEXFEM_LT3D_HPP

#include "feXFEM_LT3d_card.hpp"
#include "xfemSupportLT3d.hpp"


/*! The XFEM finite element */
template<Int R, typename PMAPTYPE>
class feXFEM_LT3d
{
    /*! @name Typedefs */ //@{
  public:
    typedef PMAPTYPE                    FETYPE_PMAPTYPE;
    typedef Real                        BASETYPE;
    typedef linearTetra                 GEOSHAPE;
    typedef feXFEM_LT3d_card            FECARD;
    typedef elCard3d<GEOSHAPE,PMAPTYPE> ELCARD;
    typedef feDynamicDofCard3d          DOFCARD;
    //@}
    
    /*! @name Static flags */ //@{
  public:
    static const FEClass     feClass     = feDynamic;
    static const FEGridType  feGridType  = primal;
    static const FEBaseType  feBaseType  = scalar;
    static const FELabel     feLabel     = FE_XF_3d;
    static const FEBaseLabel feBaseLabel = BL_XF_3d;
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
    feXFEM_LT3d();
    feXFEM_LT3d(const FECARD & FECard, const ELCARD & ELCard);
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
template<Int R, typename PMAPTYPE>
feXFEM_LT3d<R,PMAPTYPE>::
feXFEM_LT3d()
{
  elCardLoaded = false;
  feCardLoaded = false;
}

template<Int R, typename PMAPTYPE>
feXFEM_LT3d<R,PMAPTYPE>::
feXFEM_LT3d(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  FeCard = FECard;
  ElCard = ELCard;
}

template<Int R, typename PMAPTYPE>
void
feXFEM_LT3d<R,PMAPTYPE>::
setCards(const FECARD & FECard, const ELCARD & ELCard)
{
  elCardLoaded = true;
  feCardLoaded = true;
  
  ElCard = ELCard;
  FeCard = FECard;
}

template<Int R, typename PMAPTYPE>
UInt
feXFEM_LT3d<R,PMAPTYPE>::
getNumBasis() const
{
  assert( (FeCard.getNumEnrichedBasis() <= fePr3d<R,PMAPTYPE>::numBasis) );
  return( FeCard.getIsActive() * (fePr3d<R,PMAPTYPE>::numBasis + FeCard.getNumEnrichedBasis()) );
}

template<Int R, typename PMAPTYPE>
sVect<typename feXFEM_LT3d<R,PMAPTYPE>::DOFCARD>
feXFEM_LT3d<R,PMAPTYPE>::
getDofCards() const
{
  typedef typename fePr3d<R,PMAPTYPE>::DOFCARD DOFCARD_STATIC;
  sVect<DOFCARD> out(getNumBasis());
  DOFCARD_STATIC staticCard;
  
  //Static cards
  for(UInt i=1; i <= fePr3d<R,PMAPTYPE>::numBasis; ++i)
  {
    staticCard = fePr3d<R,PMAPTYPE>::getDofCard(i);
    out(i) = DOFCARD(staticCard.getGeoType(),
                     staticCard.getLevel(),
		     staticCard.getLocalId() );
  }
  
  //Dynamic cards
  for(UInt i=1; i <= FeCard.getNumEnrichedBasis(); ++i)
  {
    assert( FeCard.getEnrichedDof(i) >= 1 );
    assert( (FeCard.getEnrichedDof(i) <= fePr3d<R,PMAPTYPE>::numBasis) );
    
    out(fePr3d<R,PMAPTYPE>::numBasis + i) = out(FeCard.getEnrichedDof(i));
    out(fePr3d<R,PMAPTYPE>::numBasis + i).setLevel(2);
  }
  
  return(out);
}

template<Int R, typename PMAPTYPE>
UInt
feXFEM_LT3d<R,PMAPTYPE>::
getMaxLevVertex() const
{
  
  return( (UInt(FeCard.getNumEnrichedBasis() == 0) * fePr3d<R,PMAPTYPE >::dofPerVertex) +
          (UInt(FeCard.getNumEnrichedBasis() != 0) * fePr3d<R,PMAPTYPE >::dofPerVertex * 2) );
}

template<Int R, typename PMAPTYPE>
UInt
feXFEM_LT3d<R,PMAPTYPE>::
getMaxLevEdge() const
{
  return( (UInt(FeCard.getNumEnrichedBasis() == 0) * fePr3d<R,PMAPTYPE >::dofPerEdge) +
          (UInt(FeCard.getNumEnrichedBasis() != 0) * fePr3d<R,PMAPTYPE >::dofPerEdge * 2) );
}

template<Int R, typename PMAPTYPE>
UInt
feXFEM_LT3d<R,PMAPTYPE>::
getMaxLevFace() const
{
  return( (UInt(FeCard.getNumEnrichedBasis() == 0) * fePr3d<R,PMAPTYPE >::dofPerFace) +
          (UInt(FeCard.getNumEnrichedBasis() != 0) * fePr3d<R,PMAPTYPE >::dofPerFace * 2) );
}

template<Int R, typename PMAPTYPE>
UInt
feXFEM_LT3d<R,PMAPTYPE>::
getMaxLevVolume() const
{
  return( (UInt(FeCard.getNumEnrichedBasis() == 0) * fePr3d<R,PMAPTYPE >::dofPerVolume) +
          (UInt(FeCard.getNumEnrichedBasis() != 0) * fePr3d<R,PMAPTYPE >::dofPerVolume * 2) );
}

template<Int R, typename PMAPTYPE>
void
feXFEM_LT3d<R,PMAPTYPE>::
localEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == getNumBasis());
  
  //Static polynomial evaluation  
  feStaticRealEval< fePr3d<R,PMAPTYPE>::feLabel, fePr3d<R,PMAPTYPE >::numPoly >::eval(val,Y);
  
  //Enrichment
  Real x = Y.getX();
  Real y = Y.getY();
  Real z = Y.getZ();
  
  Real phi = (1.0 - x - y - z) * FeCard.getPhi(1) +
                            x  * FeCard.getPhi(2) +
		            y  * FeCard.getPhi(3) +
		            z  * FeCard.getPhi(4);
  
  for(UInt i=1; i <= FeCard.getNumEnrichedBasis(); ++i)
  {
    val(fePr3d<R,PMAPTYPE>::numBasis + i) = val(FeCard.getEnrichedDof(i)) * abs(phi - 0.5);
  }
}

template<Int R, typename PMAPTYPE>
void
feXFEM_LT3d<R,PMAPTYPE>::
localGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == getNumBasis());
  assert(gradY.size() == getNumBasis());
  assert(gradZ.size() == getNumBasis());
  
  //Static data
  static const Real Hx[4] = {-1.0, 1.0, 0.0, 0.0};
  static const Real Hy[4] = {-1.0, 0.0, 1.0, 0.0};
  static const Real Hz[4] = {-1.0, 0.0, 0.0, 1.0};
  
  //Static polynomial derivarive evaluation
  feStaticRealDx< fePr3d<R,PMAPTYPE>::feLabel, fePr3d<R,PMAPTYPE >::numPoly >::eval(gradX,Y);
  feStaticRealDy< fePr3d<R,PMAPTYPE>::feLabel, fePr3d<R,PMAPTYPE >::numPoly >::eval(gradY,Y);
  feStaticRealDz< fePr3d<R,PMAPTYPE>::feLabel, fePr3d<R,PMAPTYPE >::numPoly >::eval(gradZ,Y);
  
  //Static polynomial evaluation  
   sVect<BASETYPE> val(fePr3d<R,PMAPTYPE>::numBasis);
  feStaticRealEval< fePr3d<R,PMAPTYPE>::feLabel, fePr3d<R,PMAPTYPE >::numPoly >::eval(val,Y);
  
  //Enrichment
  Real x = Y.getX();
  Real y = Y.getY();
  Real z = Y.getZ();
  
  Real phi = (1.0 - x - y - z) * FeCard.getPhi(1) +
                            x  * FeCard.getPhi(2) +
		            y  * FeCard.getPhi(3) +
		            z  * FeCard.getPhi(4);
			    
  Real phi_x = Hx[0] *  FeCard.getPhi(1) + Hx[1] *  FeCard.getPhi(2) + Hx[2] *  FeCard.getPhi(3) + Hx[3] *  FeCard.getPhi(4);
  Real phi_y = Hy[0] *  FeCard.getPhi(1) + Hy[1] *  FeCard.getPhi(2) + Hy[2] *  FeCard.getPhi(3) + Hy[3] *  FeCard.getPhi(4);
  Real phi_z = Hz[0] *  FeCard.getPhi(1) + Hz[1] *  FeCard.getPhi(2) + Hz[2] *  FeCard.getPhi(3) + Hz[3] *  FeCard.getPhi(4);
			    
  for(UInt i=1; i <= FeCard.getNumEnrichedBasis(); ++i)
  {
    gradX(fePr3d<R,PMAPTYPE>::numBasis + i) = gradX(FeCard.getEnrichedDof(i)) * abs(phi - 0.5)  +  val(i) * ((phi >= 0.5) - (phi  < 0.5)) * phi_x;
    gradY(fePr3d<R,PMAPTYPE>::numBasis + i) = gradY(FeCard.getEnrichedDof(i)) * abs(phi - 0.5)  +  val(i) * ((phi >= 0.5) - (phi  < 0.5)) * phi_y;
    gradZ(fePr3d<R,PMAPTYPE>::numBasis + i) = gradZ(FeCard.getEnrichedDof(i)) * abs(phi - 0.5)  +  val(i) * ((phi >= 0.5) - (phi  < 0.5)) * phi_z;
  }
}

template<Int R, typename PMAPTYPE>
void
feXFEM_LT3d<R,PMAPTYPE>::
globalEval(const point3d & Y, sVect<BASETYPE> & val) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(val.size() == getNumBasis());
  
  localEval(Y,val);
}

template<Int R, typename PMAPTYPE>
void
feXFEM_LT3d<R,PMAPTYPE>::
globalGrad(const point3d & Y, sVect<BASETYPE> & gradX, sVect<BASETYPE> & gradY, sVect<BASETYPE> & gradZ) const
{
  assert(elCardLoaded);
  assert(feCardLoaded);
  assert(gradX.size() == getNumBasis());
  assert(gradY.size() == getNumBasis());
  assert(gradZ.size() == getNumBasis());
  
  //Compute tensor
  geoMapInterface<GEOSHAPE> geoInterface;
  tensor3d T = geoInterface.getGradient(ElCard.getNodes(),Y);
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
